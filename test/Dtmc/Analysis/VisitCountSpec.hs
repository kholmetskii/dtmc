{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.VisitCountSpec (
    spec,
) where

import Data.Finite (
    Finite,
    finites,
 )
import Dtmc.Analysis.Classification (
    accessible,
    recurrentState,
 )
import Dtmc.Analysis.Event (
    DiscreteEvent (..),
 )
import Dtmc.Analysis.FiniteTime qualified as FT
import Dtmc.Analysis.HittingTime qualified as Hit
import Dtmc.Analysis.ReturnTime qualified as Return
import Dtmc.Analysis.VisitCount (
    Expectation (..),
 )
import Dtmc.Analysis.VisitCount qualified as Visit
import Dtmc.Distribution (
    distributionWeights,
 )
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.TestSupport (
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.Transition.Kernel qualified as Kernel
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    mkTransitionMatrix,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
    shouldSatisfy,
 )
import Test.Hspec.QuickCheck (
    prop,
 )
import Test.QuickCheck (
    choose,
    conjoin,
    counterexample,
    forAll,
    property,
    (===),
 )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

twoCycle :: TransitionMatrix (Finite 2)
twoCycle =
    checked $
        mkTransitionMatrix
            ( S.matrix
                [ 0
                , 1
                , 1
                , 0
                ] ::
                S.Sq 2
            )

-- Target 0 returns with probability 1/4. State 1 first reaches it with
-- probability 1/2, while absorbing state 2 cannot reach it.
transientVisitChain :: TransitionMatrix (Finite 3)
transientVisitChain =
    checked $
        mkTransitionMatrix
            ( S.matrix
                [ 1 / 4
                , 0
                , 3 / 4
                , 1 / 2
                , 0
                , 1 / 2
                , 0
                , 0
                , 1
                ] ::
                S.Sq 3
            )

-- States 0 and 1 may enter absorbing target 2; absorbing state 3 cannot.
recurrentVisitChain :: TransitionMatrix (Finite 4)
recurrentVisitChain =
    checked $
        mkTransitionMatrix
            ( S.matrix
                [ 0
                , 1 / 2
                , 1 / 2
                , 0
                , 1 / 2
                , 0
                , 0
                , 1 / 2
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 1
                ] ::
                S.Sq 4
            )

mixedInitial :: DistributionMap.DistributionMap (Finite 2)
mixedInitial =
    checked (DistributionMap.mkDistributionMap [(0, 0.25), (1, 0.75)])

asKernel ::
    TransitionMatrix (Finite 2) ->
    Kernel.TransitionKernel (Finite 2)
asKernel matrix =
    Kernel.transitionKernel $ \source ->
        checked $
            DistributionMap.mkDistributionMap
                [ (destination, FT.stepProbability matrix source destination)
                | destination <- finites
                ]

simpleRandomWalk :: Kernel.TransitionKernel Integer
simpleRandomWalk =
    Kernel.transitionKernel $ \state ->
        checked $
            DistributionMap.mkDistributionMap
                [(state - 1, 0.5), (state + 1, 0.5)]

closeTo :: Double -> Double -> Bool
closeTo expected actual = abs (actual - expected) <= testTolerance

entries :: (KnownNat n) => S.R n -> [Double]
entries = LA.toList . S.extract

expectationCloseTo :: Double -> Expectation -> Bool
expectationCloseTo expected (FiniteExpectation actual) = closeTo expected actual
expectationCloseTo _ InfiniteExpectation = False

spec :: Spec
spec = do
    describe "totalProbabilityByState" $ do
        it "matches the geometric law for a transient target" $ do
            let probabilities count =
                    entries (checked ((Visit.totalProbabilityByState . EqualTo) count transientVisitChain 0))
            sequence_
                [ actual `shouldSatisfy` closeTo expected
                | (actual, expected) <-
                    zip (probabilities 0) [0, 1 / 2, 1]
                ]
            sequence_
                [ actual `shouldSatisfy` closeTo expected
                | (actual, expected) <-
                    zip (probabilities 1) [3 / 4, 3 / 8, 0]
                ]
            sequence_
                [ actual `shouldSatisfy` closeTo expected
                | (actual, expected) <-
                    zip (probabilities 3) [3 / 64, 3 / 128, 0]
                ]
            entries (checked (Visit.infiniteProbabilityByState transientVisitChain 0))
                `shouldBe` [0, 0, 0]

        it "puts all positive recurrent-target mass at infinity" $ do
            entries
                (checked ((Visit.totalProbabilityByState . EqualTo) 1 recurrentVisitChain 2))
                `shouldBe` [0, 0, 0, 0]
            sequence_
                [ actual `shouldSatisfy` closeTo expected
                | (actual, expected) <-
                    zip
                        (entries (checked (Visit.infiniteProbabilityByState recurrentVisitChain 2)))
                        [2 / 3, 1 / 3, 1, 0]
                ]
            sequence_
                [ actual `shouldSatisfy` closeTo expected
                | (actual, expected) <-
                    zip
                        (entries (checked ((Visit.totalProbabilityByState . EqualTo) 0 recurrentVisitChain 2)))
                        [1 / 3, 2 / 3, 0, 1]
                ]

        it "counts the target at time zero" $ do
            checked ((Visit.totalProbability . EqualTo) 0 transientVisitChain 0 0)
                `shouldBe` 0
            checked ((Visit.totalProbability . EqualTo) 1 transientVisitChain 0 0)
                `shouldSatisfy` closeTo (3 / 4)

        prop "scalar queries look up the all-state result (random @3)" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left err -> counterexample (show err) False
                    Right matrix ->
                        conjoin
                            [ conjoin
                                [ case (Visit.totalProbabilityByState . EqualTo) count matrix 0 of
                                    Left err -> counterexample (show err) False
                                    Right probabilities ->
                                        conjoin
                                            [ (Visit.totalProbability . EqualTo) count matrix 0 initial
                                                === Right probability
                                            | (initial, probability) <-
                                                zip (finites :: [Finite 3]) (entries probabilities)
                                            ]
                                | count <- [0, 1, 3]
                                ]
                            , case Visit.infiniteProbabilityByState matrix 0 of
                                Left err -> counterexample (show err) False
                                Right probabilities ->
                                    conjoin
                                        [ Visit.infiniteProbability matrix 0 initial
                                            === Right probability
                                        | (initial, probability) <-
                                            zip (finites :: [Finite 3]) (entries probabilities)
                                        ]
                            ]

    describe "totalExpectationByState" $ do
        it "matches h / (1 - f) for a transient target" $ do
            let expectations = checked (Visit.totalExpectationByState transientVisitChain 0)
            sequence_
                [ expectation `shouldSatisfy` expectationCloseTo expected
                | (expectation, expected) <- zip expectations [4 / 3, 2 / 3, 0]
                ]
            checked (Visit.totalExpectation transientVisitChain 0 1)
                `shouldSatisfy` expectationCloseTo (2 / 3)

        it "is infinite exactly where a recurrent target is reachable" $ do
            checked (Visit.totalExpectationByState recurrentVisitChain 2)
                `shouldBe` [InfiniteExpectation, InfiniteExpectation, InfiniteExpectation, FiniteExpectation 0]
            checked (Visit.totalExpectation recurrentVisitChain 2 3)
                `shouldBe` FiniteExpectation 0

        prop "agrees with hitting, return, recurrence, and reachability (random @3)" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left err -> counterexample (show err) False
                    Right matrix ->
                        case do
                            hits <- Hit.eventualProbabilityByState matrix [0]
                            returning <- Return.eventualProbability matrix 0
                            zeroVisits <- (Visit.totalProbabilityByState . EqualTo) 0 matrix 0
                            oneVisit <- (Visit.totalProbabilityByState . EqualTo) 1 matrix 0
                            twoVisits <- (Visit.totalProbabilityByState . EqualTo) 2 matrix 0
                            infiniteVisits <- Visit.infiniteProbabilityByState matrix 0
                            expectations <- Visit.totalExpectationByState matrix 0
                            pure
                                ( hits
                                , returning
                                , zeroVisits
                                , oneVisit
                                , twoVisits
                                , infiniteVisits
                                , expectations
                                ) of
                            Left err -> counterexample (show err) False
                            Right
                                ( hits
                                    , returning
                                    , zeroVisits
                                    , oneVisit
                                    , twoVisits
                                    , infiniteVisits
                                    , expectations
                                    ) ->
                                    let hitValues = entries hits
                                        zeroValues = entries zeroVisits
                                        oneValues = entries oneVisit
                                        twoValues = entries twoVisits
                                        infiniteValues = entries infiniteVisits
                                        states = finites :: [Finite 3]
                                        structuralExpectations =
                                            [ if accessible matrix initial 0
                                                then InfiniteExpectation
                                                else FiniteExpectation 0
                                            | initial <- states
                                            ]
                                     in conjoin
                                            [ conjoin
                                                [ property (closeTo (1 - hit) zero)
                                                | (hit, zero) <- zip hitValues zeroValues
                                                ]
                                            , if recurrentState matrix 0
                                                then
                                                    conjoin
                                                        [ oneValues === [0, 0, 0]
                                                        , twoValues === [0, 0, 0]
                                                        , infiniteValues === hitValues
                                                        , expectations === structuralExpectations
                                                        ]
                                                else
                                                    conjoin
                                                        [ infiniteValues === [0, 0, 0]
                                                        , conjoin
                                                            [ property (closeTo (hit * (1 - returning)) one)
                                                            | (hit, one) <- zip hitValues oneValues
                                                            ]
                                                        , conjoin
                                                            [ property (closeTo (one * returning) two)
                                                            | (one, two) <- zip oneValues twoValues
                                                            ]
                                                        , conjoin
                                                            [ case expectation of
                                                                FiniteExpectation value ->
                                                                    property
                                                                        (closeTo (hit / (1 - returning)) value)
                                                                InfiniteExpectation -> property False
                                                            | (hit, expectation) <- zip hitValues expectations
                                                            ]
                                                        ]
                                            ]

    describe "boundedLaw" $ do
        it "is a point mass at zero for bound zero" $
            distributionWeights
                ( Visit.boundedLaw
                    0
                    (DistributionMap.pointMass (0 :: Finite 2))
                    twoCycle
                    (== 0)
                )
                `shouldBe` [(0, 1)]

        it "counts the initial state at a positive bound" $
            distributionWeights
                (Visit.boundedLaw 1 mixedInitial twoCycle (== 0))
                `shouldBe` [(0, 0.75), (1, 0.25)]

        it "counts deterministic visits at times zero through bound minus one" $ do
            let initial = DistributionMap.pointMass (0 :: Finite 2)
            distributionWeights (Visit.boundedLaw 1 initial twoCycle (== 0))
                `shouldBe` [(1, 1)]
            distributionWeights (Visit.boundedLaw 2 initial twoCycle (== 0))
                `shouldBe` [(1, 1)]
            distributionWeights (Visit.boundedLaw 3 initial twoCycle (== 0))
                `shouldBe` [(2, 1)]

        it "computes an exact law on an infinite random walk" $
            distributionWeights
                ( Visit.boundedLaw
                    3
                    (DistributionMap.pointMass (0 :: Integer))
                    simpleRandomWalk
                    (== 0)
                )
                `shouldBe` [(1, 0.5), (2, 0.5)]

        prop "has total mass one and no count above the bound (random @3)" $
            forAll (choose (0, 5 :: Int)) $ \rawBound ->
                forAll (genTransitionMatrix @3) $ \rawMatrix ->
                    case mkTransitionMatrix rawMatrix of
                        Left err -> counterexample (show err) False
                        Right matrix ->
                            let bound = fromIntegral rawBound
                                law =
                                    Visit.boundedLaw
                                        bound
                                        (DistributionMap.pointMass (0 :: Finite 3))
                                        matrix
                                        (== 0)
                                weights = distributionWeights law
                             in conjoin
                                    [ property (closeTo 1 (sum (map snd weights)))
                                    , property (all ((<= bound) . fst) weights)
                                    ]

    describe "boundedProbability" $ do
        it "looks up one coordinate of the count distribution" $ do
            let initial = DistributionMap.pointMass (0 :: Integer)
            Visit.boundedProbability 3 (EqualTo 1) initial simpleRandomWalk (== 0)
                `shouldSatisfy` closeTo 0.5
            Visit.boundedProbability 3 (EqualTo 2) initial simpleRandomWalk (== 0)
                `shouldSatisfy` closeTo 0.5
            Visit.boundedProbability 3 (EqualTo 3) initial simpleRandomWalk (== 0)
                `shouldBe` 0

        it "agrees for a matrix and its equivalent kernel" $ do
            let initial = DistributionMap.pointMass (0 :: Finite 2)
                kernel = asKernel twoCycle
            sequence_
                [ Visit.boundedProbability bound (EqualTo count) initial twoCycle (== 0)
                    `shouldBe` Visit.boundedProbability bound (EqualTo count) initial kernel (== 0)
                | bound <- [0 .. 5]
                , count <- [0 .. bound]
                ]

    describe "boundedExpectation" $ do
        it "is the expectation of the random-walk count law" $
            Visit.boundedExpectation
                3
                (DistributionMap.pointMass (0 :: Integer))
                simpleRandomWalk
                (== 0)
                `shouldSatisfy` closeTo 1.5

        prop "equals the sum of finite-time visit probabilities (random @3)" $
            forAll (choose (0, 5 :: Int)) $ \rawBound ->
                forAll (genTransitionMatrix @3) $ \rawMatrix ->
                    case mkTransitionMatrix rawMatrix of
                        Left err -> counterexample (show err) False
                        Right matrix ->
                            let bound = fromIntegral rawBound
                                initial = DistributionMap.pointMass (0 :: Finite 3)
                                expectation =
                                    Visit.boundedExpectation bound initial matrix (== 0)
                                marginalSum =
                                    sum
                                        [ FT.stateProbability (fromIntegral time) initial matrix 0
                                        | time <- [0 .. rawBound - 1]
                                        ]
                             in property (closeTo marginalSum expectation)

    describe "occupationMatrix" $ do
        it "matches the closed form of the transient block" $
            Visit.occupationMatrix transientVisitChain
                `shouldSatisfy` matchesOccupation
                    [ [Just (4 / 3), Just 0, Nothing]
                    , [Just (2 / 3), Just 1, Nothing]
                    , [Just 0, Just 0, Nothing]
                    ]

        it "is infinite everywhere when no state is transient" $
            Visit.occupationMatrix twoCycle
                `shouldSatisfy` matchesOccupation
                    [ [Nothing, Nothing]
                    , [Nothing, Nothing]
                    ]

        prop "agrees with totalExpectation entry by entry" $
            forAll (genTransitionMatrix @3) $ \m ->
                case mkTransitionMatrix m of
                    Left err -> counterexample (show err) False
                    Right p ->
                        case Visit.occupationMatrix p of
                            -- A refused solve is a documented outcome.
                            Left _ -> property True
                            Right rows ->
                                conjoin
                                    [ counterexample
                                        (show (i, j, entry))
                                        (agreesWithSingle entry (Visit.totalExpectation p j i))
                                    | (i, row) <- zip (finites :: [Finite 3]) rows
                                    , (j, entry) <- zip (finites :: [Finite 3]) row
                                    ]

-- Nothing stands for InfiniteExpectation; Just v for a finite entry near v.
matchesOccupation ::
    [[Maybe Double]] ->
    Either error [[Expectation]] ->
    Bool
matchesOccupation _ (Left _) = False
matchesOccupation expected (Right actual) =
    length expected == length actual
        && and (zipWith matchesRow expected actual)
  where
    matchesRow e a = length e == length a && and (zipWith matchesEntry e a)
    matchesEntry (Just v) x = expectationCloseTo v x
    matchesEntry Nothing InfiniteExpectation = True
    matchesEntry Nothing _ = False

agreesWithSingle :: Expectation -> Either error Expectation -> Bool
agreesWithSingle _ (Left _) = True
agreesWithSingle InfiniteExpectation (Right InfiniteExpectation) = True
agreesWithSingle (FiniteExpectation x) (Right (FiniteExpectation y)) = closeTo x y
agreesWithSingle _ _ = False
