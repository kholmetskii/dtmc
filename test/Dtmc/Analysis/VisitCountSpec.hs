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
import Dtmc.Analysis.FiniteTime (
    probabilityAtTime,
    transitionProbability,
 )
import Dtmc.Analysis.HittingTime (
    hittingProbabilities,
 )
import Dtmc.Analysis.ReturnTime (
    returnProbability,
 )
import Dtmc.Analysis.VisitCount (
    MeanCount (..),
    VisitCountOutcome (..),
    visitCountDistributionBefore,
    visitCountExpectation,
    visitCountExpectationBefore,
    visitCountExpectations,
    visitCountProbabilities,
    visitCountProbability,
    visitCountProbabilityBefore,
 )
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
                [ (destination, transitionProbability matrix source destination)
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

meanCloseTo :: Double -> MeanCount -> Bool
meanCloseTo expected (FiniteMeanCount actual) = closeTo expected actual
meanCloseTo _ InfiniteMeanCount = False

spec :: Spec
spec = do
    describe "visitCountProbabilities" $ do
        it "matches the geometric law for a transient target" $ do
            let probabilities outcome =
                    entries (checked (visitCountProbabilities transientVisitChain 0 outcome))
            sequence_
                [ actual `shouldSatisfy` closeTo expected
                | (actual, expected) <-
                    zip (probabilities (FiniteVisits 0)) [0, 1 / 2, 1]
                ]
            sequence_
                [ actual `shouldSatisfy` closeTo expected
                | (actual, expected) <-
                    zip (probabilities (FiniteVisits 1)) [3 / 4, 3 / 8, 0]
                ]
            sequence_
                [ actual `shouldSatisfy` closeTo expected
                | (actual, expected) <-
                    zip (probabilities (FiniteVisits 3)) [3 / 64, 3 / 128, 0]
                ]
            probabilities InfiniteVisits `shouldBe` [0, 0, 0]

        it "puts all positive recurrent-target mass at infinity" $ do
            entries
                (checked (visitCountProbabilities recurrentVisitChain 2 (FiniteVisits 1)))
                `shouldBe` [0, 0, 0, 0]
            sequence_
                [ actual `shouldSatisfy` closeTo expected
                | (actual, expected) <-
                    zip
                        (entries (checked (visitCountProbabilities recurrentVisitChain 2 InfiniteVisits)))
                        [2 / 3, 1 / 3, 1, 0]
                ]
            sequence_
                [ actual `shouldSatisfy` closeTo expected
                | (actual, expected) <-
                    zip
                        (entries (checked (visitCountProbabilities recurrentVisitChain 2 (FiniteVisits 0))))
                        [1 / 3, 2 / 3, 0, 1]
                ]

        it "counts the target at time zero" $ do
            checked (visitCountProbability transientVisitChain 0 (FiniteVisits 0) 0)
                `shouldBe` 0
            checked (visitCountProbability transientVisitChain 0 (FiniteVisits 1) 0)
                `shouldSatisfy` closeTo (3 / 4)

        prop "scalar queries look up the all-state result (random @3)" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left err -> counterexample (show err) False
                    Right matrix ->
                        conjoin
                            [ case visitCountProbabilities matrix 0 outcome of
                                Left err -> counterexample (show err) False
                                Right probabilities ->
                                    conjoin
                                        [ visitCountProbability matrix 0 outcome initial
                                            === Right probability
                                        | (initial, probability) <-
                                            zip (finites :: [Finite 3]) (entries probabilities)
                                        ]
                            | outcome <-
                                [ FiniteVisits 0
                                , FiniteVisits 1
                                , FiniteVisits 3
                                , InfiniteVisits
                                ]
                            ]

    describe "visitCountExpectations" $ do
        it "matches h / (1 - f) for a transient target" $ do
            let means = checked (visitCountExpectations transientVisitChain 0)
            sequence_
                [ mean `shouldSatisfy` meanCloseTo expected
                | (mean, expected) <- zip means [4 / 3, 2 / 3, 0]
                ]
            checked (visitCountExpectation transientVisitChain 0 1)
                `shouldSatisfy` meanCloseTo (2 / 3)

        it "is infinite exactly where a recurrent target is reachable" $ do
            checked (visitCountExpectations recurrentVisitChain 2)
                `shouldBe` [InfiniteMeanCount, InfiniteMeanCount, InfiniteMeanCount, FiniteMeanCount 0]
            checked (visitCountExpectation recurrentVisitChain 2 3)
                `shouldBe` FiniteMeanCount 0

        prop "agrees with hitting, return, recurrence, and reachability (random @3)" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left err -> counterexample (show err) False
                    Right matrix ->
                        case do
                            hits <- hittingProbabilities matrix [0]
                            returning <- returnProbability matrix 0
                            zeroVisits <- visitCountProbabilities matrix 0 (FiniteVisits 0)
                            oneVisit <- visitCountProbabilities matrix 0 (FiniteVisits 1)
                            twoVisits <- visitCountProbabilities matrix 0 (FiniteVisits 2)
                            infiniteVisits <- visitCountProbabilities matrix 0 InfiniteVisits
                            means <- visitCountExpectations matrix 0
                            pure
                                ( hits
                                , returning
                                , zeroVisits
                                , oneVisit
                                , twoVisits
                                , infiniteVisits
                                , means
                                ) of
                            Left err -> counterexample (show err) False
                            Right
                                ( hits
                                    , returning
                                    , zeroVisits
                                    , oneVisit
                                    , twoVisits
                                    , infiniteVisits
                                    , means
                                    ) ->
                                    let hitValues = entries hits
                                        zeroValues = entries zeroVisits
                                        oneValues = entries oneVisit
                                        twoValues = entries twoVisits
                                        infiniteValues = entries infiniteVisits
                                        states = finites :: [Finite 3]
                                        structuralMeans =
                                            [ if accessible matrix initial 0
                                                then InfiniteMeanCount
                                                else FiniteMeanCount 0
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
                                                        , means === structuralMeans
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
                                                            [ case mean of
                                                                FiniteMeanCount value ->
                                                                    property
                                                                        (closeTo (hit / (1 - returning)) value)
                                                                InfiniteMeanCount -> property False
                                                            | (hit, mean) <- zip hitValues means
                                                            ]
                                                        ]
                                            ]

    describe "visitCountDistributionBefore" $ do
        it "is a point mass at zero for bound zero" $
            distributionWeights
                ( visitCountDistributionBefore
                    0
                    (DistributionMap.pointMass (0 :: Finite 2))
                    twoCycle
                    (== 0)
                )
                `shouldBe` [(0, 1)]

        it "counts the initial state at a positive bound" $
            distributionWeights
                (visitCountDistributionBefore 1 mixedInitial twoCycle (== 0))
                `shouldBe` [(0, 0.75), (1, 0.25)]

        it "counts deterministic visits at times zero through bound minus one" $ do
            let initial = DistributionMap.pointMass (0 :: Finite 2)
            distributionWeights (visitCountDistributionBefore 1 initial twoCycle (== 0))
                `shouldBe` [(1, 1)]
            distributionWeights (visitCountDistributionBefore 2 initial twoCycle (== 0))
                `shouldBe` [(1, 1)]
            distributionWeights (visitCountDistributionBefore 3 initial twoCycle (== 0))
                `shouldBe` [(2, 1)]

        it "computes an exact law on an infinite random walk" $
            distributionWeights
                ( visitCountDistributionBefore
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
                                    visitCountDistributionBefore
                                        bound
                                        (DistributionMap.pointMass (0 :: Finite 3))
                                        matrix
                                        (== 0)
                                weights = distributionWeights law
                             in conjoin
                                    [ property (closeTo 1 (sum (map snd weights)))
                                    , property (all ((<= bound) . fst) weights)
                                    ]

    describe "visitCountProbabilityBefore" $ do
        it "looks up one coordinate of the count distribution" $ do
            let initial = DistributionMap.pointMass (0 :: Integer)
            visitCountProbabilityBefore 3 initial simpleRandomWalk (== 0) 1
                `shouldSatisfy` closeTo 0.5
            visitCountProbabilityBefore 3 initial simpleRandomWalk (== 0) 2
                `shouldSatisfy` closeTo 0.5
            visitCountProbabilityBefore 3 initial simpleRandomWalk (== 0) 3
                `shouldBe` 0

        it "agrees for a matrix and its equivalent kernel" $ do
            let initial = DistributionMap.pointMass (0 :: Finite 2)
                kernel = asKernel twoCycle
            sequence_
                [ visitCountProbabilityBefore bound initial twoCycle (== 0) count
                    `shouldBe` visitCountProbabilityBefore bound initial kernel (== 0) count
                | bound <- [0 .. 5]
                , count <- [0 .. bound]
                ]

    describe "visitCountExpectationBefore" $ do
        it "is the expectation of the random-walk count law" $
            visitCountExpectationBefore
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
                                    visitCountExpectationBefore bound initial matrix (== 0)
                                marginalSum =
                                    sum
                                        [ probabilityAtTime (fromIntegral time) initial matrix 0
                                        | time <- [0 .. rawBound - 1]
                                        ]
                             in property (closeTo marginalSum expectation)
