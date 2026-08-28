{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.TimeSpecSupport (
    hittingTimeSpec,
    returnTimeSpec,
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
import Dtmc.Analysis.HittingTime (
    Expectation (..),
    LinearSystemError (..),
 )
import Dtmc.Analysis.HittingTime qualified as Hit
import Dtmc.Analysis.ReturnTime qualified as Return
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.State (
    FiniteState,
    finiteStates,
 )
import Dtmc.TestSupport (
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.Transition.Kernel qualified as Kernel
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    identityMatrix,
    mkTransitionMatrix,
    unTransitionMatrix,
 )
import GHC.Generics (
    Generic,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    expectationFailure,
    it,
    shouldBe,
    shouldSatisfy,
 )
import Test.Hspec.QuickCheck (
    prop,
 )
import Test.QuickCheck (
    Property,
    conjoin,
    counterexample,
    forAll,
    property,
    (===),
 )

data NamedRuinState = Ruined | One | Two | Three | Won
    deriving (Eq, Ord, Show, Generic)

instance FiniteState NamedRuinState

fromRows :: (Show e) => Either e (TransitionMatrix (Finite n)) -> TransitionMatrix (Finite n)
fromRows = either (error . show) id

asTransitionKernel ::
    (FiniteState state) =>
    TransitionMatrix state ->
    Kernel.TransitionKernel state
asTransitionKernel matrix =
    Kernel.transitionKernel $ \source ->
        either (error . show) id $
            DistributionMap.mkDistributionMap
                [ (destination, FT.stepProbability matrix source destination)
                | destination <- finiteStates
                ]

simpleRandomWalk :: Kernel.TransitionKernel Integer
simpleRandomWalk =
    Kernel.transitionKernel $ \state ->
        either (error . show) id $
            DistributionMap.mkDistributionMap [(state - 1, 0.5), (state + 1, 0.5)]

-- Gambler's ruin on {0..4}: win 1 with probability p, lose 1 with
-- probability 1-p; 0 (ruin) and 4 (goal) are absorbing.
gambler :: Double -> TransitionMatrix (Finite 5)
gambler p =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 1
                , 0
                , 0
                , 0
                , 0
                , 1 - p
                , 0
                , p
                , 0
                , 0
                , 0
                , 1 - p
                , 0
                , p
                , 0
                , 0
                , 0
                , 1 - p
                , 0
                , p
                , 0
                , 0
                , 0
                , 0
                , 1
                ]
            )

namedGambler :: TransitionMatrix NamedRuinState
namedGambler =
    either (error . show) id $
        mkTransitionMatrix @NamedRuinState
            ( S.matrix
                [ 1
                , 0
                , 0
                , 0
                , 0
                , 0.5
                , 0
                , 0.5
                , 0
                , 0
                , 0
                , 0.5
                , 0
                , 0.5
                , 0
                , 0
                , 0
                , 0.5
                , 0
                , 0.5
                , 0
                , 0
                , 0
                , 0
                , 1
                ] ::
                S.Sq 5
            )

-- Oscillator: states 0 and 1 swap with probability 1/2 or exit to
-- their own absorbing state (0 -> 2, 1 -> 3).
oscillator :: TransitionMatrix (Finite 4)
oscillator =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 0
                , 0.5
                , 0.5
                , 0
                , 0.5
                , 0
                , 0
                , 0.5
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 1
                ]
            )

twoCycle :: TransitionMatrix (Finite 2)
twoCycle =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 0
                , 1
                , 1
                , 0
                ]
            )

nonUniformRecurrent :: TransitionMatrix (Finite 2)
nonUniformRecurrent =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 0.9
                , 0.1
                , 0.4
                , 0.6
                ]
            )

-- 0 -> 1 -> 2 (absorbing): reaching 2 requires passing through 1 first.
pathChain :: TransitionMatrix (Finite 3)
pathChain =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 0
                , 1
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 1
                ]
            )

-- Two transient equations with very different scales. The system is
-- nonsingular in exact arithmetic but too ill-conditioned for the public
-- Double-precision numerical contract.
illConditionedChain :: TransitionMatrix (Finite 3)
illConditionedChain =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 1 - epsilon
                , 0
                , epsilon
                , 0
                , 0
                , 1
                , 0
                , 0
                , 1
                ]
            )
  where
    epsilon = 1e-14

-- Ruin probability from i with N = 4: (r^i - r^N) / (1 - r^N), r = (1-p)/p.
-- Only for p /= 1/2 (the symmetric case is 1 - i/N).
ruinProbability :: Double -> Int -> Double
ruinProbability p i = (r ^^ i - r ^^ n) / (1 - r ^^ n)
  where
    r = (1 - p) / p
    n = 4 :: Int

-- Expected duration until absorption at 0 or 4, for p /= 1/2:
-- i/(q-p) - (N/(q-p)) (1 - r^i) / (1 - r^N), q = 1-p, r = q/p.
ruinDuration :: Double -> Int -> Double
ruinDuration p i =
    fromIntegral i / (q - p)
        - (fromIntegral n / (q - p)) * (1 - r ^^ i) / (1 - r ^^ n)
  where
    q = 1 - p
    r = q / p
    n = 4 :: Int

entries :: (KnownNat n) => S.R n -> [Double]
entries = LA.toList . S.extract

closeTo :: Double -> Double -> Bool
closeTo expected x = abs (x - expected) <= testTolerance

expectationCloseTo :: Double -> Expectation -> Bool
expectationCloseTo expected (FiniteExpectation v) = closeTo expected v
expectationCloseTo _ InfiniteExpectation = False

checkedChain ::
    (KnownNat n) =>
    S.Sq n ->
    (TransitionMatrix (Finite n) -> Property) ->
    Property
checkedChain matrix check =
    case mkTransitionMatrix matrix of
        Right p -> check p
        Left err ->
            counterexample ("generated matrix was rejected: " <> show err) False

hittingTimeSpec :: Spec
hittingTimeSpec = do
    describe "numerical analysis errors" $
        it "rejects an ill-conditioned eventual-hitting system explicitly" $
            Hit.eventualProbabilityByState illConditionedChain [2]
                `shouldSatisfy` isIllConditioned

    describe "eventual hitting probability" $ do
        it "matches the gambler's ruin closed form (p = 0.4)" $ do
            case Hit.eventualProbabilityByState (gambler 0.4) [0] of
                Left err -> expectationFailure (show err)
                Right result -> do
                    let h = entries result
                    length h `shouldBe` 5
                    sequence_
                        [ x `shouldSatisfy` closeTo (ruinProbability 0.4 i)
                        | (i, x) <- zip [0 ..] h
                        ]

        it "matches the symmetric closed form 1 - i/4 (p = 0.5)" $ do
            case Hit.eventualProbabilityByState (gambler 0.5) [0] of
                Left err -> expectationFailure (show err)
                Right result ->
                    sequence_
                        [ x `shouldSatisfy` closeTo (1 - fromIntegral i / 4)
                        | (i, x) <- zip [0 :: Int ..] (entries result)
                        ]

        it "solves the oscillator race to a single absorbing state" $ do
            case Hit.eventualProbabilityByState oscillator [2] of
                Left err -> expectationFailure (show err)
                Right result ->
                    sequence_
                        [ x `shouldSatisfy` closeTo v
                        | (x, v) <- zip (entries result) [2 / 3, 1 / 3, 1, 0]
                        ]

        it "is all zero for an empty target" $
            (entries <$> Hit.eventualProbabilityByState oscillator [])
                `shouldBe` Right [0, 0, 0, 0]

        it "supports a single-state lookup without changing the result" $
            Hit.eventualProbability oscillator [2] 0
                `shouldSatisfy` either (const False) (closeTo (2 / 3))

        prop "is exactly one on the target and zero off its basin (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    case Hit.eventualProbabilityByState p [0] of
                        Left err -> counterexample (show err) False
                        Right result ->
                            conjoin
                                [ counterexample (show (i, x)) $
                                    if
                                        | i == 0 -> x === 1
                                        | accessible p i 0 ->
                                            property
                                                (x >= -testTolerance && x <= 1 + testTolerance)
                                        | otherwise -> x === 0
                                | (i, x) <-
                                    zip (finites :: [Finite 4]) (entries result)
                                ]

        prop "satisfies the first-step equations off the target (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    case Hit.eventualProbabilityByState p [0] of
                        Left err -> counterexample (show err) False
                        Right h ->
                            let pushed =
                                    LA.toList
                                        ( S.extract (unTransitionMatrix p)
                                            LA.#> S.extract h
                                        )
                             in conjoin
                                    [ property (closeTo hi pi_)
                                    | (i, hi, pi_) <-
                                        zip3 (finites :: [Finite 4]) (entries h) pushed
                                    , i /= 0
                                    ]

    describe "bounded hitting times" $ do
        it "returns an empty result for the empty chain" $
            entries ((Hit.probabilityByState . LessThan) 3 (identityMatrix @(Finite 0)) [])
                `shouldBe` []

        it "places all time-zero mass on the target" $
            entries ((Hit.probabilityByState . EqualTo) 0 oscillator [2])
                `shouldBe` [0, 0, 1, 0]

        it "gives zero exact-time mass for an empty target" $
            entries ((Hit.probabilityByState . EqualTo) 5 oscillator [])
                `shouldBe` [0, 0, 0, 0]

        it "matches a one-step gambler's-ruin hit" $
            entries ((Hit.probabilityByState . EqualTo) 1 (gambler 0.5) [0])
                `shouldBe` [0, 0.5, 0, 0, 0]

        it "uses a strict time bound" $ do
            entries ((Hit.probabilityByState . LessThan) 0 oscillator [2])
                `shouldBe` [0, 0, 0, 0]
            entries ((Hit.probabilityByState . LessThan) 1 oscillator [2])
                `shouldBe` [0, 0, 1, 0]
            (Hit.probability . LessThan) 2 (gambler 0.5) (== 0) 1
                `shouldSatisfy` closeTo 0.5

        it "ignores duplicate and reordered targets" $
            entries ((Hit.probabilityByState . LessThan) 4 oscillator [2, 3, 2])
                `shouldBe` entries ((Hit.probabilityByState . LessThan) 4 oscillator [3, 2])

        it "single-state queries look up the all-state results" $ do
            let exact = entries ((Hit.probabilityByState . EqualTo) 3 oscillator [2])
                bounded = entries ((Hit.probabilityByState . LessThan) 4 oscillator [2])
            sequence_
                [ (Hit.probability . EqualTo) 3 oscillator (== 2) i
                    `shouldSatisfy` closeTo exactAt
                | (i, exactAt) <- zip (finites :: [Finite 4]) exact
                ]
            sequence_
                [ (Hit.probability . LessThan) 4 oscillator (== 2) i
                    `shouldSatisfy` closeTo boundedAt
                | (i, boundedAt) <- zip (finites :: [Finite 4]) bounded
                ]

        prop "bounded increments equal exact-time mass (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    conjoin
                        [ counterexample (show (t, i, before, after, mass)) $
                            property (closeTo mass (after - before))
                        | t <- [0 .. 4]
                        , i <- finites :: [Finite 4]
                        , let before = (Hit.probability . LessThan) t p (== 0) i
                        , let after = (Hit.probability . LessThan) (t + 1) p (== 0) i
                        , let mass = (Hit.probability . EqualTo) t p (== 0) i
                        ]

        prop "bounded probabilities increase toward the eventual value (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    conjoin
                        [ counterexample (show (bound, i, current, next, eventual)) $
                            case eventual of
                                Left err -> counterexample (show err) False
                                Right value ->
                                    property
                                        ( current >= -testTolerance
                                            && current <= next + testTolerance
                                            && next <= value + testTolerance
                                        )
                        | bound <- [0 .. 4]
                        , i <- finites :: [Finite 4]
                        , let current = (Hit.probability . LessThan) bound p (== 0) i
                        , let next = (Hit.probability . LessThan) (bound + 1) p (== 0) i
                        , let eventual = Hit.eventualProbability p [0] i
                        ]

    describe "hitting race probability" $ do
        it "is exactly one on an effective successful state" $
            Hit.raceProbability (gambler 0.5) [4] [0] 4 `shouldBe` Right 1

        it "is exactly zero on a competing state" $
            Hit.raceProbability (gambler 0.5) [4] [0] 0 `shouldBe` Right 0

        it "is exactly zero on an overlapping (tied) state" $
            -- State 2 is in both boundaries, so the tie loses: value zero.
            Hit.raceProbability oscillator [2] [2, 3] 2
                `shouldBe` Right 0

        it "gives all zeros for identical successful and competing sets" $
            ( entries
                <$> Hit.raceProbabilityByState
                    oscillator
                    [2, 3]
                    [2, 3]
            )
                `shouldBe` Right (replicate 4 0)

        it "gives all zeros for an empty successful set" $
            (entries <$> Hit.raceProbabilityByState oscillator [] [2, 3])
                `shouldBe` Right (replicate 4 0)

        it "agrees with eventual hitting for an empty competing set" $ do
            case ( Hit.raceProbabilityByState
                    oscillator
                    [2, 3]
                    []
                 , Hit.eventualProbabilityByState oscillator [2, 3]
                 ) of
                (Left err, _) -> expectationFailure (show err)
                (_, Left err) -> expectationFailure (show err)
                (Right before, Right plain) ->
                    sequence_
                        [ x `shouldSatisfy` closeTo y
                        | (x, y) <- zip (entries before) (entries plain)
                        ]

        it "is exactly zero when the successful set is unreachable" $
            -- Absorbing state 3 cannot reach absorbing state 2.
            Hit.raceProbability oscillator [2] [] 3
                `shouldBe` Right 0

        it "is exactly zero when success needs a competitor first" $
            -- 0 -> 1 -> 2 with 1 competing: 2 is reachable only through 1.
            Hit.raceProbability pathChain [2] [1] 0 `shouldBe` Right 0

        it "ignores duplicate targets" $ do
            case ( Hit.raceProbabilityByState
                    oscillator
                    [2, 2]
                    [3, 3]
                 , Hit.raceProbabilityByState
                    oscillator
                    [2]
                    [3]
                 ) of
                (Left err, _) -> expectationFailure (show err)
                (_, Left err) -> expectationFailure (show err)
                (Right withDuplicates, Right once) ->
                    sequence_
                        [ x `shouldSatisfy` closeTo y
                        | (x, y) <- zip (entries withDuplicates) (entries once)
                        ]

        it "ignores target order" $ do
            case ( Hit.raceProbabilityByState
                    oscillator
                    [2, 0]
                    [3, 1]
                 , Hit.raceProbabilityByState
                    oscillator
                    [0, 2]
                    [1, 3]
                 ) of
                (Left err, _) -> expectationFailure (show err)
                (_, Left err) -> expectationFailure (show err)
                (Right reordered, Right ordered) ->
                    sequence_
                        [ x `shouldSatisfy` closeTo y
                        | (x, y) <- zip (entries reordered) (entries ordered)
                        ]

        it "single-state lookups match the all-state vector" $
            case Hit.raceProbabilityByState
                oscillator
                [2]
                [3] of
                Left err -> expectationFailure (show err)
                Right result ->
                    sequence_
                        [ Hit.raceProbability
                            oscillator
                            [2]
                            [3]
                            i
                            `shouldSatisfy` either (const False) (closeTo x)
                        | (i, x) <-
                            zip (finites :: [Finite 4]) (entries result)
                        ]

        it "solves the oscillator race against a competing absorber" $ do
            case Hit.raceProbabilityByState oscillator [2] [3] of
                Left err -> expectationFailure (show err)
                Right result ->
                    sequence_
                        [ x `shouldSatisfy` closeTo v
                        | (x, v) <- zip (entries result) [2 / 3, 1 / 3, 1, 0]
                        ]

        it "matches a hand-computed symmetric race (gambler p = 0.5)" $ do
            case Hit.raceProbabilityByState (gambler 0.5) [4] [0] of
                Left err -> expectationFailure (show err)
                Right result ->
                    sequence_
                        [ x `shouldSatisfy` closeTo (fromIntegral i / 4)
                        | (i, x) <- zip [0 :: Int ..] (entries result)
                        ]

        it "disjoint races sum to one when the union is hit almost surely" $
            sequence_
                [ case ( Hit.raceProbabilityByState g [4] [0]
                       , Hit.raceProbabilityByState g [0] [4]
                       ) of
                    (Left err, _) -> expectationFailure (show err)
                    (_, Left err) -> expectationFailure (show err)
                    (Right wins, Right losses) ->
                        sequence_
                            [ (x + y) `shouldSatisfy` closeTo 1
                            | (x, y) <- zip (entries wins) (entries losses)
                            ]
                | pp <- [0.3, 0.5, 0.7]
                , let g = gambler pp
                ]

    describe "expected hitting time" $ do
        it "returns one entry per state" $ do
            -- The transient entries come from the linear solve, so they are
            -- compared within tolerance; the target entries are assigned
            -- exactly and checked exactly.
            case Hit.expectationByState oscillator [2, 3] of
                Left err -> expectationFailure (show err)
                Right eta -> do
                    sequence_
                        [ e `shouldSatisfy` expectationCloseTo 2
                        | e <- take 2 eta
                        ]
                    drop 2 eta `shouldBe` [FiniteExpectation 0, FiniteExpectation 0]

        it "matches the gambler duration closed form (p = 0.4)" $ do
            let eta = Hit.expectation (gambler 0.4) [0, 4]
            sequence_
                [ eta i
                    `shouldSatisfy` either
                        (const False)
                        (expectationCloseTo (ruinDuration 0.4 (fromIntegral i)))
                | i <- finites :: [Finite 5]
                ]

        it "matches the symmetric duration i (4 - i) (p = 0.5)" $ do
            let eta = Hit.expectation (gambler 0.5) [0, 4]
            sequence_
                [ eta i
                    `shouldSatisfy` either
                        (const False)
                        (expectationCloseTo (fromIntegral i * (4 - fromIntegral i)))
                | i <- finites :: [Finite 5]
                ]

        it "expects two steps to absorption from either oscillator state" $ do
            let eta = Hit.expectation oscillator [2, 3]
            eta 0 `shouldSatisfy` either (const False) (expectationCloseTo 2)
            eta 1 `shouldSatisfy` either (const False) (expectationCloseTo 2)
            eta 2 `shouldBe` Right (FiniteExpectation 0)
            eta 3 `shouldBe` Right (FiniteExpectation 0)

        it "is infinite when a competing absorbing state is reachable" $ do
            let eta = Hit.expectation oscillator [2]
            eta 0 `shouldBe` Right InfiniteExpectation
            eta 1 `shouldBe` Right InfiniteExpectation
            eta 2 `shouldBe` Right (FiniteExpectation 0)
            eta 3 `shouldBe` Right InfiniteExpectation

        prop "finite entries satisfy the first-step equations (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    case Hit.expectationByState p [0] of
                        Left err -> counterexample (show err) False
                        Right times ->
                            let eta i = times !! fromIntegral i
                                rows = LA.toLists (S.extract (unTransitionMatrix p))
                                firstStep i row =
                                    case eta i of
                                        InfiniteExpectation -> property True
                                        FiniteExpectation e ->
                                            case successorExpectations row of
                                                Nothing ->
                                                    counterexample
                                                        "finite state with doomed successor"
                                                        False
                                                Just total ->
                                                    property (closeTo e (1 + total))
                                successorExpectations row =
                                    sum
                                        <$> sequence
                                            [ case eta j of
                                                FiniteExpectation e -> Just (pij * e)
                                                InfiniteExpectation -> Nothing
                                            | (j, pij) <-
                                                zip (finites :: [Finite 4]) row
                                            , pij > 0
                                            , j /= 0
                                            ]
                             in conjoin
                                    [ firstStep i row
                                    | (i, row) <-
                                        zip (finites :: [Finite 4]) rows
                                    , i /= 0
                                    ]

returnTimeSpec :: Spec
returnTimeSpec = do
    describe "bounded first-return times" $ do
        it "returns an empty result for the empty chain" $
            entries ((Return.probabilityByState . LessThan) 3 (identityMatrix @(Finite 0)))
                `shouldBe` []

        it "has no return mass at time zero" $
            entries ((Return.probabilityByState . EqualTo) 0 oscillator)
                `shouldBe` [0, 0, 0, 0]

        it "uses the transition diagonal at time one" $
            entries ((Return.probabilityByState . EqualTo) 1 nonUniformRecurrent)
                `shouldBe` [0.9, 0.6]

        it "counts only the first return" $ do
            entries ((Return.probabilityByState . EqualTo) 1 oscillator)
                `shouldBe` [0, 0, 1, 1]
            entries ((Return.probabilityByState . EqualTo) 2 oscillator)
                `shouldBe` [0.25, 0.25, 0, 0]
            entries ((Return.probabilityByState . EqualTo) 2 twoCycle)
                `shouldBe` [1, 1]

        it "uses a strict time bound" $ do
            entries ((Return.probabilityByState . LessThan) 0 oscillator)
                `shouldBe` [0, 0, 0, 0]
            entries ((Return.probabilityByState . LessThan) 1 oscillator)
                `shouldBe` [0, 0, 0, 0]
            entries ((Return.probabilityByState . LessThan) 2 oscillator)
                `shouldBe` [0, 0, 1, 1]
            entries ((Return.probabilityByState . LessThan) 3 twoCycle)
                `shouldBe` [1, 1]

        it "single-state queries look up the all-state results" $ do
            let exact = entries ((Return.probabilityByState . EqualTo) 3 oscillator)
                bounded = entries ((Return.probabilityByState . LessThan) 4 oscillator)
            sequence_
                [ (Return.probability . EqualTo) 3 oscillator i
                    `shouldSatisfy` closeTo exactAt
                | (i, exactAt) <- zip (finites :: [Finite 4]) exact
                ]
            sequence_
                [ (Return.probability . LessThan) 4 oscillator i
                    `shouldSatisfy` closeTo boundedAt
                | (i, boundedAt) <- zip (finites :: [Finite 4]) bounded
                ]

        prop "bounded increments equal exact-time mass (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    conjoin
                        [ counterexample (show (t, i, before, after, mass)) $
                            property (closeTo mass (after - before))
                        | t <- [0 .. 4]
                        , i <- finites :: [Finite 4]
                        , let before = (Return.probability . LessThan) t p i
                        , let after = (Return.probability . LessThan) (t + 1) p i
                        , let mass = (Return.probability . EqualTo) t p i
                        ]

        prop "bounded probabilities increase toward the eventual value (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    conjoin
                        [ counterexample (show (bound, i, current, next, eventual)) $
                            case eventual of
                                Left err -> counterexample (show err) False
                                Right value ->
                                    property
                                        ( current >= -testTolerance
                                            && current <= next + testTolerance
                                            && next <= value + testTolerance
                                        )
                        | bound <- [0 .. 4]
                        , i <- finites :: [Finite 4]
                        , let current = (Return.probability . LessThan) bound p i
                        , let next = (Return.probability . LessThan) (bound + 1) p i
                        , let eventual = Return.eventualProbability p i
                        ]

    describe "eventual return probability" $ do
        it "returns all state values in one solve" $ do
            -- The transient entries come from the fundamental-matrix solve,
            -- so they are compared within tolerance; the recurrent entries
            -- are assigned exactly one by the classification and checked
            -- exactly.
            case Return.eventualProbabilityByState oscillator of
                Left err -> expectationFailure (show err)
                Right result -> do
                    let f = entries result
                    sequence_
                        [ x `shouldSatisfy` closeTo 0.25
                        | x <- take 2 f
                        ]
                    drop 2 f `shouldBe` [1, 1]

        prop "agrees with the first-step decomposition (random @4)" $
            -- Two independent theorems for the same quantity: the
            -- implementation computes f_i = 1 - 1/N_ii from the renewal
            -- identity, while conditioning on the first step gives
            -- f_i = sum_j P_ij h_j{i}.
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    case Return.eventualProbabilityByState p of
                        Left err -> counterexample (show err) False
                        Right returns ->
                            let rows = LA.toLists (S.extract (unTransitionMatrix p))
                             in conjoin
                                    [ case Hit.eventualProbabilityByState p [i] of
                                        Left err -> counterexample (show err) False
                                        Right hits ->
                                            let firstStep =
                                                    sum
                                                        ( zipWith
                                                            (*)
                                                            row
                                                            (entries hits)
                                                        )
                                             in counterexample
                                                    (show (i, f, firstStep))
                                                    (property (closeTo firstStep f))
                                    | (i, row, f) <-
                                        zip3
                                            (finites :: [Finite 4])
                                            rows
                                            (entries returns)
                                    ]

        it "is one for an absorbing state" $
            Return.eventualProbability (gambler 0.5) 0
                `shouldSatisfy` either (const False) (closeTo 1)

        it "is one quarter for an oscillator state" $
            -- From 0: half the time exit to 2 (never return); otherwise reach
            -- 1, whence the return probability to 0 is 1/2. So f = 1/4.
            Return.eventualProbability oscillator 0
                `shouldSatisfy` either (const False) (closeTo 0.25)

        it "is one for both states of the two-cycle" $ do
            Return.eventualProbability twoCycle 0
                `shouldSatisfy` either (const False) (closeTo 1)
            Return.eventualProbability twoCycle 1
                `shouldSatisfy` either (const False) (closeTo 1)

        prop "is close to one on recurrent states and within [0, 1] (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    conjoin
                        [ counterexample (show (i, f)) $
                            case f of
                                Left err -> counterexample (show err) False
                                Right value ->
                                    property
                                        ( value >= -testTolerance
                                            && value <= 1 + testTolerance
                                            && ( not (recurrentState p i)
                                                    || closeTo 1 value
                                               )
                                        )
                        | i <- finites :: [Finite 4]
                        , let f = Return.eventualProbability p i
                        ]

    describe "expected return time" $ do
        it "returns all state values in one table" $
            Return.expectationByState oscillator
                `shouldBe` Right [InfiniteExpectation, InfiniteExpectation, FiniteExpectation 1, FiniteExpectation 1]

        it "is one for an absorbing state" $
            Return.expectation oscillator 2 `shouldBe` Right (FiniteExpectation 1)

        it "is two for either state of the two-cycle" $ do
            Return.expectation twoCycle 0
                `shouldSatisfy` either (const False) (expectationCloseTo 2)
            Return.expectation twoCycle 1
                `shouldSatisfy` either (const False) (expectationCloseTo 2)

        it "handles a non-uniform recurrent class" $ do
            Return.expectation nonUniformRecurrent 0
                `shouldSatisfy` either (const False) (expectationCloseTo 1.25)
            Return.expectation nonUniformRecurrent 1
                `shouldSatisfy` either (const False) (expectationCloseTo 5)

        it "is infinite for the oscillator's transient states" $ do
            Return.expectation oscillator 0 `shouldBe` Right InfiniteExpectation
            Return.expectation oscillator 1 `shouldBe` Right InfiniteExpectation

        prop "is finite exactly on recurrent states (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    conjoin
                        [ counterexample (show i) $
                            case Return.expectation p i of
                                Left err -> counterexample (show err) False
                                Right result ->
                                    isFinite result === recurrentState p i
                        | i <- finites :: [Finite 4]
                        ]

    describe "Transition realization independence" $ do
        it "uses strict hitting bounds on an infinite random walk" $ do
            (Hit.probability . EqualTo) 2 simpleRandomWalk (== 2) 0
                `shouldSatisfy` closeTo 0.25
            (Hit.probability . LessThan) 2 simpleRandomWalk (== 2) 0
                `shouldBe` 0
            (Hit.probability . LessThan) 3 simpleRandomWalk (== 2) 0
                `shouldSatisfy` closeTo 0.25

        it "distinguishes return time from time-zero hitting" $ do
            (Hit.probability . EqualTo) 0 simpleRandomWalk (== 0) 0
                `shouldBe` 1
            (Return.probability . EqualTo) 0 simpleRandomWalk 0
                `shouldBe` 0
            (Return.probability . EqualTo) 2 simpleRandomWalk 0
                `shouldSatisfy` closeTo 0.5
            (Return.probability . LessThan) 2 simpleRandomWalk 0
                `shouldBe` 0
            (Return.probability . LessThan) 3 simpleRandomWalk 0
                `shouldSatisfy` closeTo 0.5

        prop "matches matrix and equivalent-kernel bounded queries" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left problem -> counterexample (show problem) False
                    Right matrix ->
                        let kernel = asTransitionKernel matrix
                            target state = state == (2 :: Finite 3)
                         in conjoin
                                [ counterexample (show (state, time)) $
                                    property $
                                        and
                                            [ closeTo
                                                ((Hit.probability . EqualTo) time matrix target state)
                                                ((Hit.probability . EqualTo) time kernel target state)
                                            , closeTo
                                                ((Hit.probability . LessThan) time matrix target state)
                                                ((Hit.probability . LessThan) time kernel target state)
                                            , closeTo
                                                ((Return.probability . EqualTo) time matrix state)
                                                ((Return.probability . EqualTo) time kernel state)
                                            , closeTo
                                                ((Return.probability . LessThan) time matrix state)
                                                ((Return.probability . LessThan) time kernel state)
                                            ]
                                | state <- finites :: [Finite 3]
                                , time <- [0 .. 4]
                                ]

    describe "named finite states" $ do
        it "solves eventual and competing hitting queries by constructor" $ do
            case Hit.eventualProbabilityByState namedGambler [Ruined, Won] of
                Left err -> expectationFailure (show err)
                Right result ->
                    sequence_
                        [ probability `shouldSatisfy` closeTo 1
                        | probability <- entries result
                        ]
            Hit.raceProbability namedGambler [Won] [Ruined] Two
                `shouldSatisfy` either (const False) (closeTo 0.5)

        it "solves bounded hitting queries in named state order" $
            entries ((Hit.probabilityByState . LessThan) 3 namedGambler [Won])
                `shouldBe` [0, 0, 0.25, 0.5, 1]

        it "solves named expected hitting and return times" $ do
            Hit.expectation namedGambler [Ruined, Won] Two
                `shouldSatisfy` either (const False) (expectationCloseTo 4)
            Return.expectation namedGambler Ruined
                `shouldBe` Right (FiniteExpectation 1)

isFinite :: Expectation -> Bool
isFinite (FiniteExpectation _) = True
isFinite InfiniteExpectation = False

isIllConditioned :: Either LinearSystemError value -> Bool
isIllConditioned (Left (IllConditionedSystem estimate)) =
    estimate < 1e-12
isIllConditioned _ = False
