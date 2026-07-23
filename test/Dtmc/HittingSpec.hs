{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.HittingSpec (
    spec,
) where

import Data.Finite (
    Finite,
    finites,
 )
import Dtmc.Classification (
    accessible,
    recurrentState,
 )
import Dtmc.Hitting (
    MeanTime (..),
    expectedHittingTime,
    expectedHittingTimes,
    expectedReturnTime,
    expectedReturnTimes,
    hittingProbability,
    hittingProbabilities,
    returnProbability,
    returnProbabilities,
 )
import Dtmc.TestSupport (
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    mkTransitionMatrix,
    unTransitionMatrix,
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
    Property,
    conjoin,
    counterexample,
    forAll,
    property,
    (===),
 )

fromRows :: (Show e) => Either e (TransitionMatrix n) -> TransitionMatrix n
fromRows = either (error . show) id

-- Gambler's ruin on {0..4}: win 1 with probability p, lose 1 with
-- probability 1-p; 0 (ruin) and 4 (goal) are absorbing.
gambler :: Double -> TransitionMatrix 5
gambler p =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 1, 0, 0, 0, 0
                , 1 - p, 0, p, 0, 0
                , 0, 1 - p, 0, p, 0
                , 0, 0, 1 - p, 0, p
                , 0, 0, 0, 0, 1
                ]
            )

-- Oscillator: states 0 and 1 swap with probability 1/2 or exit to
-- their own absorbing state (0 -> 2, 1 -> 3).
oscillator :: TransitionMatrix 4
oscillator =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 0, 0.5, 0.5, 0
                , 0.5, 0, 0, 0.5
                , 0, 0, 1, 0
                , 0, 0, 0, 1
                ]
            )

twoCycle :: TransitionMatrix 2
twoCycle =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 0, 1
                , 1, 0
                ]
            )

nonUniformRecurrent :: TransitionMatrix 2
nonUniformRecurrent =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 0.9, 0.1
                , 0.4, 0.6
                ]
            )

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

meanCloseTo :: Double -> MeanTime -> Bool
meanCloseTo expected (FiniteMean v) = closeTo expected v
meanCloseTo _ InfiniteMean = False

checkedChain ::
    (KnownNat n) =>
    S.Sq n ->
    (TransitionMatrix n -> Property) ->
    Property
checkedChain matrix check =
    case mkTransitionMatrix matrix of
        Right p -> check p
        Left err ->
            counterexample ("generated matrix was rejected: " <> show err) False

spec :: Spec
spec = do
    describe "hittingProbabilities" $ do
        it "matches the gambler's ruin closed form (p = 0.4)" $ do
            let h = entries (hittingProbabilities (gambler 0.4) [0])
            length h `shouldBe` 5
            sequence_
                [ x `shouldSatisfy` closeTo (ruinProbability 0.4 i)
                | (i, x) <- zip [0 ..] h
                ]

        it "matches the symmetric closed form 1 - i/4 (p = 0.5)" $ do
            let h = entries (hittingProbabilities (gambler 0.5) [0])
            sequence_
                [ x `shouldSatisfy` closeTo (1 - fromIntegral i / 4)
                | (i, x) <- zip [0 :: Int ..] h
                ]

        it "solves the oscillator race to a single absorbing state" $ do
            let h = entries (hittingProbabilities oscillator [2])
            sequence_
                [ x `shouldSatisfy` closeTo v
                | (x, v) <- zip h [2 / 3, 1 / 3, 1, 0]
                ]

        it "is all zero for an empty target" $
            entries (hittingProbabilities oscillator [])
                `shouldBe` [0, 0, 0, 0]

        it "supports a single-state lookup without changing the result" $
            hittingProbability oscillator [2] 0
                `shouldSatisfy` closeTo (2 / 3)

        prop "is exactly one on the target and zero off its basin (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    let h = entries (hittingProbabilities p [0])
                     in conjoin
                            [ counterexample (show (i, x)) $
                                if
                                    | i == 0 -> x === 1
                                    | accessible p i 0 ->
                                        property
                                            (x >= -testTolerance && x <= 1 + testTolerance)
                                    | otherwise -> x === 0
                            | (i, x) <- zip (finites :: [Finite 4]) h
                            ]

        prop "satisfies the first-step equations off the target (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    let h = hittingProbabilities p [0]
                        pushed =
                            LA.toList
                                (S.extract (unTransitionMatrix p) LA.#> S.extract h)
                     in conjoin
                            [ property (closeTo hi pi_)
                            | (i, hi, pi_) <-
                                zip3 (finites :: [Finite 4]) (entries h) pushed
                            , i /= 0
                            ]

    describe "expectedHittingTimes" $ do
        it "returns one entry per state" $ do
            -- The transient entries come from the linear solve, so they are
            -- compared within tolerance; the target entries are assigned
            -- exactly and checked exactly.
            let eta = expectedHittingTimes oscillator [2, 3]
            sequence_
                [ e `shouldSatisfy` meanCloseTo 2
                | e <- take 2 eta
                ]
            drop 2 eta `shouldBe` [FiniteMean 0, FiniteMean 0]

        it "matches the gambler duration closed form (p = 0.4)" $ do
            let eta = expectedHittingTime (gambler 0.4) [0, 4]
            sequence_
                [ eta i `shouldSatisfy` meanCloseTo (ruinDuration 0.4 (fromIntegral i))
                | i <- finites :: [Finite 5]
                ]

        it "matches the symmetric duration i (4 - i) (p = 0.5)" $ do
            let eta = expectedHittingTime (gambler 0.5) [0, 4]
            sequence_
                [ eta i `shouldSatisfy` meanCloseTo (fromIntegral i * (4 - fromIntegral i))
                | i <- finites :: [Finite 5]
                ]

        it "expects two steps to absorption from either oscillator state" $ do
            let eta = expectedHittingTime oscillator [2, 3]
            eta 0 `shouldSatisfy` meanCloseTo 2
            eta 1 `shouldSatisfy` meanCloseTo 2
            eta 2 `shouldBe` FiniteMean 0
            eta 3 `shouldBe` FiniteMean 0

        it "is infinite when a competing absorbing state is reachable" $ do
            let eta = expectedHittingTime oscillator [2]
            eta 0 `shouldBe` InfiniteMean
            eta 1 `shouldBe` InfiniteMean
            eta 2 `shouldBe` FiniteMean 0
            eta 3 `shouldBe` InfiniteMean

        prop "finite entries satisfy the first-step equations (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    let eta = expectedHittingTime p [0]
                        rows = LA.toLists (S.extract (unTransitionMatrix p))
                        firstStep i row =
                            case eta i of
                                InfiniteMean -> property True
                                FiniteMean e ->
                                    case successorMeans row of
                                        Nothing ->
                                            counterexample
                                                "finite state with doomed successor"
                                                False
                                        Just total ->
                                            property (closeTo e (1 + total))
                        successorMeans row =
                            sum
                                <$> sequence
                                    [ case eta j of
                                        FiniteMean e -> Just (pij * e)
                                        InfiniteMean -> Nothing
                                    | (j, pij) <- zip (finites :: [Finite 4]) row
                                    , pij > 0
                                    , j /= 0
                                    ]
                     in conjoin
                            [ firstStep i row
                            | (i, row) <- zip (finites :: [Finite 4]) rows
                            , i /= 0
                            ]

    describe "returnProbabilities" $ do
        it "returns all state values in one solve" $ do
            -- The transient entries come from the fundamental-matrix solve,
            -- so they are compared within tolerance; the recurrent entries
            -- are assigned exactly one by the classification and checked
            -- exactly.
            let f = entries (returnProbabilities oscillator)
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
                    let rows = LA.toLists (S.extract (unTransitionMatrix p))
                     in conjoin
                            [ counterexample (show (i, f, firstStep)) $
                                property (closeTo firstStep f)
                            | (i, row, f) <-
                                zip3
                                    (finites :: [Finite 4])
                                    rows
                                    (entries (returnProbabilities p))
                            , let firstStep =
                                    sum
                                        ( zipWith
                                            (*)
                                            row
                                            (entries (hittingProbabilities p [i]))
                                        )
                            ]

        it "is one for an absorbing state" $
            returnProbability (gambler 0.5) 0 `shouldSatisfy` closeTo 1

        it "is one quarter for an oscillator state" $
            -- From 0: half the time exit to 2 (never return); otherwise reach
            -- 1, whence the return probability to 0 is 1/2. So f = 1/4.
            returnProbability oscillator 0 `shouldSatisfy` closeTo 0.25

        it "is one for both states of the two-cycle" $ do
            returnProbability twoCycle 0 `shouldSatisfy` closeTo 1
            returnProbability twoCycle 1 `shouldSatisfy` closeTo 1

        prop "is close to one on recurrent states and within [0, 1] (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    conjoin
                        [ counterexample (show (i, f)) $
                            property
                                ( f >= -testTolerance
                                    && f <= 1 + testTolerance
                                    && (not (recurrentState p i) || closeTo 1 f)
                                )
                        | i <- finites :: [Finite 4]
                        , let f = returnProbability p i
                        ]

    describe "expectedReturnTimes" $ do
        it "returns all state values in one table" $
            expectedReturnTimes oscillator
                `shouldBe` [InfiniteMean, InfiniteMean, FiniteMean 1, FiniteMean 1]

        it "is one for an absorbing state" $
            expectedReturnTime oscillator 2 `shouldBe` FiniteMean 1

        it "is two for either state of the two-cycle" $ do
            expectedReturnTime twoCycle 0 `shouldSatisfy` meanCloseTo 2
            expectedReturnTime twoCycle 1 `shouldSatisfy` meanCloseTo 2

        it "handles a non-uniform recurrent class" $ do
            expectedReturnTime nonUniformRecurrent 0 `shouldSatisfy` meanCloseTo 1.25
            expectedReturnTime nonUniformRecurrent 1 `shouldSatisfy` meanCloseTo 5

        it "is infinite for the oscillator's transient states" $ do
            expectedReturnTime oscillator 0 `shouldBe` InfiniteMean
            expectedReturnTime oscillator 1 `shouldBe` InfiniteMean

        prop "is finite exactly on recurrent states (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                checkedChain matrix $ \p ->
                    conjoin
                        [ counterexample (show i) $
                            isFinite (expectedReturnTime p i)
                                === recurrentState p i
                        | i <- finites :: [Finite 4]
                        ]
  where
    isFinite (FiniteMean _) = True
    isFinite InfiniteMean = False
