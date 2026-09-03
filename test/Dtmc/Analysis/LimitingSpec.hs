{-# LANGUAGE DataKinds #-}

module Dtmc.Analysis.LimitingSpec (
    spec,
) where

import Data.Finite (
    Finite,
 )
import Dtmc.Analysis.Limiting (
    converges,
    cyclicLimits,
    limitingMatrix,
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.TestSupport (
    approxEq,
    testTolerance,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    identityMatrix,
    matrixPower,
    mkTransitionMatrix,
    unTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )
import Test.Hspec (
    Spec,
    describe,
    expectationFailure,
    it,
    shouldBe,
    shouldSatisfy,
 )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

-- Section 4.2: closed classes {0} and {1,2}, both aperiodic.
twoClosedClasses :: TransitionMatrix (Finite 3)
twoClosedClasses =
    checked
        ( mkTransitionMatrix
            (S.matrix [1, 0, 0, 0, 0.4, 0.6, 0, 0.5, 0.5] :: S.Sq 3)
        )

-- State 0 is transient; {1,2} is the only recurrent class.
withTransient :: TransitionMatrix (Finite 3)
withTransient =
    checked
        ( mkTransitionMatrix
            (S.matrix [0, 0.5, 0.5, 0, 0.4, 0.6, 0, 0.5, 0.5] :: S.Sq 3)
        )

-- Irreducible, aperiodic, stationary distribution (0.8, 0.2).
twoState :: TransitionMatrix (Finite 2)
twoState =
    checked
        ( mkTransitionMatrix
            (S.matrix [0.9, 0.1, 0.4, 0.6] :: S.Sq 2)
        )

-- Irreducible with period 3, so P^n never settles.
threeCycle :: TransitionMatrix (Finite 3)
threeCycle =
    checked
        ( mkTransitionMatrix
            (S.matrix [0, 1, 0, 0, 0, 1, 1, 0, 0] :: S.Sq 3)
        )

-- Reducible with disjoint recurrent cycles of periods 2 and 3.
mixedPeriods :: TransitionMatrix (Finite 5)
mixedPeriods =
    checked
        ( mkTransitionMatrix
            ( S.matrix
                [ 0
                , 1
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 1
                , 0
                , 0
                ] ::
                S.Sq 5
            )
        )

-- State 0 is transient and enters the recurrent period-2 class {1,2}.
withTransientCycle :: TransitionMatrix (Finite 3)
withTransientCycle =
    checked
        ( mkTransitionMatrix
            (S.matrix [0, 1, 0, 0, 0, 1, 0, 1, 0] :: S.Sq 3)
        )

powerRows :: (FiniteState state) => Natural -> TransitionMatrix state -> [[Double]]
powerRows steps p =
    LA.toLists (S.extract (unTransitionMatrix (matrixPower steps p)))

matrixCloseTo :: [[Double]] -> [[Double]] -> Bool
matrixCloseTo expected actual =
    length expected == length actual
        && and (zipWith rowCloseTo expected actual)
  where
    rowCloseTo e a =
        length e == length a && and (zipWith (approxEq testTolerance) e a)

spec :: Spec
spec = do
    describe "converges" $ do
        it "accepts an aperiodic irreducible chain" $
            converges twoState `shouldBe` True

        it "accepts several aperiodic recurrent classes" $
            converges twoClosedClasses `shouldBe` True

        it "rejects a periodic class" $
            converges threeCycle `shouldBe` False

    describe "limitingMatrix" $ do
        it "matches the closed form of the notes" $
            case limitingMatrix twoClosedClasses of
                Right (Just rows) ->
                    rows
                        `shouldSatisfy` matrixCloseTo
                            [ [1, 0, 0]
                            , [0, 5 / 11, 6 / 11]
                            , [0, 5 / 11, 6 / 11]
                            ]
                other -> expectationFailure ("unexpected result: " ++ show other)

        it "repeats the stationary distribution in every row of an ergodic chain" $
            case limitingMatrix twoState of
                Right (Just rows) ->
                    rows `shouldSatisfy` matrixCloseTo [[0.8, 0.2], [0.8, 0.2]]
                other -> expectationFailure ("unexpected result: " ++ show other)

        it "is exactly zero on a transient column" $
            case limitingMatrix withTransient of
                Right (Just rows) -> do
                    map (take 1) rows `shouldBe` [[0], [0], [0]]
                    rows
                        `shouldSatisfy` matrixCloseTo
                            [ [0, 5 / 11, 6 / 11]
                            , [0, 5 / 11, 6 / 11]
                            , [0, 5 / 11, 6 / 11]
                            ]
                other -> expectationFailure ("unexpected result: " ++ show other)

        it "reports that a periodic chain has no limit" $
            limitingMatrix threeCycle `shouldBe` Right Nothing

        it "agrees with a high matrix power" $ do
            -- An independent route: repeated squaring rather than the
            -- hitting/stationary decomposition.
            case limitingMatrix twoState of
                Right (Just rows) ->
                    rows `shouldSatisfy` matrixCloseTo (powerRows 256 twoState)
                other -> expectationFailure ("unexpected result: " ++ show other)
            case limitingMatrix twoClosedClasses of
                Right (Just rows) ->
                    rows `shouldSatisfy` matrixCloseTo (powerRows 256 twoClosedClasses)
                other -> expectationFailure ("unexpected result: " ++ show other)

    describe "cyclicLimits" $ do
        it "returns one limit per period and reproduces the powers" $
            case cyclicLimits threeCycle of
                Right [atZero, atOne, atTwo] -> do
                    atZero `shouldSatisfy` matrixCloseTo (powerRows 3 threeCycle)
                    atOne `shouldSatisfy` matrixCloseTo (powerRows 4 threeCycle)
                    atTwo `shouldSatisfy` matrixCloseTo (powerRows 5 threeCycle)
                other -> expectationFailure ("expected three limits: " ++ show other)

        it "collapses to the ordinary limit when aperiodic" $
            case (cyclicLimits twoState, limitingMatrix twoState) of
                (Right [only], Right (Just rows)) ->
                    only `shouldSatisfy` matrixCloseTo rows
                other -> expectationFailure ("unexpected result: " ++ show other)

        it "collapses to the ordinary limit for a reducible aperiodic chain" $
            case (cyclicLimits twoClosedClasses, limitingMatrix twoClosedClasses) of
                (Right [only], Right (Just rows)) ->
                    only `shouldSatisfy` matrixCloseTo rows
                other -> expectationFailure ("unexpected result: " ++ show other)

        it "uses the least common multiple of recurrent class periods" $
            case cyclicLimits mixedPeriods of
                Right limits -> do
                    length limits `shouldBe` 6
                    and
                        ( zipWith
                            matrixCloseTo
                            [powerRows r mixedPeriods | r <- [0 .. 5]]
                            limits
                        )
                        `shouldBe` True
                other -> expectationFailure ("unexpected result: " ++ show other)

        it "accounts for the entry phase of transient states" $
            case cyclicLimits withTransientCycle of
                Right [atZero, atOne] -> do
                    atZero `shouldSatisfy` matrixCloseTo (powerRows 100 withTransientCycle)
                    atOne `shouldSatisfy` matrixCloseTo (powerRows 101 withTransientCycle)
                other -> expectationFailure ("expected two limits: " ++ show other)

        it "returns one empty limit for the empty chain" $
            cyclicLimits (identityMatrix @(Finite 0)) `shouldBe` Right [[]]
