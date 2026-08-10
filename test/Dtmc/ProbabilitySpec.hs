{-# LANGUAGE TypeApplications #-}

module Dtmc.ProbabilitySpec (
    spec,
) where

import Data.Finite (
    Finite,
 )
import Data.List.NonEmpty (
    NonEmpty ((:|)),
 )
import Dtmc.Distribution (
    probabilityAt,
 )
import Dtmc.Distribution.Vector (
    DistributionVector,
    mkDistributionVector,
 )
import Dtmc.Probability (
    Observation (..),
    ProbabilityError (..),
    conditionalProbability,
    pathProbability,
    probability,
    probabilityAtTime,
    transitionProbability,
 )
import Dtmc.TestSupport (
    approxEq,
    assignment1Lambda,
    assignment1Matrix,
    genSimplexPoint,
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    mkTransitionMatrix,
 )
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
    conjoin,
    counterexample,
    forAll,
    property,
    (===),
 )

-- A three-state chain with several impossible one-step transitions.
chain :: TransitionMatrix (Finite 3)
chain =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 0.5
                , 0.5
                , 0.0
                , 0.0
                , 0.2
                , 0.8
                , 1.0
                , 0.0
                , 0.0
                ] ::
                S.Sq 3
            )

initial :: DistributionVector (Finite 3)
initial =
    either (error . show) id $
        mkDistributionVector (S.vector [0.6, 0.3, 0.1] :: S.R 3)

{- | Hold for a @Right@ whose 'Double' is within 'testTolerance' of the
expected value; fail for any @Left@ or out-of-tolerance value.
-}
rightCloseTo :: Double -> Either ProbabilityError Double -> Bool
rightCloseTo expected (Right actual) = approxEq testTolerance actual expected
rightCloseTo _ (Left _) = False

spec :: Spec
spec = do
    describe "Observation" $ do
        it "is polymorphic in the state type" $
            (At 2 "rain" :: Observation String) `shouldBe` At 2 "rain"

    describe "pathProbability" $ do
        it "returns the initial probability for a one-state path" $
            approxEq
                testTolerance
                (pathProbability initial chain (0 :| []))
                (probabilityAt initial 0)
                `shouldBe` True

        it "is lambda_i * P(i, j) for a two-state path" $
            approxEq
                testTolerance
                (pathProbability initial chain (0 :| [1]))
                (0.6 * 0.5)
                `shouldBe` True

        it "is the product of initial and transition probabilities" $
            approxEq
                testTolerance
                (pathProbability initial chain (0 :| [1, 2]))
                (0.6 * 0.5 * 0.8)
                `shouldBe` True

        it "is zero for a path with an impossible transition" $ do
            approxEq
                testTolerance
                (pathProbability initial chain (0 :| [2]))
                0
                `shouldBe` True
            approxEq
                testTolerance
                (pathProbability initial chain (0 :| [1, 0]))
                0
                `shouldBe` True

        prop "a one-state path equals the initial probability" $
            forAll ((,) <$> genSimplexPoint 3 <*> genTransitionMatrix @3) $
                \(entries, matrix) ->
                    case ( mkDistributionVector @(Finite 3) (S.vector entries :: S.R 3)
                         , mkTransitionMatrix matrix
                         ) of
                        (Right mu, Right p) ->
                            conjoin
                                [ pathProbability mu p (i :| [])
                                    === probabilityAt mu i
                                | i <- [0, 1, 2]
                                ]
                        result ->
                            counterexample
                                ("generated input was rejected: " <> show result)
                                False

        prop "a two-state path equals lambda_i * P(i, j)" $
            forAll ((,) <$> genSimplexPoint 3 <*> genTransitionMatrix @3) $
                \(entries, matrix) ->
                    case ( mkDistributionVector @(Finite 3) (S.vector entries :: S.R 3)
                         , mkTransitionMatrix matrix
                         ) of
                        (Right mu, Right p) ->
                            conjoin
                                [ property $
                                    approxEq
                                        testTolerance
                                        (pathProbability mu p (i :| [j]))
                                        ( probabilityAt mu i
                                            * transitionProbability p i j
                                        )
                                | i <- [0, 1, 2]
                                , j <- [0, 1, 2]
                                ]
                        result ->
                            counterexample
                                ("generated input was rejected: " <> show result)
                                False

    describe "probability" $ do
        it "returns exactly one for no observations" $
            probability initial chain [] `shouldBe` 1

        it "equals probabilityAtTime for a single observation" $
            approxEq
                testTolerance
                (probability initial chain [At 1 1])
                (probabilityAtTime 1 initial chain 1)
                `shouldBe` True

        it "is unchanged by observation order" $
            approxEq
                testTolerance
                (probability initial chain [At 0 0, At 1 1])
                (probability initial chain [At 1 1, At 0 0])
                `shouldBe` True

        it "is unchanged by duplicate observations" $
            approxEq
                testTolerance
                (probability initial chain [At 0 0, At 0 0, At 1 1])
                (probability initial chain [At 0 0, At 1 1])
                `shouldBe` True

        it "is exactly zero for conflicting states at one time" $
            probability initial chain [At 0 0, At 0 1] `shouldBe` 0

        it "agrees with pathProbability over times 0, 1, 2" $
            approxEq
                testTolerance
                (probability initial chain [At 0 0, At 1 1, At 2 2])
                (pathProbability initial chain (0 :| [1, 2]))
                `shouldBe` True

        it "is exactly zero through an impossible transition" $
            probability initial chain [At 0 0, At 1 2] `shouldBe` 0

        it "matches a hand-computed multi-gap example" $
            approxEq
                testTolerance
                ( probability
                    assignment1Lambda
                    assignment1Matrix
                    [At 2 2, At 3 4, At 6 3]
                )
                (5 / 96)
                `shouldBe` True

        prop "a single observation equals probabilityAtTime" $
            forAll ((,) <$> genSimplexPoint 3 <*> genTransitionMatrix @3) $
                \(entries, matrix) ->
                    case ( mkDistributionVector @(Finite 3) (S.vector entries :: S.R 3)
                         , mkTransitionMatrix matrix
                         ) of
                        (Right mu, Right p) ->
                            conjoin
                                [ property $
                                    approxEq
                                        testTolerance
                                        (probability mu p [At t i])
                                        (probabilityAtTime t mu p i)
                                | t <- [0, 1, 2]
                                , i <- [0, 1, 2]
                                ]
                        result ->
                            counterexample
                                ("generated input was rejected: " <> show result)
                                False

        prop "is invariant under observation order" $
            forAll ((,) <$> genSimplexPoint 3 <*> genTransitionMatrix @3) $
                \(entries, matrix) ->
                    case ( mkDistributionVector @(Finite 3) (S.vector entries :: S.R 3)
                         , mkTransitionMatrix matrix
                         ) of
                        (Right mu, Right p) ->
                            property $
                                approxEq
                                    testTolerance
                                    (probability mu p [At 1 1, At 3 2])
                                    (probability mu p [At 3 2, At 1 1])
                        result ->
                            counterexample
                                ("generated input was rejected: " <> show result)
                                False

    describe "conditionalProbability" $ do
        it "returns the event probability for an empty condition" $
            conditionalProbability initial chain [At 1 1] []
                `shouldSatisfy` rightCloseTo
                    (probability initial chain [At 1 1])

        it "returns one for an empty event and a positive condition" $
            conditionalProbability initial chain [] [At 0 0]
                `shouldSatisfy` rightCloseTo 1

        it "returns one when conditioning an observation on itself" $
            conditionalProbability initial chain [At 1 1] [At 1 1]
                `shouldSatisfy` rightCloseTo 1

        it "ignores observations shared by event and condition" $
            conditionalProbability initial chain [At 1 1] [At 1 1, At 2 2]
                `shouldSatisfy` rightCloseTo 1

        it "returns zero for a conflict against a possible condition" $
            conditionalProbability initial chain [At 1 0] [At 1 1]
                `shouldSatisfy` rightCloseTo 0

        it "reports a zero-probability condition" $
            conditionalProbability initial chain [At 0 0] [At 0 0, At 1 2]
                `shouldBe` Left ZeroProbabilityCondition

        it "reports a contradictory condition" $
            conditionalProbability initial chain [At 0 0] [At 1 1, At 1 2]
                `shouldBe` Left ZeroProbabilityCondition

        it "is unaffected by event and condition ordering" $ do
            conditionalProbability initial chain [At 2 2, At 1 1] [At 0 0]
                `shouldSatisfy` rightCloseTo 0.4
            conditionalProbability initial chain [At 1 1, At 2 2] [At 0 0]
                `shouldSatisfy` rightCloseTo 0.4

    describe "conditionalProbability assignment 1 regressions" $ do
        -- 1(b): P(X_10 = D, X_11 = D | X_3 = A, X_7 = E) = 1/4.
        it "1(b): P(X10=D, X11=D | X3=A, X7=E) = 1/4" $
            conditionalProbability
                assignment1Lambda
                assignment1Matrix
                [At 10 3, At 11 3]
                [At 3 0, At 7 4]
                `shouldSatisfy` rightCloseTo (1 / 4)

        -- 1(e): P(X_6 = D, X_2 = C | X_3 = E) = 15/92, event out of order.
        it "1(e): P(X6=D, X2=C | X3=E) = 15/92" $
            conditionalProbability
                assignment1Lambda
                assignment1Matrix
                [At 6 3, At 2 2]
                [At 3 4]
                `shouldSatisfy` rightCloseTo (15 / 92)

        it "1(e) supporting: P(X2=C) = 5/36" $
            approxEq
                testTolerance
                (probabilityAtTime 2 assignment1Lambda assignment1Matrix 2)
                (5 / 36)
                `shouldBe` True

        it "1(e) supporting: P(X3=E) = 23/72" $
            approxEq
                testTolerance
                (probabilityAtTime 3 assignment1Lambda assignment1Matrix 4)
                (23 / 72)
                `shouldBe` True

        it "1(e) supporting: P(X2=C, X3=E, X6=D) = 5/96" $
            approxEq
                testTolerance
                ( probability
                    assignment1Lambda
                    assignment1Matrix
                    [At 2 2, At 3 4, At 6 3]
                )
                (5 / 96)
                `shouldBe` True
