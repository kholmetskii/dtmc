{-# LANGUAGE TypeApplications #-}

module Dtmc.ProbabilitySpec (
    spec,
) where

import Data.List.NonEmpty (
    NonEmpty ((:|)),
 )
import Dtmc.Distribution (
    Distribution,
    mkDistribution,
    probabilityAt,
 )
import Dtmc.Probability (
    pathProbability,
 )
import Dtmc.TestSupport (
    approxEq,
    genSimplexPoint,
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    mkTransitionMatrix,
    transitionProbability,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
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
chain :: TransitionMatrix 3
chain =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 0.5, 0.5, 0.0
                , 0.0, 0.2, 0.8
                , 1.0, 0.0, 0.0
                ] ::
                S.Sq 3
            )

initial :: Distribution 3
initial =
    either (error . show) id $
        mkDistribution (S.vector [0.6, 0.3, 0.1] :: S.R 3)

spec :: Spec
spec = do
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

        prop "a one-state path equals the initial probability"
            $ forAll ((,) <$> genSimplexPoint 3 <*> genTransitionMatrix @3)
            $ \(entries, matrix) ->
                case
                    ( mkDistribution (S.vector entries :: S.R 3)
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

        prop "a two-state path equals lambda_i * P(i, j)"
            $ forAll ((,) <$> genSimplexPoint 3 <*> genTransitionMatrix @3)
            $ \(entries, matrix) ->
                case
                    ( mkDistribution (S.vector entries :: S.R 3)
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
