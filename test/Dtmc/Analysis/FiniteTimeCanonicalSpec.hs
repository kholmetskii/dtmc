{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.FiniteTimeCanonicalSpec (
    spec,
) where

import Data.Finite (
    Finite,
    finites,
 )
import Dtmc.Analysis.FiniteTime qualified as FT
import Dtmc.Analysis.ProbabilityOracle qualified as Oracle
import Dtmc.Distribution.Map (
    fromList,
    pointMass,
 )
import Dtmc.Distribution.Vector (
    DistributionVector,
 )
import Dtmc.Distribution.Vector.HMatrix (
    mkDistributionVector,
 )
import Dtmc.TestSupport (
    approxEq,
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.Transition.Kernel (
    TransitionKernel,
    transitionKernel,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
 )
import Dtmc.Transition.Matrix.HMatrix (
    mkTransitionMatrix,
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
    counterexample,
    forAll,
    property,
 )

initialWeights :: [(Finite 3, Double)]
initialWeights = zip finites [0.2, 0.3, 0.5]

initialDistribution :: DistributionVector (Finite 3)
initialDistribution =
    checked (mkDistributionVector (S.vector [0.2, 0.3, 0.5] :: S.R 3))

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

simpleRandomWalk :: TransitionKernel Integer
simpleRandomWalk =
    transitionKernel $ \state ->
        checked
            ( fromList
                [(state - 1, 0.5), (state + 1, 0.5)]
            )

close :: Double -> Double -> Bool
close = approxEq testTolerance

canonicalMatchesOracle :: TransitionMatrix (Finite 3) -> Bool
canonicalMatchesOracle matrix =
    and
        [ and
            [ close
                (FT.stepProbability matrix source destination)
                (Oracle.transitionWeight matrix source destination)
            | source <- finites
            , destination <- finites
            ]
        , and
            [ close
                (FT.nStepProbability time matrix source destination)
                (Oracle.stateProbability time [(source, 1)] matrix destination)
            | time <- [0 .. 4]
            , source <- finites
            , destination <- finites
            ]
        , and
            [ close
                (FT.probability initialDistribution matrix [FT.At time destination])
                (Oracle.stateProbability time initialWeights matrix destination)
            | time <- [0 .. 4]
            , destination <- finites
            ]
        , close
            ( FT.probability
                initialDistribution
                matrix
                [FT.At 0 0, FT.At 1 1, FT.At 2 2]
            )
            (Oracle.trajectoryProbability initialWeights matrix [0, 1, 2])
        , close
            ( FT.probability
                initialDistribution
                matrix
                [FT.At 1 1, FT.At 3 2]
            )
            ( Oracle.observationProbability
                3
                initialWeights
                matrix
                [(1, 1), (3, 2)]
            )
        , conditionalMatchesOracle
        ]
  where
    denominator =
        Oracle.observationProbability 1 initialWeights matrix [(1, 1)]
    numerator =
        Oracle.observationProbability
            3
            initialWeights
            matrix
            [(1, 1), (3, 2)]
    conditionalMatchesOracle =
        case FT.probabilityGiven
            initialDistribution
            matrix
            [FT.At 3 2]
            [FT.At 1 1] of
            Left FT.ZeroProbabilityCondition -> denominator == 0
            Right actual -> denominator /= 0 && close actual (numerator / denominator)

spec :: Spec
spec = do
    describe "canonical finite-time namespace" $ do
        it "uses the four grammar-compliant names together" $ do
            let matrix :: TransitionMatrix (Finite 2)
                matrix =
                    checked
                        ( mkTransitionMatrix
                            (S.matrix [0.5, 0.5, 0, 1] :: S.Sq 2)
                        )
                initial :: DistributionVector (Finite 2)
                initial =
                    checked
                        (mkDistributionVector (S.vector [1, 0] :: S.R 2))
            FT.stepProbability matrix 0 1 `shouldBe` 0.5
            FT.nStepProbability 2 matrix 0 1 `shouldBe` 0.75
            FT.probability initial matrix [FT.At 1 1] `shouldBe` 0.5
            FT.probability initial matrix [FT.At 0 0, FT.At 1 1]
                `shouldBe` 0.5
            FT.probabilityGiven
                initial
                matrix
                [FT.At 1 1]
                [FT.At 0 0]
                `shouldBe` Right 0.5

        it "preserves locally finite countable-state support" $ do
            FT.stepProbability simpleRandomWalk 0 1 `shouldBe` 0.5
            FT.nStepProbability 2 simpleRandomWalk 0 0 `shouldBe` 0.5
            FT.probability (pointMass 0) simpleRandomWalk [FT.At 2 0]
                `shouldBe` 0.5
            FT.probability
                (pointMass 0)
                simpleRandomWalk
                [FT.At 0 0, FT.At 1 1, FT.At 2 0]
                `shouldBe` 0.25
            FT.probability
                (pointMass 0)
                simpleRandomWalk
                [FT.At 0 0, FT.At 2 0]
                `shouldBe` 0.5
            FT.probabilityGiven
                (pointMass 0)
                simpleRandomWalk
                [FT.At 2 0]
                [FT.At 1 1]
                `shouldBe` Right 0.5

        prop "matches independent path enumeration (random @3)" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left problem -> counterexample (show problem) False
                    Right matrix -> property (canonicalMatchesOracle matrix)
