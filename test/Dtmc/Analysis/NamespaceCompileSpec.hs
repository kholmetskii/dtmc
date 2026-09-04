{-# LANGUAGE DataKinds #-}

module Dtmc.Analysis.NamespaceCompileSpec (
    spec,
) where

import Data.Finite (
    Finite,
 )
import Dtmc.Analysis.Absorption qualified as Absorption
import Dtmc.Analysis.Event (
    DiscreteEvent (..),
    matchesDiscreteEvent,
 )
import Dtmc.Analysis.Expectation (
    Expectation (..),
 )
import Dtmc.Analysis.FiniteTime qualified as FT
import Dtmc.Analysis.HittingTime qualified as Hit
import Dtmc.Analysis.ReturnTime qualified as Return
import Dtmc.Analysis.VisitCount qualified as Visit
import Dtmc.Distribution.Vector (
    DistributionVector,
 )
import Dtmc.Distribution.Vector.HMatrix (
    mkDistributionVector,
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

matrix :: TransitionMatrix (Finite 2)
matrix =
    checked
        ( mkTransitionMatrix
            (S.matrix [0.5, 0.5, 0, 1] :: S.Sq 2)
        )

initial :: DistributionVector (Finite 2)
initial =
    checked
        (mkDistributionVector (S.vector [1, 0] :: S.R 2))

mixedInitial :: DistributionVector (Finite 2)
mixedInitial =
    checked
        (mkDistributionVector (S.vector [0.25, 0.75] :: S.R 2))

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

spec :: Spec
spec = do
    describe "qualified analysis namespaces" $ do
        it "coexist under the documented aliases" $ do
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
            Hit.probability (EqualTo 1) matrix (== 1) initial `shouldBe` 0.5
            length
                [ Hit.probabilityGivenInitialState (AtMost 1) matrix (== 1) 0 `seq` ()
                , Hit.probability (AtMost 1) matrix (== 1) initial `seq` ()
                , Hit.eventualProbabilityGivenInitialState matrix [1] 0 `seq` ()
                , Hit.eventualProbability matrix [1] initial `seq` ()
                , Hit.raceProbabilityGivenInitialState matrix [1] [] 0 `seq` ()
                , Hit.raceProbability matrix [1] [] initial `seq` ()
                , Hit.expectationGivenInitialState matrix [1] 0 `seq` ()
                , Hit.expectation matrix [1] initial `seq` ()
                ]
                `shouldBe` 8
            Return.probability (EqualTo 1) matrix initial `shouldBe` 0.5
            length
                [ Return.probabilityGivenInitialState (AtMost 1) matrix 0 `seq` ()
                , Return.probability (AtMost 1) matrix initial `seq` ()
                , Return.eventualProbabilityGivenInitialState matrix 0 `seq` ()
                , Return.eventualProbability matrix initial `seq` ()
                , Return.expectationGivenInitialState matrix 0 `seq` ()
                , Return.expectation matrix initial `seq` ()
                ]
                `shouldBe` 6
            Visit.boundedProbability 2 (EqualTo 1) initial matrix (== 1)
                `shouldBe` 0.5
            length
                [ Visit.totalProbabilityGivenInitialState (EqualTo 1) matrix 1 0 `seq` ()
                , Visit.totalProbability (AtMost 1) matrix 1 initial `seq` ()
                , Visit.infiniteProbabilityGivenInitialState matrix 1 0 `seq` ()
                , Visit.infiniteProbability matrix 1 initial `seq` ()
                , Visit.totalExpectationGivenInitialState matrix 1 0 `seq` ()
                , Visit.totalExpectation matrix 1 initial `seq` ()
                , Visit.boundedLaw 2 initial matrix (== 1) `seq` ()
                , Visit.boundedProbability 2 (AtMost 1) initial matrix (== 1) `seq` ()
                , Visit.boundedProbabilityGivenInitialState 2 (AtMost 1) 0 matrix (== 1) `seq` ()
                , Visit.boundedExpectation 2 initial matrix (== 1) `seq` ()
                , Visit.boundedExpectationGivenInitialState 2 0 matrix (== 1) `seq` ()
                ]
                `shouldBe` 11
            matchesDiscreteEvent (AtMost 1) 1 `shouldBe` True

        it "distinguishes distribution and initial-state forms" $ do
            Hit.probability (EqualTo 1) matrix (== 1) mixedInitial
                `shouldBe` 0.125
            Hit.probabilityGivenInitialState (EqualTo 1) matrix (== 1) 0
                `shouldBe` 0.5
            Hit.expectation matrix [1] mixedInitial
                `shouldBe` Right (FiniteExpectation 0.5)
            Return.probability (EqualTo 1) matrix mixedInitial
                `shouldBe` 0.875
            Return.eventualProbability matrix mixedInitial
                `shouldBe` Right 0.875
            Return.expectation matrix mixedInitial
                `shouldBe` Right InfiniteExpectation
            Visit.totalProbability (EqualTo 0) matrix 0 mixedInitial
                `shouldBe` Right 0.75
            Visit.totalExpectation matrix 0 mixedInitial
                `shouldBe` Right (FiniteExpectation 0.5)
            Absorption.probability matrix 1 mixedInitial
                `shouldBe` Right 1
            Absorption.expectation matrix mixedInitial
                `shouldBe` Right (FiniteExpectation 0.5)
