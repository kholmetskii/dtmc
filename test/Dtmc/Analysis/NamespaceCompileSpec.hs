{-# LANGUAGE DataKinds #-}

module Dtmc.Analysis.NamespaceCompileSpec (
    spec,
) where

import Data.Finite (
    Finite,
 )
import Data.List.NonEmpty (
    NonEmpty ((:|)),
 )
import Dtmc qualified as Core
import Dtmc.Analysis.Event (
    DiscreteEvent (..),
    matchesDiscreteEvent,
 )
import Dtmc.Analysis.FiniteTime qualified as FT
import Dtmc.Analysis.HittingTime qualified as Hit
import Dtmc.Analysis.ReturnTime qualified as Return
import Dtmc.Analysis.VisitCount qualified as Visit
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
 )

matrix :: Core.TransitionMatrix (Finite 2)
matrix =
    checked
        ( Core.mkTransitionMatrix
            (S.matrix [0.5, 0.5, 0, 1] :: S.Sq 2)
        )

initial :: Core.DistributionVector (Finite 2)
initial =
    checked
        (Core.mkDistributionVector (S.vector [1, 0] :: S.R 2))

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

spec :: Spec
spec = do
    describe "qualified analysis namespaces" $ do
        it "coexist under the documented aliases" $ do
            FT.stepProbability matrix 0 1 `shouldBe` 0.5
            FT.nStepProbability 2 matrix 0 1 `shouldBe` 0.75
            FT.stateProbability 1 initial matrix 1 `shouldBe` 0.5
            FT.pathProbability initial matrix (0 :| [1]) `shouldBe` 0.5
            FT.observationProbability initial matrix [FT.At 1 1]
                `shouldBe` 0.5
            FT.conditionalObservationProbability
                initial
                matrix
                [FT.At 1 1]
                [FT.At 0 0]
                `shouldBe` Right 0.5
            Hit.probability (EqualTo 1) matrix (== 1) 0 `shouldBe` 0.5
            length
                [ Hit.probability (AtMost 1) matrix (== 1) 0 `seq` ()
                , Hit.probabilityByState (AtMost 1) matrix [1] `seq` ()
                , Hit.eventualProbability matrix [1] 0 `seq` ()
                , Hit.eventualProbabilityByState matrix [1] `seq` ()
                , Hit.raceProbability matrix [1] [] 0 `seq` ()
                , Hit.raceProbabilityByState matrix [1] [] `seq` ()
                , Hit.expectation matrix [1] 0 `seq` ()
                , Hit.expectationByState matrix [1] `seq` ()
                ]
                `shouldBe` 8
            Return.probability (EqualTo 1) matrix 0 `shouldBe` 0.5
            length
                [ Return.probability (AtMost 1) matrix 0 `seq` ()
                , Return.probabilityByState (AtMost 1) matrix `seq` ()
                , Return.eventualProbability matrix 0 `seq` ()
                , Return.eventualProbabilityByState matrix `seq` ()
                , Return.expectation matrix 0 `seq` ()
                , Return.expectationByState matrix `seq` ()
                ]
                `shouldBe` 6
            Visit.boundedProbability 2 (EqualTo 1) initial matrix (== 1)
                `shouldBe` 0.5
            length
                [ Visit.totalProbability (EqualTo 1) matrix 1 0 `seq` ()
                , Visit.totalProbabilityByState (AtMost 1) matrix 1 `seq` ()
                , Visit.infiniteProbability matrix 1 0 `seq` ()
                , Visit.infiniteProbabilityByState matrix 1 `seq` ()
                , Visit.totalExpectation matrix 1 0 `seq` ()
                , Visit.totalExpectationByState matrix 1 `seq` ()
                , Visit.boundedLaw 2 initial matrix (== 1) `seq` ()
                , Visit.boundedProbability 2 (AtMost 1) initial matrix (== 1) `seq` ()
                , Visit.boundedExpectation 2 initial matrix (== 1) `seq` ()
                ]
                `shouldBe` 9
            matchesDiscreteEvent (AtMost 1) 1 `shouldBe` True
