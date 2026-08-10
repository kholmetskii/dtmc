{-# LANGUAGE DeriveGeneric #-}

module Dtmc.FacadeSpec (
    spec,
) where

import Data.Finite (
    Finite,
 )
import Dtmc
import GHC.Generics (
    Generic,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
 )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

data CafeState
    = Thinking
    | Menu
    | Drink
    | Food
    | PlainWaffle
    | ChocolateWaffle
    | Leave
    deriving (Eq, Ord, Show, Generic)

instance FiniteState CafeState

finiteInitial :: DistributionVector (Finite 2)
finiteInitial =
    checked (mkDistributionVector (S.vector [1, 0] :: S.R 2))

finiteTransition :: TransitionMatrix (Finite 2)
finiteTransition =
    checked
        ( mkTransitionMatrix
            ( S.matrix
                [ 0
                , 1
                , 1
                , 0
                ] ::
                S.Sq 2
            )
        )

integerInitial :: DistributionMap Integer
integerInitial = pointMass 0

integerTransition :: TransitionKernel Integer
integerTransition = deterministicKernel (+ 1)

cafeInitial :: DistributionVector CafeState
cafeInitial =
    checked
        ( mkDistributionVector
            (S.vector [1, 0, 0, 0, 0, 0, 0] :: S.R 7)
        )

cafeTransition :: TransitionMatrix CafeState
cafeTransition =
    checked
        ( mkTransitionMatrix
            ( S.matrix
                [ 0
                , 1 / 5
                , 0
                , 1 / 5
                , 1 / 5
                , 1 / 5
                , 1 / 5
                , 1 / 5
                , 0
                , 2 / 5
                , 0
                , 2 / 5
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1 / 2
                , 0
                , 0
                , 1 / 2
                , 1 / 2
                , 0
                , 0
                , 0
                , 0
                , 1 / 2
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
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1
                ] ::
                S.Sq 7
            )
        )

spec :: Spec
spec =
    describe "Dtmc facade" $ do
        it "runs a finite vector and matrix through the shared probability API" $
            probabilityAtTime 1 finiteInitial finiteTransition 1
                `shouldBe` 1

        it "runs a map and infinite-state kernel through the same API" $
            probabilityAtTime 3 integerInitial integerTransition 3
                `shouldBe` 1

        it "runs the seven-state cafe analysis entirely with named states" $ do
            probabilityAt cafeInitial Thinking `shouldBe` 1
            absorbingStates (classify cafeTransition) `shouldBe` [Leave]
            abs
                ( hittingBeforeProbability
                    cafeTransition
                    [PlainWaffle, ChocolateWaffle]
                    [Drink, Leave]
                    Thinking
                    - 29 / 43
                )
                < 1e-12
                `shouldBe` True
