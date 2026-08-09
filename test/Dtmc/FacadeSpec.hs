module Dtmc.FacadeSpec (
    spec,
) where

import Dtmc
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
 )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

finiteInitial :: DistributionVector 2
finiteInitial =
    checked (mkDistributionVector (S.vector [1, 0] :: S.R 2))

finiteTransition :: TransitionMatrix 2
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

spec :: Spec
spec =
    describe "Dtmc facade" $ do
        it "runs a finite vector and matrix through the shared probability API" $
            probabilityAtTime 1 finiteInitial finiteTransition 1
                `shouldBe` 1

        it "runs a map and infinite-state kernel through the same API" $
            probabilityAtTime 3 integerInitial integerTransition 3
                `shouldBe` 1
