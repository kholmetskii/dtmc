module Dtmc.SimulationSpec (
    spec,
) where

import Control.Monad (
    replicateM,
 )
import Control.Monad.ST (
    runST,
 )
import Data.Finite (
    Finite,
 )
import Dtmc.Distribution ( mkDistribution )
import Dtmc.Simulation ( sampleFrom, step )
import Dtmc.TransitionMatrix
    ( TransitionMatrix, mkTransitionMatrix )
import Numeric.LinearAlgebra.Static qualified as S
import System.Random.MWC qualified as MWC
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
 )

cyclicThree :: TransitionMatrix 3
cyclicThree =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 0, 1, 0
                , 0, 0, 1
                , 1, 0, 0
                ]
            )

absorbingTwo :: TransitionMatrix 2
absorbingTwo =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 1, 0
                , 0.3, 0.7
                ]
            )

threeCycleOrbit :: [Finite 3]
threeCycleOrbit = runST $ do
    generator <- MWC.create
    first <- step cyclicThree 0 generator
    second <- step cyclicThree first generator
    third <- step cyclicThree second generator
    pure [first, second, third]

absorbingSamples :: [Finite 2]
absorbingSamples = runST $ do
    generator <- MWC.create
    replicateM 50 (step absorbingTwo 0 generator)

pointMassSamples :: [Finite 3]
pointMassSamples = runST $ do
    generator <- MWC.create
    let distribution =
            either (error . show) id $
                mkDistribution (S.vector [0, 1, 0] :: S.R 3)
    replicateM 20 (sampleFrom distribution generator)

spec :: Spec
spec = do
    describe "sampleFrom" $
        it "always samples the support of a point mass" $
            pointMassSamples `shouldBe` replicate 20 1

    describe "step" $ do
        it "follows a deterministic three-cycle" $
            threeCycleOrbit `shouldBe` [1, 2, 0]

        it "never leaves an absorbing state" $
            absorbingSamples `shouldBe` replicate 50 0

        it "runs in ST through PrimMonad" $
            length threeCycleOrbit `shouldBe` 3
