{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TypeApplications #-}

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
import Dtmc.Distribution.Map (pointMass)
import Dtmc.Distribution.Vector (mkDistributionVector)
import Dtmc.Simulation (sample, step)
import Dtmc.State (FiniteState)
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    mkTransitionMatrix,
 )
import GHC.Generics (Generic)
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )
import System.Random.MWC qualified as MWC
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
 )

data NamedSample = FirstSample | SecondSample | ThirdSample
    deriving (Eq, Ord, Show, Generic)

instance FiniteState NamedSample

cyclicThree :: TransitionMatrix (Finite 3)
cyclicThree =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 0
                , 1
                , 0
                , 0
                , 0
                , 1
                , 1
                , 0
                , 0
                ]
            )

namedCyclicThree :: TransitionMatrix NamedSample
namedCyclicThree =
    either (error . show) id $
        mkTransitionMatrix @NamedSample
            (S.matrix [0, 1, 0, 0, 0, 1, 1, 0, 0] :: S.Sq 3)

absorbingTwo :: TransitionMatrix (Finite 2)
absorbingTwo =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 1
                , 0
                , 0.3
                , 0.7
                ]
            )

threeCycleOrbit :: [Finite 3]
threeCycleOrbit = runST $ do
    generator <- MWC.create
    first <- step cyclicThree 0 generator
    second <- step cyclicThree first generator
    third <- step cyclicThree second generator
    pure [first, second, third]

namedThreeCycleOrbit :: [NamedSample]
namedThreeCycleOrbit = runST $ do
    generator <- MWC.create
    first <- step namedCyclicThree FirstSample generator
    second <- step namedCyclicThree first generator
    third <- step namedCyclicThree second generator
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
                mkDistributionVector @(Finite 3) (S.vector [0, 1, 0] :: S.R 3)
    replicateM 20 (sample distribution generator)

namedPointMassSamples :: [NamedSample]
namedPointMassSamples = runST $ do
    generator <- MWC.create
    let distribution =
            either (error . show) id $
                mkDistributionVector @NamedSample
                    (S.vector [0, 1, 0] :: S.R 3)
    replicateM 20 (sample distribution generator)

mapPointMassSamples :: [Natural]
mapPointMassSamples = runST $ do
    generator <- MWC.create
    replicateM 20 (sample (pointMass 7) generator)

spec :: Spec
spec = do
    describe "sample" $ do
        it "samples a dense point mass" $
            pointMassSamples `shouldBe` replicate 20 1

        it "samples a dense point mass as a named state" $
            namedPointMassSamples `shouldBe` replicate 20 SecondSample

        it "samples a map-backed point mass through the same function" $
            mapPointMassSamples `shouldBe` replicate 20 7

    describe "step" $ do
        it "follows a deterministic three-cycle" $
            threeCycleOrbit `shouldBe` [1, 2, 0]

        it "follows a deterministic cycle over named states" $
            namedThreeCycleOrbit
                `shouldBe` [SecondSample, ThirdSample, FirstSample]

        it "never leaves an absorbing state" $
            absorbingSamples `shouldBe` replicate 50 0

        it "runs in ST through PrimMonad" $
            length threeCycleOrbit `shouldBe` 3
