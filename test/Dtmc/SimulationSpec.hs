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
import Dtmc.Distribution (
    Distribution (..),
 )
import Dtmc.Distribution.Map (
    pointMass,
    toDistributionMap,
 )
import Dtmc.Distribution.Vector.HMatrix (mkDistributionVector)
import Dtmc.Simulation (
    SimulationError (..),
    sample,
    simulateN,
    step,
 )
import Dtmc.State (FiniteState)
import Dtmc.Transition.Kernel qualified as Kernel
import Dtmc.Transition.Matrix (
    TransitionMatrix,
 )
import Dtmc.Transition.Matrix.HMatrix (
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

newtype UncheckedDistribution
    = UncheckedDistribution [(Int, Double)]

instance Distribution UncheckedDistribution where
    type DistributionState UncheckedDistribution = Int

    probabilityAt (UncheckedDistribution entries) state =
        sum
            [ weight
            | (storedState, weight) <- entries
            , storedState == state
            ]

    distributionWeights (UncheckedDistribution entries) = entries

checkedSimulation :: (Monad m) => m (Either SimulationError value) -> m value
checkedSimulation action = do
    result <- action
    pure (either (error . show) id result)

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
    first <- checkedSimulation (step cyclicThree 0 generator)
    second <- checkedSimulation (step cyclicThree first generator)
    third <- checkedSimulation (step cyclicThree second generator)
    pure [first, second, third]

namedThreeCycleOrbit :: [NamedSample]
namedThreeCycleOrbit = runST $ do
    generator <- MWC.create
    first <- checkedSimulation (step namedCyclicThree FirstSample generator)
    second <- checkedSimulation (step namedCyclicThree first generator)
    third <- checkedSimulation (step namedCyclicThree second generator)
    pure [first, second, third]

absorbingSamples :: [Finite 2]
absorbingSamples = runST $ do
    generator <- MWC.create
    replicateM 50 (checkedSimulation (step absorbingTwo 0 generator))

pointMassSamples :: [Finite 3]
pointMassSamples = runST $ do
    generator <- MWC.create
    let distribution =
            either (error . show) id $
                mkDistributionVector @(Finite 3) (S.vector [0, 1, 0] :: S.R 3)
    replicateM 20 (checkedSimulation (sample distribution generator))

namedPointMassSamples :: [NamedSample]
namedPointMassSamples = runST $ do
    generator <- MWC.create
    let distribution =
            either (error . show) id $
                mkDistributionVector @NamedSample
                    (S.vector [0, 1, 0] :: S.R 3)
    replicateM 20 (checkedSimulation (sample distribution generator))

mapPointMassSamples :: [Natural]
mapPointMassSamples = runST $ do
    generator <- MWC.create
    replicateM 20 (checkedSimulation (sample (pointMass 7) generator))

sampleUnchecked :: [(Int, Double)] -> Either SimulationError Int
sampleUnchecked entries = runST $ do
    generator <- MWC.create
    sample (UncheckedDistribution entries) generator

invalidSampleAndGeneratorState :: (Either SimulationError Int, Bool)
invalidSampleAndGeneratorState = runST $ do
    generator <- MWC.create
    before <- MWC.save generator
    result <- sample (UncheckedDistribution []) generator
    after <- MWC.save generator
    pure (result, before == after)

zeroStepAndGeneratorState :: (Either SimulationError [Finite 3], Bool)
zeroStepAndGeneratorState = runST $ do
    generator <- MWC.create
    before <- MWC.save generator
    result <-
        simulateN
            0
            ( error "zero-step simulation evaluated its kernel" ::
                TransitionMatrix (Finite 3)
            )
            0
            generator
    after <- MWC.save generator
    pure (result, before == after)

emptyKernel :: Kernel.TransitionKernel Int
emptyKernel =
    Kernel.transitionKernel
        (const (toDistributionMap (UncheckedDistribution [])))

spec :: Spec
spec = do
    describe "sample" $ do
        it "samples a dense point mass" $
            pointMassSamples `shouldBe` replicate 20 1

        it "samples a dense point mass as a named state" $
            namedPointMassSamples `shouldBe` replicate 20 SecondSample

        it "samples a map-backed point mass through the same function" $
            mapPointMassSamples `shouldBe` replicate 20 7

        it "repairs a tolerated negative weight" $
            sampleUnchecked [(1, -1e-12), (2, 1)]
                `shouldBe` Right 2

        it "rejects empty stored support" $
            sampleUnchecked []
                `shouldBe` Left EmptySupport

        it "rejects a non-finite weight by index" $ do
            sampleUnchecked [(1, 0 / 0), (2, 1)]
                `shouldBe` Left (NonFiniteWeight 0)
            sampleUnchecked [(1, 1), (2, 1 / 0)]
                `shouldBe` Left (NonFiniteWeight 1)

        it "rejects a weight below the repair tolerance" $
            sampleUnchecked [(1, -1e-6), (2, 1)]
                `shouldBe` Left (NegativeWeight 0 (-1e-6))

        it "rejects a non-positive repaired total" $
            sampleUnchecked [(1, 0), (2, -1e-12)]
                `shouldBe` Left (NonPositiveTotal 0)

        it "rejects overflow in the total" $
            sampleUnchecked
                [(1, 1.7976931348623157e308), (2, 1.7976931348623157e308)]
                `shouldBe` Left NonFiniteTotal

        it "does not advance the generator when validation fails" $ do
            let (result, unchanged) = invalidSampleAndGeneratorState
            result `shouldBe` Left EmptySupport
            unchanged `shouldBe` True

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

        it "returns a transition-law validation failure" $
            runST
                ( do
                    generator <- MWC.create
                    step emptyKernel 0 generator
                )
                `shouldBe` Left EmptySupport

    describe "simulateN" $ do
        it "simulates a finite matrix through the shared interface" $
            let trajectory = runST $ do
                    generator <- MWC.create
                    checkedSimulation (simulateN 4 cyclicThree 0 generator)
             in trajectory `shouldBe` [0, 1, 2, 0, 1]

        it "returns the initial state plus the requested kernel transitions" $
            let trajectory = runST $ do
                    generator <- MWC.create
                    checkedSimulation
                        ( simulateN
                            4
                            ( Kernel.deterministicKernel
                                (\state -> (state + 1) `mod` (3 :: Int))
                            )
                            0
                            generator
                        )
             in trajectory `shouldBe` [0, 1, 2, 0, 1]

        it "does not inspect the kernel or advance the generator at zero steps" $ do
            let (result, unchanged) = zeroStepAndGeneratorState
            result `shouldBe` Right [0]
            unchanged `shouldBe` True

        it "stops at the first invalid transition law" $
            runST
                ( do
                    generator <- MWC.create
                    simulateN 3 emptyKernel 0 generator
                )
                `shouldBe` Left EmptySupport
