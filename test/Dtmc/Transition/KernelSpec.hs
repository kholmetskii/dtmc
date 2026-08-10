module Dtmc.Transition.KernelSpec (
    spec,
) where

import Control.Monad.ST (
    runST,
 )
import Dtmc.Distribution qualified as Distribution
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.Dynamics qualified as Dynamics
import Dtmc.Hitting qualified as Hitting
import Dtmc.Probability qualified as Probability
import Dtmc.Simulation qualified as Simulation
import Dtmc.TestSupport (
    approxEq,
    testTolerance,
 )
import Dtmc.Transition.Kernel qualified as Kernel
import System.Random.MWC qualified as MWC
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
    shouldSatisfy,
 )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

simpleRandomWalk :: Kernel.TransitionKernel Integer
simpleRandomWalk =
    Kernel.transitionKernel $ \state ->
        checked (DistributionMap.mkDistributionMap [(state - 1, 0.5), (state + 1, 0.5)])

closeTo :: Double -> Double -> Bool
closeTo = approxEq testTolerance

spec :: Spec
spec = do
    describe "locally finite dynamics" $
        it "evolves an infinite-state random walk without enumerating its state space" $
            Distribution.distributionWeights
                (Dynamics.evolveN 2 (DistributionMap.pointMass 0) simpleRandomWalk)
                `shouldBe` [(-2, 0.25), (0, 0.5), (2, 0.25)]

    describe "finite-horizon probabilities" $
        it "computes transition probabilities on an infinite state type" $ do
            Probability.transitionProbabilityN 2 simpleRandomWalk 0 0
                `shouldSatisfy` closeTo 0.5
            Probability.transitionProbabilityN 3 simpleRandomWalk 0 0
                `shouldBe` 0

    describe "bounded hitting and first return" $ do
        it "uses strict hitting bounds on the infinite random walk" $ do
            Hitting.hittingTimeProbabilityAt simpleRandomWalk (== 2) 0 2
                `shouldSatisfy` closeTo 0.25
            Hitting.hittingTimeProbabilityBefore simpleRandomWalk (== 2) 0 2
                `shouldBe` 0
            Hitting.hittingTimeProbabilityBefore simpleRandomWalk (== 2) 0 3
                `shouldSatisfy` closeTo 0.25

        it "distinguishes return time from time-zero hitting" $ do
            Hitting.hittingTimeProbabilityAt simpleRandomWalk (== 0) 0 0
                `shouldBe` 1
            Hitting.returnTimeProbabilityAt simpleRandomWalk 0 0
                `shouldBe` 0
            Hitting.returnTimeProbabilityAt simpleRandomWalk 0 2
                `shouldSatisfy` closeTo 0.5
            Hitting.returnTimeProbabilityBefore simpleRandomWalk 0 2
                `shouldBe` 0
            Hitting.returnTimeProbabilityBefore simpleRandomWalk 0 3
                `shouldSatisfy` closeTo 0.5

    describe "simulation" $
        it "returns the initial state plus exactly the requested transitions" $
            let trajectory = runST $ do
                    generator <- MWC.create
                    Simulation.simulateN
                        4
                        (Kernel.deterministicKernel (\state -> (state + 1) `mod` (3 :: Int)))
                        0
                        generator
             in trajectory `shouldBe` [0, 1, 2, 0, 1]
