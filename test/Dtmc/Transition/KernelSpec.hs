module Dtmc.Transition.KernelSpec (
    spec,
) where

import Dtmc.Distribution qualified as Distribution
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.Transition qualified as Transition
import Dtmc.Transition.Kernel qualified as Kernel
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
 )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

simpleRandomWalk :: Kernel.TransitionKernel Integer
simpleRandomWalk =
    Kernel.transitionKernel $ \state ->
        checked
            (DistributionMap.mkDistributionMap [(state - 1, 0.5), (state + 1, 0.5)])

spec :: Spec
spec =
    describe "TransitionKernel" $ do
        it "preserves source-dependent transition laws" $
            let kernel =
                    Kernel.transitionKernel $ \source ->
                        checked $
                            DistributionMap.mkDistributionMap
                                [(source - 1, 0.4), (source + 1, 0.6 :: Double)]
             in Distribution.distributionWeights (Transition.transitionLaw kernel (10 :: Integer))
                    `shouldBe` [(9, 0.4), (11, 0.6)]

        it "constructs deterministic point-mass transitions" $
            Distribution.distributionWeights
                ( Transition.transitionLaw
                    (Kernel.deterministicKernel (* 2))
                    (6 :: Integer)
                )
                `shouldBe` [(12, 1)]

        it "supports locally finite laws on an infinite state type" $
            Distribution.distributionWeights (Transition.transitionLaw simpleRandomWalk 0)
                `shouldBe` [(-1, 0.5), (1, 0.5)]
