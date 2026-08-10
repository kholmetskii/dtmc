module Dtmc.Transition.KernelSpec (
    spec,
) where

import Dtmc.Distribution qualified as Distribution
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.Transition qualified as Transition
import Dtmc.Transition.Kernel qualified as Kernel
import Dtmc.Transition.TestSupport (
    checked,
    simpleRandomWalk,
 )
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
 )

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
