{-# LANGUAGE TypeApplications #-}

module Dtmc.Transition.InterfaceSpec (
    spec,
) where

import Data.Finite (
    Finite,
    finites,
 )
import Dtmc.Analysis.FiniteTime (
    transitionProbability,
 )
import Dtmc.Distribution qualified as Distribution
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.TestSupport (
    genTransitionMatrix,
 )
import Dtmc.Transition qualified as Transition
import Dtmc.Transition.Kernel qualified as Kernel
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    mkTransitionMatrix,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
 )
import Test.Hspec.QuickCheck (
    prop,
 )
import Test.QuickCheck (
    conjoin,
    counterexample,
    forAll,
    (===),
 )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

finiteChain :: TransitionMatrix (Finite 3)
finiteChain =
    checked $
        mkTransitionMatrix
            ( S.matrix
                [ 0.5
                , 0.5
                , 0
                , 0
                , 0.2
                , 0.8
                , 1
                , 0
                , 0
                ] ::
                S.Sq 3
            )

asTransitionKernel ::
    TransitionMatrix (Finite 3) ->
    Kernel.TransitionKernel (Finite 3)
asTransitionKernel matrix =
    Kernel.transitionKernel $ \source ->
        checked $
            DistributionMap.mkDistributionMap
                [ (destination, transitionProbability matrix source destination)
                | destination <- finites
                ]

spec :: Spec
spec =
    describe "Transition interface" $ do
        it "exposes a matrix row as a finite-support transition law" $
            Distribution.distributionWeights (Transition.transitionLaw finiteChain 1)
                `shouldBe` [(1, 0.2), (2, 0.8)]

        it "exposes a source-dependent kernel through the same operation" $
            let kernel =
                    Kernel.transitionKernel $ \source ->
                        checked $
                            DistributionMap.mkDistributionMap
                                [(source, 0.25), (source + 1, 0.75 :: Double)]
             in Distribution.distributionWeights (Transition.transitionLaw kernel (4 :: Int))
                    `shouldBe` [(4, 0.25), (5, 0.75)]

        it "exposes deterministic kernels as point-mass laws" $
            Distribution.distributionWeights
                (Transition.transitionLaw (Kernel.deterministicKernel (+ 1)) (4 :: Int))
                `shouldBe` [(5, 1)]

        prop "gives matrices and equivalent kernels the same transition laws" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left problem -> counterexample (show problem) False
                    Right matrix ->
                        let kernel = asTransitionKernel matrix
                         in conjoin
                                [ Distribution.distributionWeights
                                    (Transition.transitionLaw matrix source)
                                    === Distribution.distributionWeights
                                        (Transition.transitionLaw kernel source)
                                | source <- finites :: [Finite 3]
                                ]
