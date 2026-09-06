{-# LANGUAGE TypeApplications #-}

module Dtmc.Transition.InterfaceSpec (
    spec,
) where

import Data.Finite (
    Finite,
    finites,
 )
import Dtmc.Analysis.FiniteTime (
    stepProbability,
 )
import Dtmc.Distribution qualified as Distribution
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.TestSupport (
    approxEq,
    chunksOf,
    genTransitionRows,
 )
import Dtmc.Transition qualified as Transition
import Dtmc.Transition.Kernel qualified as Kernel
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    TransitionMatrixError,
    fromRows,
 )
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
 )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

finiteChain :: TransitionMatrix (Finite 3)
finiteChain =
    checked $
        fromRows
            ( chunksOf
                3
                [ 0.5
                , 0.5
                , 0
                , 0
                , 0.2
                , 0.8
                , 1
                , 0
                , 0
                ]
            )

asTransitionKernel ::
    TransitionMatrix (Finite 3) ->
    Kernel.TransitionKernel (Finite 3)
asTransitionKernel matrix =
    Kernel.fromLaws $ \source ->
        checked $
            DistributionMap.fromList
                [ (destination, stepProbability matrix source destination)
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
                    Kernel.fromLaws $ \source ->
                        checked $
                            DistributionMap.fromList
                                [(source, 0.25), (source + 1, 0.75 :: Double)]
             in Distribution.distributionWeights (Transition.transitionLaw kernel (4 :: Int))
                    `shouldBe` [(4, 0.25), (5, 0.75)]

        it "exposes deterministic kernels as point-mass laws" $
            Distribution.distributionWeights
                ( Transition.transitionLaw
                    (Kernel.fromLaws (DistributionMap.pointMass . (+ 1)))
                    (4 :: Int)
                )
                `shouldBe` [(5, 1)]

        prop "gives matrices and equivalent kernels approximately equal laws" $
            forAll (genTransitionRows 3) $ \rawMatrix ->
                case fromRows rawMatrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 3)) of
                    Left problem -> counterexample (show problem) False
                    Right matrix ->
                        let kernel = asTransitionKernel matrix
                         in conjoin
                                [ let matrixLaw =
                                        Transition.transitionLaw matrix source
                                      kernelLaw =
                                        Transition.transitionLaw kernel source
                                   in counterexample ("source: " <> show source) $
                                        Distribution.support matrixLaw
                                            == Distribution.support kernelLaw
                                            && and
                                                [ approxEq
                                                    1e-12
                                                    ( Distribution.probabilityAt
                                                        matrixLaw
                                                        destination
                                                    )
                                                    ( Distribution.probabilityAt
                                                        kernelLaw
                                                        destination
                                                    )
                                                | destination <- finites :: [Finite 3]
                                                ]
                                | source <- finites :: [Finite 3]
                                ]
