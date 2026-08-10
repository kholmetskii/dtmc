{-# LANGUAGE ScopedTypeVariables #-}

module Dtmc.Transition.TestSupport (
    checked,
    finiteChain,
    finiteInitial,
    finiteCycle,
    asTransitionKernel,
    asDistributionMap,
    kernelChain,
    mapInitial,
    simpleRandomWalk,
    closeTo,
) where

import Data.Finite (
    Finite,
 )
import Dtmc.Distribution qualified as Distribution
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.Distribution.Vector (
    DistributionVector,
    mkDistributionVector,
 )
import Dtmc.Probability qualified as Probability
import Dtmc.State qualified as State
import Dtmc.TestSupport (
    approxEq,
    testTolerance,
 )
import Dtmc.Transition.Kernel qualified as Kernel
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    mkTransitionMatrix,
 )
import Numeric.LinearAlgebra.Static qualified as S

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

finiteInitial :: DistributionVector (Finite 3)
finiteInitial =
    checked (mkDistributionVector (S.vector [0.6, 0.3, 0.1] :: S.R 3))

finiteCycle :: TransitionMatrix (Finite 3)
finiteCycle =
    checked $
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
                ] ::
                S.Sq 3
            )

asTransitionKernel ::
    forall state.
    (State.FiniteState state) =>
    TransitionMatrix state ->
    Kernel.TransitionKernel state
asTransitionKernel matrix =
    Kernel.transitionKernel $ \source ->
        checked $
            DistributionMap.mkDistributionMap
                [ (destination, Probability.transitionProbability matrix source destination)
                | destination <- State.finiteStates
                ]

asDistributionMap ::
    forall state.
    (State.FiniteState state) =>
    DistributionVector state ->
    DistributionMap.DistributionMap state
asDistributionMap distribution =
    checked $
        DistributionMap.mkDistributionMap
            [ (state, Distribution.probabilityAt distribution state)
            | state <- State.finiteStates
            ]

kernelChain :: Kernel.TransitionKernel (Finite 3)
kernelChain = asTransitionKernel finiteChain

mapInitial :: DistributionMap.DistributionMap (Finite 3)
mapInitial = asDistributionMap finiteInitial

simpleRandomWalk :: Kernel.TransitionKernel Integer
simpleRandomWalk =
    Kernel.transitionKernel $ \state ->
        checked (DistributionMap.mkDistributionMap [(state - 1, 0.5), (state + 1, 0.5)])

closeTo :: Double -> Double -> Bool
closeTo = approxEq testTolerance
