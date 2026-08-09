{- |
Module      : Dtmc.Transition.Kernel
Description : Locally finite transition kernels over unrestricted state types.

A t'TransitionKernel' represents a transition rule directly as a function from
each state to its validated map-backed next-state distribution. No global
state-space enumeration is required or attempted.
-}
module Dtmc.Transition.Kernel (
    TransitionKernel,
    transitionKernel,
    deterministicKernel,
) where

import Dtmc.Distribution.Map (
    DistributionMap,
    pointMass,
 )
import Dtmc.Transition (
    Transition (..),
 )

-- | A locally finite transition kernel over a potentially infinite state type.
newtype TransitionKernel state
    = TransitionKernel (state -> DistributionMap state)

type role TransitionKernel nominal

instance Transition (TransitionKernel state) where
    type TransitionState (TransitionKernel state) = state

    transitionLaw (TransitionKernel kernel) = kernel

{- | Build a kernel from a function returning an already-validated map-backed
law. No global traversal is required or attempted.
-}
transitionKernel ::
    (state -> DistributionMap state) ->
    TransitionKernel state
transitionKernel = TransitionKernel

-- | Build the deterministic kernel that maps each state to one successor.
deterministicKernel :: (state -> state) -> TransitionKernel state
deterministicKernel successor =
    transitionKernel (pointMass . successor)
