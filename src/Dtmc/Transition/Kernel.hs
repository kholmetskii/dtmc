{- |
Module      : Dtmc.Transition.Kernel
Description : Locally finite transition kernels over unrestricted state types.

A t'TransitionKernel' represents a transition rule directly as a function from
each state to its validated map-backed next-state distribution. No global
state-space enumeration is required or attempted.
-}
module Dtmc.Transition.Kernel (
    TransitionKernel,
    fromLaws,
) where

import Dtmc.Distribution.Map (
    DistributionMap,
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

{- | Construct a kernel from the function that supplies its transition laws.
Each law must already be a validated t'DistributionMap'; no global state-space
traversal is required or attempted.
'Dtmc.Transition.transitionLaw' reads those laws back.

Complexity: @O(1)@ time and @O(1)@ space.
-}
fromLaws ::
    (state -> DistributionMap state) ->
    TransitionKernel state
fromLaws = TransitionKernel
