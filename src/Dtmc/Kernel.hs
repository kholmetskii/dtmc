{- |
Module      : Dtmc.Kernel
Description : Shared locally finite transition-kernel capability.

'Transition' is the common finite-horizon boundary for finite matrices and
locally finite countable-state chains. It exposes the validated finite law of
one transition from a supplied state; it does not require global state-space
enumeration.

The instance for 'TransitionMatrix' converts a dense
t'Dtmc.Distribution.Vector.DistributionVector' row to finite support.
This enables shared sparse algorithms, while the specialised finite API keeps
its existing dense implementations for performance and global analyses.
-}
module Dtmc.Kernel (
    Transition (..),
    TransitionKernel,
    transitionKernel,
    deterministicKernel,
) where

import Data.Finite (
    Finite,
 )
import Dtmc.Distribution.Map (
    DistributionMap,
    pointMass,
    toDistributionMap,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    rowAt,
 )
import GHC.TypeNats (
    KnownNat,
 )

{- | A time-homogeneous chain whose transition from any supplied state has
finite support. The complete state space may be finite or infinite.

This capability is sufficient for exact finite-horizon sparse algorithms. It
does not imply that states can be enumerated, so it cannot by itself support
generic classification, stationary, eventual-hitting, or expectation
algorithms.
-}
class Transition transition where
    -- | State type governed by this transition representation.
    type TransitionState transition

    -- | Validated finite-support law of the next state.
    transitionLaw ::
        transition ->
        TransitionState transition ->
        DistributionMap (TransitionState transition)

instance (KnownNat n) => Transition (TransitionMatrix n) where
    type TransitionState (TransitionMatrix n) = Finite n

    transitionLaw matrix =
        toDistributionMap . rowAt matrix

-- | A locally finite kernel over a potentially infinite state type.
newtype TransitionKernel state
    = TransitionKernel (state -> DistributionMap state)

type role TransitionKernel nominal

instance Transition (TransitionKernel state) where
    type TransitionState (TransitionKernel state) = state

    transitionLaw (TransitionKernel kernel) = kernel

{- | Build a kernel from a function returning an already-validated sparse law.
No global traversal is required or attempted.
-}
transitionKernel ::
    (state -> DistributionMap state) ->
    TransitionKernel state
transitionKernel = TransitionKernel

-- | Build the deterministic kernel that maps each state to one successor.
deterministicKernel :: (state -> state) -> TransitionKernel state
deterministicKernel successor =
    transitionKernel (pointMass . successor)
