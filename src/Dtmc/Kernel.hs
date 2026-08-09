{- |
Module      : Dtmc.Kernel
Description : Shared locally finite transition-kernel capability.

'MarkovKernel' is the common finite-horizon boundary for finite matrices and
locally finite countable-state chains. It exposes the validated finite law of
one transition from a supplied state; it does not require global state-space
enumeration.

The instance for 'TransitionMatrix' converts a dense
t'Dtmc.Distribution.DistributionVector' row to finite support.
This enables shared sparse algorithms, while the specialised finite API keeps
its existing dense implementations for performance and global analyses.
-}
module Dtmc.Kernel (
    MarkovKernel (..),
    TransitionKernel,
    transitionKernel,
    transitionsFrom,
    deterministicKernel,
) where

import Data.Finite (
    Finite,
 )
import Dtmc.Distribution (
    Distribution (toSparseDistribution),
    SparseDistribution,
    pointMass,
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
class MarkovKernel kernel where
    -- | State type governed by this kernel representation.
    type KernelState kernel

    -- | Validated finite-support law of the next state.
    transitionLaw :: kernel -> KernelState kernel -> SparseDistribution (KernelState kernel)

instance (KnownNat n) => MarkovKernel (TransitionMatrix n) where
    type KernelState (TransitionMatrix n) = Finite n

    transitionLaw matrix =
        toSparseDistribution . rowAt matrix

-- | A locally finite kernel over a potentially infinite state type.
newtype TransitionKernel state
    = TransitionKernel (state -> SparseDistribution state)

type role TransitionKernel nominal

instance MarkovKernel (TransitionKernel state) where
    type KernelState (TransitionKernel state) = state

    transitionLaw = transitionsFrom

{- | Build a kernel from a function returning an already-validated sparse law.
No global traversal is required or attempted.
-}
transitionKernel ::
    (state -> SparseDistribution state) ->
    TransitionKernel state
transitionKernel = TransitionKernel

-- | Obtain the next-state law from one state.
transitionsFrom :: TransitionKernel state -> state -> SparseDistribution state
transitionsFrom (TransitionKernel kernel) = kernel

-- | Build the deterministic kernel that maps each state to one successor.
deterministicKernel :: (state -> state) -> TransitionKernel state
deterministicKernel successor =
    transitionKernel (pointMass . successor)
