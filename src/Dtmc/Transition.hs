{- |
Module      : Dtmc.Transition
Description : Shared abstraction for locally finite transition rules.

'Transition' captures the operation shared by finite transition matrices and
locally finite kernels: obtaining the validated finite-support law of the next
state from a supplied current state. Concrete representations live in
"Dtmc.Transition.Matrix" and "Dtmc.Transition.Kernel".
-}
module Dtmc.Transition (
    Transition (..),
) where

import Dtmc.Distribution.Map (DistributionMap)

{- | A time-homogeneous transition rule whose law from any supplied state has
finite support. The complete state space may be finite or infinite.

This capability is sufficient for exact finite-horizon map-backed algorithms.
It does not imply that states can be enumerated, so it cannot by itself support
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
