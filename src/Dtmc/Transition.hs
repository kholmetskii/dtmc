{- |
Module      : Dtmc.Transition
Description : Shared abstraction for locally finite transition rules.

'Transition' captures the operation shared by finite transition matrices and
locally finite kernels: obtaining the validated finite-support law of the next
state from a supplied current state. Concrete representations live in
"Dtmc.Transition.Matrix" and "Dtmc.Transition.Kernel".
-}
module Dtmc.Transition (
    Transition (TransitionState, transitionLaw),
) where

import Dtmc.Transition.Internal (
    Transition (TransitionState, transitionLaw),
 )
