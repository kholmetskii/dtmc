{- |
Module      : Dtmc.Transition.Internal
Description : Transition class and optional dense finite-state backend.

The public "Dtmc.Transition" module exposes only 'TransitionState' and
'transitionLaw'. The remaining method exposes an optional representation
capability, kept internal so optimized library analyses do not enlarge the
public API or burden third-party instances.
-}
module Dtmc.Transition.Internal (
    DenseTransitionBackend (..),
    denseRowBranches,
    Transition (..),
) where

import Data.Vector.Storable qualified as Storable
import Dtmc.Distribution.Map (
    DistributionMap,
 )
import Numeric.LinearAlgebra qualified as LA

{- | Dense finite-state data shared by representation-specific analysis
implementations. State order, indexing, and the matrix must describe the same
canonical coordinate system.

This type is internal: public transition instances need only implement
'transitionLaw'.
-}
data DenseTransitionBackend state = DenseTransitionBackend
    { denseTransitionMatrix :: LA.Matrix Double
    , denseTransitionStates :: [state]
    , denseTransitionIndex :: state -> Int
    }

{- | Test whether one stored row has more than one positive transition. The
scan stops after its second positive entry.
-}
denseRowBranches :: DenseTransitionBackend state -> state -> Bool
denseRowBranches backend state = go 0 0
  where
    stored = denseTransitionMatrix backend
    rowIndex = denseTransitionIndex backend state
    entries =
        LA.flatten
            ( LA.subMatrix
                (rowIndex, 0)
                (1, LA.cols stored)
                stored
            )
    entryCount = Storable.length entries

    go :: Int -> Int -> Bool
    go offset positiveCount
        | positiveCount > 1 = True
        | offset == entryCount = False
        | entries Storable.! offset > 0 = go (offset + 1) (positiveCount + 1)
        | otherwise = go (offset + 1) positiveCount

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

    {- | Return the validated finite-support law of the next state.

    Complexity: implementation-dependent.
    -}
    transitionLaw ::
        transition ->
        TransitionState transition ->
        DistributionMap (TransitionState transition)

    -- Optional representation capability. Existing and third-party instances
    -- inherit the locally finite, representation-independent fallback.
    transitionDenseBackend ::
        transition ->
        Maybe (DenseTransitionBackend (TransitionState transition))
    transitionDenseBackend _ = Nothing

    {-# MINIMAL transitionLaw #-}
