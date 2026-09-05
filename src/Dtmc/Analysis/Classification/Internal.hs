{-# LANGUAGE ExplicitNamespaces #-}

{- |
Module      : Dtmc.Analysis.Classification.Internal
Description : Internal carriers and graph operations for chain classification.

Raw carrier types and solver-oriented graph operations behind
"Dtmc.Analysis.Classification": the per-class summary t'CommClass' and the
whole-chain structural report t'Classification'. This module exposes the
report constructor for trusted internal use; constructing it here may produce
summary fields inconsistent with its communicating classes.
-}
module Dtmc.Analysis.Classification.Internal (
    type CommClass (..),
    type Classification (..),
    backwardReachable,
) where

import Data.Maybe (
    fromMaybe,
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.State.Internal (
    stateFromInt,
    stateIndexInt,
 )
import Dtmc.Transition.Matrix.Internal (
    TransitionMatrix,
    tmSupport,
 )
import Dtmc.Transition.Matrix.Internal.Graph qualified as G
import Numeric.Natural (
    Natural,
 )

{- | Structural facts about one communicating class. For a finite valid DTMC,
a closed class consists of recurrent states.
-}
data CommClass state = CommClass
    { classMembers :: [state]
    -- ^ Member states in ascending order.
    , classPeriod :: Maybe Natural
    -- ^ Shared state period, or 'Nothing' when the class has no cycle.
    , classClosed :: Bool
    -- ^ Whether no positive-probability transition leaves the class.
    }

deriving instance (Eq state) => Eq (CommClass state)

deriving instance (Show state) => Show (CommClass state)

{- | A consistent structural report built by
'Dtmc.Analysis.Classification.classify'. The constructor is exposed here for
trusted internal use; "Dtmc.Analysis.Classification" keeps it hidden so its
summary fields stay aligned with its communicating classes.
-}
data Classification state = Classification
    { classesOf :: [CommClass state]
    -- ^ The communicating classes, ordered by least member.
    , isIrreducible :: Bool
    -- ^ Whether the states form a single (non-empty) communicating class.
    , isAperiodic :: Bool
    -- ^ Whether every class has period @1@ (and there is at least one class).
    , isErgodic :: Bool
    {- ^ Whether the chain is irreducible and aperiodic. For a finite DTMC this
    implies convergence to its unique stationary distribution.
    -}
    , chainPeriod :: Maybe Natural
    {- ^ The period of an irreducible chain (@Just d@), or @Nothing@ for a
    reducible chain, where period is a per-class notion, or when the single
    class has no cycles.
    -}
    , recurrentStatesOf :: [state]
    -- ^ States in closed classes, which are recurrent in a finite chain.
    , transientStatesOf :: [state]
    -- ^ States in non-closed classes, which are transient.
    , absorbingStates :: [state]
    {- ^ Singleton closed classes. For exact stochastic rows these are
    absorbing states with @P(i,i) = 1@; numerically derived or otherwise
    unchecked rows are classified only by strict-positive support.
    -}
    }

type role Classification nominal

deriving instance (Eq state) => Eq (Classification state)

deriving instance (Show state) => Show (Classification state)

toState :: (FiniteState state) => Int -> state
toState index =
    fromMaybe
        (error "Dtmc.Analysis.Classification.Internal: graph vertex out of bounds")
        (stateFromInt index)

toIndex :: (FiniteState state) => state -> Int
toIndex = stateIndexInt

{- | Return states from which an allowed seed is reachable along a support
path containing only states accepted by @allowed@. Disallowed seeds are
ignored; the result is duplicate-free and ordered by state index.

For the complexity bounds, @n@ is the state count, @E@ the support-edge count,
@s@ the number of supplied seeds, and @r@ the number of returned states.

Complexity: excluding @n@ evaluations of @allowed@, 'FiniteState' method
costs, and shared support-graph construction, @O(n + E + s)@ time,
@O(n + E + s)@ temporary space, and @O(r)@ result space. The first reverse
traversal also retains @O(n + E)@ predecessor-cache space.
-}
backwardReachable ::
    (FiniteState state) =>
    TransitionMatrix state ->
    (state -> Bool) ->
    [state] ->
    [state]
backwardReachable p allowed seeds =
    map toState (G.backwardReachable (tmSupport p) (allowed . toState) (map toIndex seeds))
