{-# LANGUAGE ExplicitNamespaces #-}

{- |
Module      : Dtmc.Classification.Internal
Description : Raw carrier types for chain classification (unsafe underbelly).

Raw carrier types behind "Dtmc.Classification": the per-class summary
'type CommClass', the whole-chain structural report 'type Classification', and
the 'type Irreducible' certificate. This module exposes their constructors so trusted
internal code can build and pattern-match on them directly.

The public "Dtmc.Classification" module re-exports 'type Classification' and
'type Irreducible' /abstractly/ (constructors hidden) and provides the only
validating way to build an 'type Irreducible' witness,
'Dtmc.Classification.witnessIrreducible'. Constructing these values here
bypasses those guarantees, so an 'type Irreducible' built directly is not
certified to wrap an irreducible matrix, and a 'type Classification' built
directly may hold summary fields inconsistent with its communicating classes.
-}
module Dtmc.Classification.Internal (
    type CommClass (..),
    type Classification (..),
    type Irreducible (Irreducible),
    unIrreducible,
) where

import Dtmc.State (
    FiniteState,
 )
import Dtmc.Transition.Matrix.Internal (
    TransitionMatrix,
 )
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

{- | A consistent structural report built by 'Dtmc.Classification.classify'.
The constructor is exposed here for trusted internal use; "Dtmc.Classification"
keeps it hidden so its summary fields stay aligned with its communicating
classes.
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
    {- ^ The period of the chain when it is irreducible (@Just d@); @Nothing@ for
    a reducible chain (where the period is a per-class notion) or when the
    single class has no cycles.
    -}
    , recurrentStatesOf :: [state]
    -- ^ States lying in closed classes -- recurrent, in the finite-chain sense.
    , transientStatesOf :: [state]
    -- ^ States lying in non-closed classes -- transient.
    , absorbingStates :: [state]
    {- ^ Singleton closed classes. For exact stochastic rows these are
    absorbing states with @P(i,i) = 1@; tolerated or unchecked rows are
    classified only by strict-positive support.
    -}
    }

type role Classification nominal

deriving instance (Eq state) => Eq (Classification state)

deriving instance (Show state) => Show (Classification state)

{- | A transition matrix certified as irreducible by
'Dtmc.Classification.witnessIrreducible'. The constructor is exposed here for
trusted internal use but hidden by "Dtmc.Classification", so user code cannot
forge the witness through the public API.
-}
newtype Irreducible state = Irreducible (TransitionMatrix state)

type role Irreducible nominal

deriving instance (FiniteState state) => Show (Irreducible state)

-- | Recover the certified transition matrix in @O(1)@ time.
unIrreducible :: Irreducible state -> TransitionMatrix state
unIrreducible (Irreducible p) = p
