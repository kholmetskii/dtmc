{-# LANGUAGE ExplicitNamespaces #-}

{- |
Module      : Dtmc.Analysis.Classification.Internal
Description : Internal carriers and construction for chain classification.

Raw carrier types behind "Dtmc.Analysis.Classification": the per-class summary
t'CommClass' and whole-chain structural report t'Classification'. The report
builder is shared with the transition-matrix cache.
-}
module Dtmc.Analysis.Classification.Internal (
    type CommClass (..),
    type Classification (..),
    classificationFromGraph,
) where

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
    , chainPeriodOf :: Maybe Natural
    {- ^ The period of an irreducible chain (@Just d@), or @Nothing@ for a
    reducible chain, where period is a per-class notion, or when the single
    class has no cycles.
    -}
    , recurrentStatesOf :: [state]
    -- ^ States in closed classes, which are recurrent in a finite chain.
    , transientStatesOf :: [state]
    -- ^ States in non-closed classes, which are transient.
    , absorbingStatesOf :: [state]
    {- ^ Singleton closed classes. For exact stochastic rows these are
    absorbing states with @P(i,i) = 1@; numerically derived or otherwise
    unchecked rows are classified only by strict-positive support.
    -}
    }

type role Classification nominal

deriving instance (Eq state) => Eq (Classification state)

deriving instance (Show state) => Show (Classification state)

{- | Build a complete typed classification from an index-to-state conversion
and a shared support graph. Supplying the canonical conversion of a
'Dtmc.State.FiniteState' instance produces the public classification.

The result is deliberately lazy so a transition matrix can carry it without
forcing graph construction. Once forced, all whole-chain projections share
the same classes and state lists.

Complexity: @O(n + E)@ time after the graph's component and period caches are
available, @O(n)@ temporary space, and @O(n)@ retained result space.
-}
classificationFromGraph :: (Int -> state) -> G.Graph -> Classification state
classificationFromGraph toState graph =
    Classification
        { classesOf = communicating
        , isIrreducible = irreducible'
        , isAperiodic = aperiodic'
        , isErgodic = irreducible' && aperiodic'
        , chainPeriodOf = case communicating of
            [communicatingClass] -> classPeriod communicatingClass
            _ -> Nothing
        , recurrentStatesOf =
            concatMap classMembers (filter classClosed communicating)
        , transientStatesOf =
            concatMap classMembers (filter (not . classClosed) communicating)
        , absorbingStatesOf =
            [ state
            | communicatingClass <- communicating
            , classClosed communicatingClass
            , [state] <- [classMembers communicatingClass]
            ]
        }
  where
    components = G.components graph
    communicating =
        [ CommClass
            { classMembers = map toState component
            , classPeriod = G.periodOf graph first
            , classClosed = G.inClosedComponent graph first
            }
        | component@(first : _) <- components
        ]
    irreducible' = case components of
        [component] -> not (null component)
        _ -> False
    aperiodic' =
        not (null components)
            && all ((== Just 1) . G.componentPeriod graph) components
{-# INLINE classificationFromGraph #-}
