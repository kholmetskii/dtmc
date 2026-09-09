{-# LANGUAGE GADTs #-}

{- |
Module      : Dtmc.Transition.Matrix.Internal
Description : Raw carrier for transition matrices (unsafe underbelly).

Raw carrier behind t'Dtmc.Transition.Matrix.TransitionMatrix': an hmatrix
matrix paired with its lazy support graph and complete classification. The
public smart constructor validates its square shape and canonicalises rows;
this internal module exposes unchecked construction.

The constructor is positional so the public matrix projection cannot act as a
record-update setter and desynchronise the matrix from its cached graph.
-}
module Dtmc.Transition.Matrix.Internal (
    TransitionMatrix (TransitionMatrix),
    unTransitionMatrix,
    tmSupport,
    tmClassification,
    unsafeTransitionMatrix,
    matrixRowAt,
) where

import Data.Maybe (fromMaybe)
import Data.Vector.Storable qualified as Storable
import Dtmc.Analysis.Classification.Internal (
    Classification,
    classificationFromGraph,
 )
import Dtmc.Distribution.Map (
    fromDistribution,
 )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector (DistributionVector),
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.State.Internal (
    stateCardinalityInt,
    stateFromInt,
    stateIndexInt,
 )
import Dtmc.Transition (
    Transition (..),
 )
import Dtmc.Transition.Matrix.Internal.Graph (
    Graph,
    fromAdjacency,
 )
import Numeric.LinearAlgebra qualified as LA

{- | A stored square matrix whose rows and columns follow the canonical order
of its finite state type. Entry @(i,j)@ is the transition probability from
state @i@ to state @j@. 'Dtmc.Transition.Matrix.fromKernel' materialises
already-validated rows, while 'Dtmc.Transition.Matrix.fromRows' applies
tolerant row validation and canonicalisation. The internal constructor and
arithmetic instances do not revalidate.

Each value also carries its support graph and complete typed classification as
/lazy/ arguments, so graph-based analyses on the same value share both the
graph build and public classification results. Construct internal values with
@unsafeTransitionMatrix@ rather than pairing these fields directly.
-}
data TransitionMatrix state where
    -- | Unchecked matrix/cache triple; both caches must match the matrix.
    TransitionMatrix ::
        (FiniteState state) =>
        LA.Matrix Double ->
        Graph ->
        Classification state ->
        TransitionMatrix state

-- Nominal role prevents coercion between distinct state types, including
-- state types with the same cardinality.
type role TransitionMatrix nominal

{- | Return the stored matrix unchanged without forcing the support graph.

Complexity: @O(1)@ time and @O(1)@ space.
-}
unTransitionMatrix ::
    TransitionMatrix state ->
    LA.Matrix Double
unTransitionMatrix (TransitionMatrix matrix _ _) = matrix

{- | Return the lazy support graph, with edge @i -> j@ exactly when the stored
entry is strictly positive. No tolerance is applied: a tiny positive rounding
value creates an edge, while zero or a negative value does not.

The result is shared by later analyses of the same value.

Complexity: @O(1)@ projection time and @O(1)@ projection space. The first
analysis that forces the graph takes @O(n^2)@ time and @O(n^2)@ temporary
space; the resulting graph occupies @O(n + E)@ space for @E@ support edges.
-}
tmSupport :: TransitionMatrix state -> Graph
tmSupport (TransitionMatrix _ support _) = support

{- | Return the lazy complete classification associated with the matrix.
Whole-chain classification queries on the same matrix therefore share their
typed classes and state lists as well as the underlying graph facts.

Complexity: @O(1)@ projection time and space. The first full evaluation takes
the classification cost documented by
'Dtmc.Analysis.Classification.communicatingClasses'; later projections reuse
the retained @O(n)@ classification.
-}
tmClassification :: TransitionMatrix state -> Classification state
tmClassification (TransitionMatrix _ _ classification) = classification
{-# INLINE tmClassification #-}

-- Manual 'Show': 'Graph' has no 'Show', and the derived cache should not
-- appear in the rendering.
instance Show (TransitionMatrix state) where
    showsPrec d p =
        showParen (d > 10) $
            showString "TransitionMatrix "
                . showsPrec 11 (unTransitionMatrix p)

{- | Pair a raw matrix with its lazy support graph and classification. This
performs no row-stochastic, finiteness, or simplex validation; internal
callers must establish the required invariant.

Complexity: @O(1)@ construction time and @O(1)@ construction space. Forcing
the support graph takes @O(n^2)@ time and @O(n^2)@ temporary space; the graph
occupies @O(n + E)@ space for @E@ support edges.
-}
unsafeTransitionMatrix ::
    (FiniteState state) =>
    LA.Matrix Double ->
    TransitionMatrix state
unsafeTransitionMatrix matrix =
    TransitionMatrix matrix support classification
  where
    support = supportGraphOf matrix
    classification = classificationFromGraph toState support
    toState index =
        fromMaybe
            (error "Dtmc.Transition.Matrix.Internal: graph vertex out of bounds")
            (stateFromInt index)
{-# INLINE unsafeTransitionMatrix #-}

{- | Wrap one stored matrix row as a distribution vector without revalidation.
The finite-state index makes the lookup total.

Complexity: excluding 'Dtmc.State.stateIndex', @O(n)@ time and @O(n)@ result
space for state cardinality @n@.
-}
matrixRowAt ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    DistributionVector state
matrixRowAt matrix state = DistributionVector row
  where
    stored = unTransitionMatrix matrix
    row =
        LA.flatten
            ( LA.subMatrix
                (stateIndexInt state, 0)
                (1, LA.cols stored)
                stored
            )

instance (FiniteState state) => Transition (TransitionMatrix state) where
    type TransitionState (TransitionMatrix state) = state

    transitionLaw matrix =
        fromDistribution . matrixRowAt matrix

-- Use strict positivity without tolerance so graph queries reflect the stored
-- matrix exactly; keep construction here so the cache cannot become stale.
supportGraphOf ::
    LA.Matrix Double ->
    Graph
supportGraphOf matrix =
    fromAdjacency dim associations
  where
    dim = LA.rows matrix
    associations = Storable.ifoldr associationFor [] (LA.flatten matrix)
    associationFor offset probability rest =
        let (row, column) = offset `quotRem` dim
         in row
                `seq` column
                `seq` ((row, column), probability > 0) : rest

{- | Matrix multiplication as transition composition: @p '<>' q@ takes a @p@
step followed by a @q@ step. Exact products preserve row-stochasticity and
associativity; 'Double' results are neither revalidated nor exactly
associative.
-}
instance Semigroup (TransitionMatrix state) where
    (<>) ::
        TransitionMatrix state ->
        TransitionMatrix state ->
        TransitionMatrix state
    p@(TransitionMatrix _ _ _) <> q =
        unsafeTransitionMatrix (unTransitionMatrix p LA.<> unTransitionMatrix q)

{- | The identity matrix represents zero transitions and is the unit of the
transition-composition monoid.
-}
instance (FiniteState state) => Monoid (TransitionMatrix state) where
    mempty :: TransitionMatrix state
    mempty = unsafeTransitionMatrix (LA.ident (stateCardinalityInt @state))
