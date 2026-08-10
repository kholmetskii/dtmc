{- |
Module      : Dtmc.Transition.Matrix.Internal
Description : Raw carrier for transition matrices (unsafe underbelly).

Raw carrier behind t'Dtmc.Transition.Matrix.TransitionMatrix': a statically
sized matrix paired with its lazy support graph. The public smart constructor
validates rows; this internal module exposes unchecked construction.

The constructor is positional so the public matrix projection cannot act as a
record-update setter and desynchronise the matrix from its cached graph.
-}
module Dtmc.Transition.Matrix.Internal (
    TransitionMatrix (TransitionMatrix),
    unTransitionMatrix,
    tmSupport,
    unsafeTransitionMatrix,
    matrixRowAt,
) where

import Data.Finite (
    getFinite,
 )
import Dtmc.Distribution.Map (
    toDistributionMap,
 )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector (DistributionVector),
 )
import Dtmc.Internal.Graph (
    Graph,
    fromAdjacency,
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
    stateIndex,
 )
import Dtmc.Transition (
    Transition (..),
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

{- | A stored square matrix whose rows and columns follow the canonical order
of its finite state type. Entry @(i,j)@ is the transition probability from
state @i@ to state @j@. The public constructor applies tolerant row validation;
the internal constructor and arithmetic instances do not revalidate.

Each value also carries its support graph as a /lazy/ second argument, so any
graph-based analyses on the same value share one build. Construct internal
values with @unsafeTransitionMatrix@ rather than pairing a matrix and graph
directly.
-}
data TransitionMatrix state
    = -- | Unchecked matrix/cache pair; the graph must match the matrix.
      TransitionMatrix (S.Sq (Cardinality state)) Graph

-- Nominal role prevents coercion between distinct state types, including
-- state types with the same cardinality.
type role TransitionMatrix nominal

{- | Return the stored matrix unchanged. This is an @O(1)@ projection and does
not force the support graph.
-}
unTransitionMatrix ::
    TransitionMatrix state ->
    S.Sq (Cardinality state)
unTransitionMatrix (TransitionMatrix matrix _) = matrix

{- | The lazy support graph, with edge @i -> j@ exactly when the stored entry
is strictly positive. No tolerance is applied: a tiny positive rounding value
creates an edge, while zero or a negative value does not.

The projection is @O(1)@. The first analysis that builds adjacency scans all
@n^2@ entries; the result is shared by later analyses of the same value.
-}
tmSupport :: TransitionMatrix state -> Graph
tmSupport (TransitionMatrix _ support) = support

-- Manual 'Show': 'Graph' has no 'Show', and the derived cache should not
-- appear in the rendering.
instance (FiniteState state) => Show (TransitionMatrix state) where
    showsPrec d p =
        showParen (d > 10) $
            showString "TransitionMatrix "
                . showsPrec 11 (unTransitionMatrix p)

{- | Pair a raw matrix with its lazy support graph. This performs no
row-stochastic, finiteness, or simplex validation; internal callers must
establish the required invariant.

Construction is @O(1)@ before the graph is forced.
-}
unsafeTransitionMatrix ::
    (FiniteState state) =>
    S.Sq (Cardinality state) ->
    TransitionMatrix state
unsafeTransitionMatrix matrix =
    TransitionMatrix matrix (supportGraphOf matrix)

{- | Wrap one stored matrix row as a distribution vector without revalidation.
The finite-state index makes the lookup total.
-}
matrixRowAt ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    DistributionVector state
matrixRowAt matrix state =
    DistributionVector
        ( S.toRows (unTransitionMatrix matrix)
            !! fromIntegral (getFinite (stateIndex state))
        )

instance (FiniteState state) => Transition (TransitionMatrix state) where
    type TransitionState (TransitionMatrix state) = state

    transitionLaw matrix =
        toDistributionMap . matrixRowAt matrix

-- Use strict positivity without tolerance so graph queries reflect the stored
-- matrix exactly; keep construction here so the cache cannot become stale.
supportGraphOf ::
    (KnownNat dimension) =>
    S.Sq dimension ->
    Graph
supportGraphOf matrix =
    fromAdjacency
        dim
        [ ((i, j), entry > 0)
        | (i, row) <- zip [0 ..] rows
        , (j, entry) <- zip [0 ..] row
        ]
  where
    rows = LA.toLists (S.extract matrix)
    dim = length rows

{- | Matrix multiplication as transition composition: @p '<>' q@ takes a @p@
step followed by a @q@ step. Exact products preserve row-stochasticity and
associativity; 'Double' results are neither revalidated nor exactly
associative.
-}
instance (FiniteState state) => Semigroup (TransitionMatrix state) where
    (<>) ::
        TransitionMatrix state ->
        TransitionMatrix state ->
        TransitionMatrix state
    p <> q = unsafeTransitionMatrix (unTransitionMatrix p S.<> unTransitionMatrix q)

{- | The identity matrix represents zero transitions and is the unit of the
transition-composition monoid.
-}
instance (FiniteState state) => Monoid (TransitionMatrix state) where
    mempty :: TransitionMatrix state
    mempty = unsafeTransitionMatrix S.eye
