{- |
Module      : Dtmc.TransitionMatrix.Internal
Description : Raw carrier for transition matrices (unsafe underbelly).

Raw carrier behind t'Dtmc.TransitionMatrix.TransitionMatrix': a statically
sized matrix paired with its lazy support graph. The public smart constructor
validates rows; this internal module exposes unchecked construction.

The constructor is positional so the public matrix projection cannot act as a
record-update setter and desynchronise the matrix from its cached graph.
-}
module Dtmc.TransitionMatrix.Internal (
    TransitionMatrix (TransitionMatrix),
    unTransitionMatrix,
    tmSupport,
    unsafeTransitionMatrix,
) where

import Dtmc.Internal.Graph (
    Graph,
    fromAdjacency,
 )
import GHC.TypeNats (
    KnownNat,
    Nat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

{- | A stored @n x n@ matrix with entry @(i,j)@ representing the transition
probability from state @i@ to state @j@. The public constructor applies
tolerant row validation; the internal constructor and arithmetic instances do
not revalidate.

Each value also carries its support graph as a /lazy/ second argument, so any
graph-based analyses on the same value share one build. Construct internal
values with @unsafeTransitionMatrix@ rather than pairing a matrix and graph
directly.
-}
data TransitionMatrix (n :: Nat)
    = -- | Unchecked matrix/cache pair; the graph must match the matrix.
      TransitionMatrix (S.Sq n) Graph

-- Nominal role on @n@, for the same reason as
-- t'Dtmc.Distribution.DistributionVector'.
type role TransitionMatrix nominal

{- | Return the stored matrix unchanged. This is an @O(1)@ projection and does
not force the support graph.
-}
unTransitionMatrix :: TransitionMatrix n -> S.Sq n
unTransitionMatrix (TransitionMatrix matrix _) = matrix

{- | The lazy support graph, with edge @i -> j@ exactly when the stored entry
is strictly positive. No tolerance is applied: a tiny positive rounding value
creates an edge, while zero or a negative value does not.

The projection is @O(1)@. The first analysis that builds adjacency scans all
@n^2@ entries; the result is shared by later analyses of the same value.
-}
tmSupport :: TransitionMatrix n -> Graph
tmSupport (TransitionMatrix _ support) = support

-- Manual 'Show' (not derived): 'Graph' has no 'Show', and the support graph is a
-- derived cache that should not appear in the rendering.
instance (KnownNat n) => Show (TransitionMatrix n) where
    showsPrec d p =
        showParen (d > 10) $
            showString "TransitionMatrix "
                . showsPrec 11 (unTransitionMatrix p)

{- | Pair a raw matrix with its lazy support graph. This performs no
row-stochastic, finiteness, or simplex validation; internal callers must
establish the required invariant.

Construction is @O(1)@ before the graph is forced.
-}
unsafeTransitionMatrix :: (KnownNat n) => S.Sq n -> TransitionMatrix n
unsafeTransitionMatrix matrix =
    TransitionMatrix matrix (supportGraphOf matrix)

-- Use strict positivity without tolerance so graph queries reflect the stored
-- matrix exactly; keep construction here so the cache cannot become stale.
supportGraphOf :: (KnownNat n) => S.Sq n -> Graph
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
instance (KnownNat n) => Semigroup (TransitionMatrix n) where
    (<>) :: TransitionMatrix n -> TransitionMatrix n -> TransitionMatrix n
    p <> q = unsafeTransitionMatrix (unTransitionMatrix p S.<> unTransitionMatrix q)

{- | The identity matrix represents zero transitions and is the unit of the
transition-composition monoid.
-}
instance (KnownNat n) => Monoid (TransitionMatrix n) where
    mempty :: TransitionMatrix n
    mempty = unsafeTransitionMatrix S.eye
