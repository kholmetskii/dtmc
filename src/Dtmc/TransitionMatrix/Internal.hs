-- |
-- Module      : Dtmc.TransitionMatrix.Internal
-- Description : Raw carrier for transition matrices (unsafe underbelly).
--
-- The dimension-indexed representation behind
-- t'Dtmc.TransitionMatrix.TransitionMatrix': the @n*n@ real matrix of
-- "Numeric.LinearAlgebra.Static" paired with its support graph, the type-level
-- 'Nat' @n@ pinning the number of states. This module fixes the /shape/ and the
-- @matrix\/support@ pairing; the row-stochastic invariant is enforced by the
-- smart constructor 'Dtmc.TransitionMatrix.mkTransitionMatrix' in the public
-- module. It carries no dependency on "Dtmc.Distribution.Internal" -- the two
-- carriers are independent at the raw level -- so the internal layer is a clean
-- DAG, not a cycle.
--
-- It exposes the raw data constructor and 'unsafeTransitionMatrix' and so is
-- __not__ part of the public API (cabal @other-modules@). "Dtmc.Classification"
-- and "Dtmc.Hitting" import it for the read-only projections; "Dtmc.TransitionMatrix"
-- re-exports only the safe surface.
--
-- The constructor is positional, not a record: a record field @unTransitionMatrix@
-- would export a setter that replaces the matrix while leaving the cached support
-- graph pointing at the old one, silently desyncing the two invariants this type
-- maintains. 'unTransitionMatrix' and 'tmSupport' are plain projections that only
-- read.
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

-- | The one-step transition matrix of a chain on @n@ states, stored as an
-- @n*n@ real matrix. A well-formed value is row-stochastic (every row is a
-- t'Dtmc.Distribution.Distribution'), guaranteed only via
-- 'Dtmc.TransitionMatrix.mkTransitionMatrix'. Entry @(i,j)@ is @P(next = j | now = i)@.
--
-- Each value also carries its support graph as a /lazy/ second argument, so any
-- number of graph-based analyses on the same value share one graph build, while
-- purely linear-algebraic uses never force it. Build values only through
-- 'unsafeTransitionMatrix' (or the instances below), which keep the support
-- graph in step with the matrix.
data TransitionMatrix (n :: Nat) = TransitionMatrix (S.Sq n) Graph

-- Nominal role on @n@, for the same reason as t'Dtmc.Distribution.Distribution'.
type role TransitionMatrix nominal

-- | Recover the underlying statically sized transition matrix. A pure
-- projection, not a record field: it cannot be used to rebuild or mutate a
-- t'TransitionMatrix', so the row-stochastic invariant is safe.
unTransitionMatrix :: TransitionMatrix n -> S.Sq n
unTransitionMatrix (TransitionMatrix matrix _) = matrix

-- | The value's support graph: a directed edge @i -> j@ for each @P(i,j) > 0@.
-- Read-only projection of the /lazy/ cache attached by 'unsafeTransitionMatrix';
-- forced on first use and shared across analyses of the same value.
tmSupport :: TransitionMatrix n -> Graph
tmSupport (TransitionMatrix _ support) = support

-- Manual 'Show' (not derived): 'Graph' has no 'Show', and the support graph is a
-- derived cache that should not appear in the rendering.
instance (KnownNat n) => Show (TransitionMatrix n) where
    showsPrec d p =
        showParen (d > 10) $
            showString "TransitionMatrix "
                . showsPrec 11 (unTransitionMatrix p)

-- | The single sanctioned builder: pair a raw matrix with its /lazy/ support
-- graph. Every other constructor in the library goes through this, so
-- 'tmSupport' is always the support of 'unTransitionMatrix'.
unsafeTransitionMatrix :: (KnownNat n) => S.Sq n -> TransitionMatrix n
unsafeTransitionMatrix matrix =
    TransitionMatrix matrix (supportGraphOf matrix)

-- The support graph of a raw matrix: edge @i -> j@ iff @P(i,j) > 0@. Lives here
-- so the builder and the instances can attach it; "Dtmc.Classification" merely
-- projects 'tmSupport'.
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

-- | Matrix product as composition of steps: @p '<>' q@ is the transition
-- matrix of "do a @p@-step, then a @q@-step". The product of two
-- row-stochastic matrices is again row-stochastic, so the invariant is
-- preserved and this is associative.
instance (KnownNat n) => Semigroup (TransitionMatrix n) where
    (<>) :: TransitionMatrix n -> TransitionMatrix n -> TransitionMatrix n
    p <> q = unsafeTransitionMatrix (unTransitionMatrix p S.<> unTransitionMatrix q)

-- | The identity matrix is the unit: the zero-step transition that leaves the
-- state unchanged. Together with '<>' this makes @Pow p k@ (via 'mconcat' /
-- 'Data.Semigroup.mtimesDefault') the @k@-step transition matrix.
instance (KnownNat n) => Monoid (TransitionMatrix n) where
    mempty :: TransitionMatrix n
    mempty = unsafeTransitionMatrix S.eye
