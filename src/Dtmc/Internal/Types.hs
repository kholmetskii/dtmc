-- |
-- Module      : Dtmc.Internal.Types
-- Description : Dimension-indexed carrier types for distributions and matrices.
--
-- Internal module holding the two representation types the rest of the library
-- is built on: thin newtypes over the statically sized vector and matrix of
-- "Numeric.LinearAlgebra.Static", with the type-level 'Nat' @n@ pinning the
-- number of states. The probability invariants (non-negativity, unit sums) are
-- enforced by the smart constructors in "Dtmc.Distribution" and
-- "Dtmc.TransitionMatrix"; this module fixes only the shape. It exposes the raw
-- constructors and so is not part of the public API.
module Dtmc.Internal.Types (
    Distribution (..),
    TransitionMatrix (..),
    transitionMatrix,
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

-- | A probability distribution over the @n@ states @{0 .. n-1}@, stored as a
-- length-@n@ real vector. A well-formed value is a point on the standard
-- simplex (entries in @[0,1]@ summing to one), guaranteed only when built via
-- 'Dtmc.Distribution.mkDistribution'.
newtype Distribution (n :: Nat) = Distribution
    { unDistribution :: S.R n
    }

-- Nominal role: forbids 'Data.Coerce.coerce' from changing @n@. 'S.R' is
-- phantom in its size at runtime, so without this a @Distribution 3@ could be
-- coerced to a @Distribution 5@, breaking the dimension invariant.
type role Distribution nominal

deriving instance (KnownNat n) => Show (Distribution n)

-- | The one-step transition matrix of a chain on @n@ states, stored as an
-- @n*n@ real matrix. A well-formed value is row-stochastic (every row is a
-- 'Distribution'), guaranteed only via 'Dtmc.TransitionMatrix.mkTransitionMatrix'.
-- Entry @(i,j)@ is @P(next = j | now = i)@.
--
-- Each value also carries its support graph as a /lazy/ field, so any number of
-- graph-based analyses on the same value share one graph build, while purely
-- linear-algebraic uses never force it. Build values only through
-- 'transitionMatrix' (or the instances below), which keep 'tmSupport' in step
-- with the matrix.
data TransitionMatrix (n :: Nat) = TransitionMatrix
    { unTransitionMatrix :: S.Sq n
    , tmSupport :: Graph
    -- ^ Lazy support graph: a directed edge @i -> j@ for each @P(i,j) > 0@.
    }

-- Nominal role on @n@, for the same reason as 'Distribution'.
type role TransitionMatrix nominal

-- Manual 'Show' (not derived): 'Graph' has no 'Show', and the support graph is a
-- derived cache that should not appear in the rendering.
instance (KnownNat n) => Show (TransitionMatrix n) where
    showsPrec d p =
        showParen (d > 10) $
            showString "TransitionMatrix {unTransitionMatrix = "
                . shows (unTransitionMatrix p)
                . showString "}"

-- | The single sanctioned builder: pair a raw matrix with its /lazy/ support
-- graph. Every other constructor in the library goes through this, so
-- 'tmSupport' is always the support of 'unTransitionMatrix'.
transitionMatrix :: (KnownNat n) => S.Sq n -> TransitionMatrix n
transitionMatrix matrix =
    TransitionMatrix
        { unTransitionMatrix = matrix
        , tmSupport = supportGraphOfMatrix matrix
        }

-- The support graph of a raw matrix: edge @i -> j@ iff @P(i,j) > 0@. Lives here
-- so the builder and the instances can attach it; "Dtmc.Classification" merely
-- projects 'tmSupport'.
supportGraphOfMatrix :: (KnownNat n) => S.Sq n -> Graph
supportGraphOfMatrix matrix =
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
    p <> q = transitionMatrix (unTransitionMatrix p S.<> unTransitionMatrix q)

-- | The identity matrix is the unit: the zero-step transition that leaves the
-- state unchanged. Together with '<>' this makes @Pow p k@ (via 'mconcat' /
-- 'mtimesDefault') the @k@-step transition matrix.
instance (KnownNat n) => Monoid (TransitionMatrix n) where
    mempty :: TransitionMatrix n
    mempty = transitionMatrix S.eye
