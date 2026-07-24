-- |
-- Module      : Dtmc.Distribution.Internal
-- Description : Raw carrier for probability distributions (unsafe underbelly).
--
-- The dimension-indexed representation behind t'Dtmc.Distribution.Distribution':
-- a thin newtype over the statically sized vector of
-- "Numeric.LinearAlgebra.Static", with the type-level 'Nat' @n@ pinning the
-- number of states. This module fixes only the /shape/; the simplex invariant
-- (non-negative entries summing to one) is enforced by the smart constructor
-- 'Dtmc.Distribution.mkDistribution' in the public module.
--
-- It exposes the raw data constructor and so is __not__ part of the public API
-- (it lives in the cabal @other-modules@). The public "Dtmc.Distribution"
-- re-exports only the safe surface -- the type, the smart constructor, and the
-- read-only projection 'unDistribution' -- while sibling library modules
-- (@rowAt@, @evolve@) import this module to build values directly, which they
-- may do because they can argue the invariant holds.
--
-- The constructor is deliberately __not__ a record field: a field label would
-- export a record-update setter (@d { unDistribution = v }@) needing only the
-- label in scope, not the hidden constructor, letting callers forge a value off
-- the simplex. 'unDistribution' is a plain projection that only reads.
module Dtmc.Distribution.Internal (
    Distribution (Distribution),
    unDistribution,
) where

import GHC.TypeNats (
    KnownNat,
    Nat,
 )
import Numeric.LinearAlgebra.Static qualified as S

-- | A probability distribution over the @n@ states @{0 .. n-1}@, stored as a
-- length-@n@ real vector. A well-formed value is a point on the standard
-- simplex (entries in @[0,1]@ summing to one), guaranteed only when built via
-- 'Dtmc.Distribution.mkDistribution'.
newtype Distribution (n :: Nat) = Distribution (S.R n)

-- Nominal role: forbids 'Data.Coerce.coerce' from changing @n@. 'S.R' is
-- phantom in its size at runtime, so without this a @Distribution 3@ could be
-- coerced to a @Distribution 5@, breaking the dimension invariant.
type role Distribution nominal

deriving instance (KnownNat n) => Show (Distribution n)

-- | Recover the underlying statically sized probability vector. A pure
-- projection (not a record field): it observes a t'Distribution' without
-- offering any way to rebuild one, so the simplex invariant cannot be bypassed
-- through it.
unDistribution :: Distribution n -> S.R n
unDistribution (Distribution vector) = vector
