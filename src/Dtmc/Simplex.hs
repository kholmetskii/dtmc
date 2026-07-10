-- |
-- Module      : Dtmc.Simplex
--
-- Membership in the standard simplex
--
-- > Δ^{n-1} = { x ∈ R^n : x_i ≥ 0, Σ_i x_i = 1 }
--
-- This is the library's ONLY probability predicate. Both 'Dtmc.Distribution'
-- and 'Dtmc.StochasticMatrix' go through it: a row of a stochastic matrix
-- lies in Δ^{n-1} in exactly the sense a distribution does (ST227 §2.4:
-- row-stochastic ⟺ every row is a distribution).
--
-- Hence a single copy of 'simplexTolerance' and a single check order.
-- Duplicating the predicate would mean two epsilons that drift apart, and
-- rejection tests exercising the wrong copy.
--
-- The error carries NO row index: the predicate does not know who called it.
-- Callers add context ('Dtmc.StochasticMatrix.InRow').
module Dtmc.Simplex
  ( SimplexError (..)
  , simplexTolerance
  , validateSimplexPoint
  ) where

import GHC.TypeNats (KnownNat)
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Static as S

-- | Why a point fails to lie in the simplex. 'Int' is the coordinate index.
data SimplexError
  = NegativeEntry Int Double
  | EntryAboveOne Int Double
  | SumOffBy Double
  deriving (Eq, Show)

-- | Validation tolerance, absolute (the target is exactly 1, so absolute and
-- relative coincide).
--
-- Bounded from both sides, see NUMERICS N2:
--
-- > generator normalisation error  <  ε  <<  meaningful invalidity
-- >         ~ n·u,  u ≈ 1.1e-16                (what a model would notice)
--
-- The lower bound is forced by the round-trip property: the generator
-- normalises rows in Double, so their sums land in 1 ± O(n·u). At n ≤ 10³ that
-- is ≈ 1e-13, four orders of margin below 1e-9.
simplexTolerance :: Double
simplexTolerance = 1e-9

-- | Checks membership in Δ^{n-1}. Builds nothing, mutates nothing.
--
-- The check ORDER is contractual: coordinatewise left-to-right first, then the
-- sum. So @[1.5, -0.5]@ yields 'EntryAboveOne' at coordinate 0, not
-- 'NegativeEntry' at coordinate 1. Determinism is what lets rejection tests
-- assert the EXACT error constructor rather than merely @isLeft@.
--
-- Edge case @n = 0@: @entries = []@, @sum [] = 0@, hence @Left (SumOffBy 0.0)@.
-- A zero-dimensional "distribution" is rejected with no special-case code.
--
-- Verified: Dtmc.SimplexSpec
validateSimplexPoint :: KnownNat n => S.R n -> Either SimplexError ()
validateSimplexPoint v =
  case firstBadEntry 0 entries of
    Just err -> Left err
    Nothing
      | abs (total - 1.0) <= simplexTolerance -> Right ()
      | otherwise                             -> Left (SumOffBy total)
  where
    entries = LA.toList (S.extract v)
    -- Naive summation: error O(n·u). Compensated summation would give O(u).
    -- At n ≤ 10³ the margin below ε is four orders. See NUMERICS N3.
    total = sum entries

-- | The first coordinate violating the pointwise bounds.
--
-- @x < -ε@, NOT @x < 0@: floating-point normalisation readily produces -1e-17
-- instead of zero, and a strict check would reject mathematically valid points.
-- The sliver @[-ε, 0)@ is rounding noise; anything below is a genuine negative.
--
-- 'EntryAboveOne' nearly follows from the other two checks: if all @x_i ≥ -ε@
-- and @|Σ - 1| ≤ ε@, then @x_j ≤ 1 + nε@. "Nearly", because a coordinate can
-- land in @(1+ε, 1+nε]@ — e.g. @[1+2ε, -ε, -ε]@ at n = 3. The check costs one
-- comparison, names the exact coordinate, and does not rely on that inference.
firstBadEntry :: Int -> [Double] -> Maybe SimplexError
firstBadEntry _ [] = Nothing
firstBadEntry i (x : xs)
  | x < negate simplexTolerance = Just (NegativeEntry i x)
  | x > 1.0 + simplexTolerance  = Just (EntryAboveOne i x)
  | otherwise                   = firstBadEntry (i + 1) xs
