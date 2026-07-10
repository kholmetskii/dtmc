-- |
-- Module      : Dtmc.Distribution
--
-- A probability distribution over a finite state space: a point of the standard
-- simplex
--
-- > Δ^{n-1} = { x ∈ R^n : x_i ≥ 0, Σ_i x_i = 1 }
--
-- This module owns the predicate. 'Dtmc.TransitionMatrix' does not copy it —
-- it calls 'validateSimplex' row by row (ST227 §2.4: a matrix is row-stochastic
-- iff every row is a distribution). One ε, one check order, no duplication.
--
-- This is the validation boundary: the 'Distribution' constructor is applied
-- here but not re-exported, so 'mkDistribution' is the only door in.
--
-- This module does NOT import 'Dtmc.TransitionMatrix'. The dependency is
-- one-way and asymmetric: Δ(S) knows nothing of kernels; kernels are built out
-- of Δ(S). See docs/DECISIONS.md D3.
module Dtmc.Distribution
  ( Distribution
  , unDistribution
  , SimplexError (..)
  , DistributionError (..)
  , simplexTolerance
  , validateSimplex
  , mkDistribution
  , approxDistributionEq
  ) where

import Data.Bifunctor (first)
import Dtmc.Errors (DistributionError (..), SimplexError (..))
import Dtmc.Internal (Distribution (..))
import GHC.TypeNats (KnownNat)
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Static as S

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
--
-- The tolerance lives next to the proof that justifies it, not in a Config
-- module where it would become a knob (docs/DECISIONS.md D3).
simplexTolerance :: Double
simplexTolerance = 1e-9

-- | The Δ^{n-1} predicate. Builds nothing, mutates nothing.
--
-- Exported because 'Dtmc.TransitionMatrix' delegates to it per row and wraps
-- the error in its own context. The wrapper is per-caller; the predicate is not.
--
-- The check ORDER is contractual: coordinatewise left-to-right first, then the
-- sum. So @[1.5, -0.5]@ yields 'EntryAboveOne' at coordinate 0, not
-- 'NegativeEntry' at coordinate 1. Determinism is what lets rejection tests
-- assert the EXACT error constructor rather than merely @isLeft@.
--
-- Edge case @n = 0@: @entries = []@, @sum [] = 0@, hence @Left (SumOffBy 0.0)@.
-- A zero-dimensional "distribution" is rejected with no special-case code.
--
-- Verified: Dtmc.DistributionSpec
validateSimplex :: KnownNat n => S.R n -> Either SimplexError ()
validateSimplex v =
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

-- | The only door in. Validates but does NOT mutate: no normalising, no
-- clamping. That is why the round-trip property is exact equality rather than
-- approximate.
--
-- Verified: Dtmc.DistributionSpec
mkDistribution :: KnownNat n => S.R n -> Either DistributionError (Distribution n)
mkDistribution v =
  Distribution v <$ first DistributionError (validateSimplex v)

-- | Coordinatewise comparison with an EXPLICIT tolerance.
--
-- The tolerance is a parameter, not 'simplexTolerance': the project has two
-- coexisting regimes six orders of magnitude apart (NUMERICS N6). An explicit
-- argument forces the caller to name which one.
--
-- @and (zipWith ...)@ rather than @LA.maxElement (abs (a - b))@: the latter
-- throws on an empty vector (n = 0).
approxDistributionEq
  :: KnownNat n => Double -> Distribution n -> Distribution n -> Bool
approxDistributionEq tol (Distribution a) (Distribution b) =
  and (zipWith close (LA.toList (S.extract a)) (LA.toList (S.extract b)))
  where
    close x y = abs (x - y) <= tol
