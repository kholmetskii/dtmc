-- |
-- Module      : Dtmc.StochasticMatrix
--
-- The transition matrix: square, row-stochastic.
--
-- The guarantee is two-tiered, and the difference matters:
--
--   * PROVED BY THE COMPILER: dimension, squareness, dimension agreement under
--     multiplication (@Numeric.LinearAlgebra.Static@).
--   * CHECKED ONCE AT A BOUNDARY, THEN CARRIED BY THE TYPE: stochasticity.
--
-- The compiler does NOT prove stochasticity: @S.matrix [1,2,-3,0.5] :: Sq 2@ is
-- a perfectly valid @Sq 2@. And @StochasticMatrix n@ does not mean "the rows
-- sum to 1"; it means "no row deviates from 1 by more than ε, and this was
-- checked once at a single audited boundary" (NUMERICS N2).
--
-- This module does NOT import 'Dtmc.Distribution'.
module Dtmc.StochasticMatrix
  ( StochasticMatrix
  , unStochasticMatrix
  , StochasticError (..)
  , mkStochasticMatrix
  , mulStochasticMatrix
  , approxStochasticMatrixEq
  ) where

import Data.Bifunctor (first)
import Data.Foldable (traverse_)
import Dtmc.Internal (StochasticMatrix (..))
import Dtmc.Simplex (SimplexError, validateSimplexPoint)
import GHC.TypeNats (KnownNat)
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Static as S

-- | Which row violated Δ^{n-1}, and how.
--
-- @data@, not @newtype@: two fields. There are no dimension constructors here —
-- @NonSquareMatrix@ and @MatrixDimensionMismatch@ became unrepresentable.
data StochasticError = InRow Int SimplexError
  deriving (Eq, Show)

-- | Validates stochasticity only. Squareness and size are already in the type.
--
-- Definition (ST227 §2.4): a matrix is row-stochastic iff every row is a
-- distribution. The code restates the definition literally.
--
-- 'traverse_' in 'Either' short-circuits on the first error:
-- @traverse_ f = foldr ((*>) . f) (pure ())@, and @Left e *> _ = Left e@.
-- Applicative, not Monad: the second check does not depend on the first
-- check's result, only on its success.
--
-- Verified: Dtmc.StochasticMatrixSpec
mkStochasticMatrix
  :: KnownNat n => S.Sq n -> Either StochasticError (StochasticMatrix n)
mkStochasticMatrix m =
  StochasticMatrix m <$ traverse_ checkRow (zip [0 ..] (S.toRows m))
  where
    -- S.toRows :: (KnownNat m, KnownNat n) => L m n -> [R n]
    -- Rows are already typed; no detour through the dynamic API.
    checkRow (i, r) = first (InRow i) (validateSimplexPoint r)

-- | The product of two row-stochastic matrices.
--
-- Proof:
--   Let A, B be row-stochastic. Entries are non-negative, so
--   (AB)_ij = Σ_k A_ik B_kj is a sum of products of non-negatives, hence ≥ 0.
--   For every row i:
--     Σ_j (AB)_ij = Σ_j Σ_k A_ik B_kj
--                 = Σ_k A_ik (Σ_j B_kj)
--                 = Σ_k A_ik · 1
--                 = 1.
--   The sums are finite (|S| < ∞), so interchanging the order of summation is
--   unconditional. ∎  (ST227 §2.4; the closure the n-step lemma rests on.)
--
-- The proof pays rent: it licenses a TOTAL signature. No @Either@ — the result
-- cannot fail to be stochastic.
--
-- The result is NOT re-validated through 'mkStochasticMatrix'. See NUMERICS N1.
--
-- Verified: Dtmc.StochasticMatrixSpec.prop_productRowStochastic
mulStochasticMatrix
  :: KnownNat n => StochasticMatrix n -> StochasticMatrix n -> StochasticMatrix n
mulStochasticMatrix (StochasticMatrix a) (StochasticMatrix b) =
  StochasticMatrix (a S.<> b)
  --                  ^^^^^ (<>) :: (KnownNat m, KnownNat k, KnownNat n)
  --                        => L m k -> L k n -> L m n
  -- Dimension agreement is checked by the compiler, not by hmatrix at runtime.
  -- That is why mulStochasticMatrix acquired KnownNat n: the constraint carries
  -- the dimension proof.

-- | Entrywise comparison with an explicit tolerance.
-- See 'Dtmc.Distribution.approxDistributionEq'.
approxStochasticMatrixEq
  :: KnownNat n => Double -> StochasticMatrix n -> StochasticMatrix n -> Bool
approxStochasticMatrixEq tol (StochasticMatrix a) (StochasticMatrix b) =
  and (zipWith close (entries a) (entries b))
  where
    entries = LA.toList . LA.flatten . S.extract
    close x y = abs (x - y) <= tol
