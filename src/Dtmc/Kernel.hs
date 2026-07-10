-- |
-- Module      : Dtmc.Kernel
--
-- Where a matrix becomes a Markov chain.
--
-- A stochastic matrix is not "a square array of numbers with a property".
-- It is a Markov kernel
--
-- > P : S → Δ(S),   i ↦ P_{i·}
--
-- so @StochasticMatrix n@ is morally isomorphic to @Finite n -> Distribution n@.
-- 'rowAt' is that direction of the isomorphism. Without it the library has no
-- chain — only a validated matrix.
--
-- 'Dtmc.Distribution' and 'Dtmc.StochasticMatrix' know nothing of each other.
-- The operation connecting them lives HERE. See docs/DECISIONS.md D5.
module Dtmc.Kernel
  ( rowAt
  ) where

import Data.Finite (Finite, getFinite)
import Dtmc.Internal (Distribution (..), StochasticMatrix (..))
import GHC.TypeNats (KnownNat)
import qualified Numeric.LinearAlgebra.Static as S

-- | Row @i@ as the transition distribution OUT OF state @i@.
--
-- Convention: distributions are ROWS and act on the left (λ ↦ λP). Hence the
-- stationary π is a LEFT eigenvector, πP = π. Confusing rows with columns here
-- breaks all of M5/M6. Caught by the 3-cycle (docs/TESTING.md T2).
--
-- Proof:
--   'Dtmc.StochasticMatrix.mkStochasticMatrix' accepted @p@, therefore every
--   row of @p@ passed 'Dtmc.Simplex.validateSimplexPoint', i.e. lies in
--   Δ^{n-1}. The index @i :: Finite n@ satisfies 0 ≤ i < n by construction, and
--   @S.toRows@ on @Sq n@ returns exactly n elements, so @(!!)@ is total.
--   Hence the raw @Distribution@ constructor is applied to a value already
--   satisfying the invariant. ∎
--
-- Note the shape of the signature: NO @Either@. The proof turns a partial
-- operation into a total one, exactly as in
-- 'Dtmc.StochasticMatrix.mulStochasticMatrix'. Re-validation is unnecessary and
-- would be a regression.
--
-- Verified: Dtmc.KernelSpec
rowAt :: KnownNat n => StochasticMatrix n -> Finite n -> Distribution n
rowAt p i =
  Distribution (S.toRows (unStochasticMatrix p) !! fromIntegral (getFinite i))
  -- Builds all n rows to get one: O(n²) per call. Deferred until M2 profiling,
  -- see docs/DECISIONS.md D6.
