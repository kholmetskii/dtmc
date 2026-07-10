-- |
-- Module      : Dtmc.TransitionMatrix
--
-- The transition matrix of a discrete-time Markov chain: square, row-stochastic.
--
-- Named for the object, not the property. A transition matrix IS a Markov kernel
--
-- > P : S → Δ(S),   i ↦ P_{i·}
--
-- so @TransitionMatrix n@ is morally isomorphic to @Finite n -> Distribution n@.
-- 'rowAt' is that direction of the isomorphism, and it is a method on the
-- matrix, not a separate module. Without it the library has no chain — only a
-- validated square array.
--
-- The dependency on 'Dtmc.Distribution' is one-way: a matrix is a family of
-- distributions, so it knows about them; they know nothing of it. That is
-- @has-a@, not a cycle. See docs/DECISIONS.md D5.
--
-- The guarantee is two-tiered, and the difference matters:
--
--   * PROVED BY THE COMPILER: dimension, squareness, dimension agreement under
--     multiplication (@Numeric.LinearAlgebra.Static@).
--   * CHECKED ONCE AT A BOUNDARY, THEN CARRIED BY THE TYPE: stochasticity.
--
-- The compiler does NOT prove stochasticity: @S.matrix [1,2,-3,0.5] :: Sq 2@ is
-- a perfectly valid @Sq 2@. And @TransitionMatrix n@ does not mean "the rows
-- sum to 1"; it means "no row deviates from 1 by more than ε, and this was
-- checked once at a single audited boundary" (NUMERICS N2).
module Dtmc.TransitionMatrix
  ( TransitionMatrix
  , unTransitionMatrix
  , SimplexError (..)
  , TransitionError (..)
  , mkTransitionMatrix
  , mulTransitionMatrix
  , rowAt
  , approxTransitionMatrixEq
  ) where

import Data.Bifunctor (first)
import Data.Finite (Finite, getFinite)
import Data.Foldable (traverse_)
import Dtmc.Distribution (validateSimplex)
import Dtmc.Errors (SimplexError (..), TransitionError (..))
import Dtmc.Internal (Distribution (..), TransitionMatrix (..))
import GHC.TypeNats (KnownNat)
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Static as S

-- | Validates stochasticity only. Squareness and size are already in the type.
--
-- Definition (ST227 §2.4): a matrix is row-stochastic iff every row is a
-- distribution. The code restates the definition literally — it calls the
-- distribution predicate rather than reimplementing it, which is why there is
-- exactly one ε and one check order in the library.
--
-- 'traverse_' in 'Either' short-circuits on the first error:
-- @traverse_ f = foldr ((*>) . f) (pure ())@, and @Left e *> _ = Left e@.
-- Applicative, not Monad: the second row's check does not depend on the first
-- check's result, only on its success.
--
-- Verified: Dtmc.TransitionMatrixSpec
mkTransitionMatrix
  :: KnownNat n => S.Sq n -> Either TransitionError (TransitionMatrix n)
mkTransitionMatrix m =
  TransitionMatrix m <$ traverse_ checkRow (zip [0 ..] (S.toRows m))
  where
    -- S.toRows :: (KnownNat m, KnownNat n) => L m n -> [R n]
    -- Rows are already typed; no detour through the dynamic API.
    checkRow (i, r) = first (InRow i) (validateSimplex r)

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
-- The result is NOT re-validated through 'mkTransitionMatrix'. See NUMERICS N1.
--
-- Verified: Dtmc.TransitionMatrixSpec.prop_productRowStochastic
mulTransitionMatrix
  :: KnownNat n => TransitionMatrix n -> TransitionMatrix n -> TransitionMatrix n
mulTransitionMatrix (TransitionMatrix a) (TransitionMatrix b) =
  TransitionMatrix (a S.<> b)
  --                  ^^^^^ (<>) :: (KnownNat m, KnownNat k, KnownNat n)
  --                        => L m k -> L k n -> L m n
  -- Dimension agreement is checked by the compiler, not by hmatrix at runtime.
  -- That is why mulTransitionMatrix acquired KnownNat n: the constraint carries
  -- the dimension proof.

-- | Row @i@ as the transition distribution OUT OF state @i@. The kernel view.
--
-- Convention: distributions are ROWS and act on the left (λ ↦ λP). Hence the
-- stationary π is a LEFT eigenvector, πP = π. Confusing rows with columns here
-- breaks all of M5/M6. Caught by the 3-cycle (docs/TESTING.md T2).
--
-- Proof:
--   'mkTransitionMatrix' accepted @p@, therefore every row of @p@ passed
--   'Dtmc.Distribution.validateSimplex', i.e. lies in Δ^{n-1}. The index
--   @i :: Finite n@ satisfies 0 ≤ i < n by construction, and @S.toRows@ on
--   @Sq n@ returns exactly n elements, so @(!!)@ is total. Hence the raw
--   @Distribution@ constructor is applied to a value already satisfying the
--   invariant. ∎
--
-- Note the shape of the signature: NO @Either@. The proof turns a partial
-- operation into a total one, exactly as in 'mulTransitionMatrix'.
-- Re-validation is unnecessary and would be a regression.
--
-- Verified: Dtmc.TransitionMatrixSpec
rowAt :: KnownNat n => TransitionMatrix n -> Finite n -> Distribution n
rowAt p i =
  Distribution (S.toRows (unTransitionMatrix p) !! fromIntegral (getFinite i))
  -- Builds all n rows to get one: O(n²) per call. Deferred until M2 profiling,
  -- see docs/DECISIONS.md D6.

-- | Entrywise comparison with an explicit tolerance.
-- See 'Dtmc.Distribution.approxDistributionEq'.
approxTransitionMatrixEq
  :: KnownNat n => Double -> TransitionMatrix n -> TransitionMatrix n -> Bool
approxTransitionMatrixEq tol (TransitionMatrix a) (TransitionMatrix b) =
  and (zipWith close (entries a) (entries b))
  where
    entries = LA.toList . LA.flatten . S.extract
    close x y = abs (x - y) <= tol
