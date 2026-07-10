module Dtmc.StochasticMatrixSpec (spec) where

import Dtmc.Generators (bumpSmallestInRow0, genStochasticSq, onLists, setEntry00)
import Dtmc.Simplex (SimplexError (..))
import Dtmc.StochasticMatrix
import qualified Numeric.LinearAlgebra.Static as S
import Test.Hspec
import Test.Hspec.QuickCheck (prop)
import Test.QuickCheck

-- DELETED during the migration to Static, and that is the measurable result:
--
--   it "rejects a non-square matrix"
--   it "rejects a square matrix with the wrong type-level dimension"
--
-- Neither compiles any more: dimension errors became unrepresentable.
--
-- n = 3, not 2: at n = 2 every permutation matrix is symmetric, so rows and
-- columns are indistinguishable (docs/TESTING.md T2).

spec :: Spec
spec = do
  prop "round-trips exactly" $
    forAll (genStochasticSq @3) $ \m ->
      case mkStochasticMatrix m of
        Right sm -> property (S.extract (unStochasticMatrix sm) == S.extract m)
        Left e   -> counterexample ("generator produced a rejected matrix: " <> show e) False

  prop "rejects a row-sum bumped by 1e-6, naming the row" $
    forAll (genStochasticSq @3) $ \m ->
      case mkStochasticMatrix (onLists (bumpSmallestInRow0 1e-6) m) of
        Left (InRow 0 (SumOffBy _)) -> property True
        Left e  -> counterexample ("expected InRow 0 (SumOffBy _), got " <> show e) False
        Right _ -> counterexample "expected rejection, got Right" False

  prop "rejects a negative entry, naming row and column" $
    forAll (genStochasticSq @3) $ \m ->
      case mkStochasticMatrix (onLists (setEntry00 (-1e-6)) m) of
        Left (InRow 0 (NegativeEntry 0 _)) -> property True
        Left e  -> counterexample ("expected InRow 0 (NegativeEntry 0 _), got " <> show e) False
        Right _ -> counterexample "expected rejection, got Right" False

  -- Re-validation is APPROPRIATE here. The implementation must not re-validate
  -- (the theorem is exact); the test must. That is "verified" paired to "proved".
  prop "prop_productRowStochastic: closure under multiplication" $
    forAll ((,) <$> genStochasticSq @3 <*> genStochasticSq @3) $ \(a, b) ->
      case (mkStochasticMatrix a, mkStochasticMatrix b) of
        (Right sa, Right sb) ->
          case mkStochasticMatrix (unStochasticMatrix (mulStochasticMatrix sa sb)) of
            Right _ -> property True
            Left e  -> counterexample ("product rejected: " <> show e) False
        _ -> counterexample "generator produced a rejected matrix" False
