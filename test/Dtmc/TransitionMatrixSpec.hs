module Dtmc.TransitionMatrixSpec (spec) where

import Dtmc.Distribution (mkDistribution, unDistribution)
import Dtmc.Generators (bumpSmallestInRow0, genStochasticSq, onLists, setEntry00)
import Dtmc.TransitionMatrix
    ( TransitionError(InRow),
      SimplexError(NegativeEntry, SumOffBy),
      TransitionMatrix(..),
      mkTransitionMatrix,
      mulTransitionMatrix,
      rowAt )
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Static as S
import Test.Hspec ( Spec, it, shouldBe )
import Test.Hspec.QuickCheck (prop)
import Test.QuickCheck
    ( Testable(property), forAll, counterexample, conjoin )

-- DELETED during the migration to Static, and that is the measurable result:
--
--   it "rejects a non-square matrix"
--   it "rejects a square matrix with the wrong type-level dimension"
--
-- Neither compiles any more: dimension errors became unrepresentable.
--
-- n = 3, not 2: at n = 2 every permutation matrix is symmetric, so rows and
-- columns are indistinguishable (docs/TESTING.md T2).

-- | The 3-cycle: row 0 = e₁, row 1 = e₂, row 2 = e₀.
--
--   under P  : 0 → 1 → 2 → 0
--   under Pᵀ : 0 → 2 → 1 → 0
--
-- Asymmetric, so it catches confused rows/columns. The identity matrix (the old
-- test) is symmetric and catches nothing: row i = column i.
cyclic3 :: TransitionMatrix 3
cyclic3 =
  either (error . show) id $
    mkTransitionMatrix (S.matrix [0, 1, 0, 0, 0, 1, 1, 0, 0])

spec :: Spec
spec = do
  ---------------------------------------------------------------------------
  -- mkTransitionMatrix: delegates to the distribution predicate, row by row
  ---------------------------------------------------------------------------
  prop "round-trips exactly" $
    forAll (genStochasticSq @3) $ \m ->
      case mkTransitionMatrix m of
        Right p -> property (S.extract (unTransitionMatrix p) == S.extract m)
        Left e  -> counterexample ("generator produced a rejected matrix: " <> show e) False

  prop "rejects a row-sum bumped by 1e-6, naming the row" $
    forAll (genStochasticSq @3) $ \m ->
      case mkTransitionMatrix (onLists (bumpSmallestInRow0 1e-6) m) of
        Left (InRow 0 (SumOffBy _)) -> property True
        Left e  -> counterexample ("expected InRow 0 (SumOffBy _), got " <> show e) False
        Right _ -> counterexample "expected rejection, got Right" False

  prop "rejects a negative entry, naming row and column" $
    forAll (genStochasticSq @3) $ \m ->
      case mkTransitionMatrix (onLists (setEntry00 (-1e-6)) m) of
        Left (InRow 0 (NegativeEntry 0 _)) -> property True
        Left e  -> counterexample ("expected InRow 0 (NegativeEntry 0 _), got " <> show e) False
        Right _ -> counterexample "expected rejection, got Right" False

  ---------------------------------------------------------------------------
  -- mulTransitionMatrix: the closure theorem
  ---------------------------------------------------------------------------
  -- Re-validation is APPROPRIATE here. The implementation must not re-validate
  -- (the theorem is exact); the test must. That is "verified" paired to "proved".
  prop "prop_productRowStochastic: closure under multiplication" $
    forAll ((,) <$> genStochasticSq @3 <*> genStochasticSq @3) $ \(a, b) ->
      case (mkTransitionMatrix a, mkTransitionMatrix b) of
        (Right pa, Right pb) ->
          case mkTransitionMatrix (unTransitionMatrix (mulTransitionMatrix pa pb)) of
            Right _ -> property True
            Left e  -> counterexample ("product rejected: " <> show e) False
        _ -> counterexample "generator produced a rejected matrix" False

  ---------------------------------------------------------------------------
  -- rowAt: the kernel view. Merged in from the old KernelSpec.
  ---------------------------------------------------------------------------
  it "rowAt reads ROWS, not columns" $
    LA.toList (S.extract (unDistribution (rowAt cyclic3 0))) `shouldBe` [0, 1, 0]

  it "rowAt agrees with the cycle on every state" $ do
    let out i = LA.toList (S.extract (unDistribution (rowAt cyclic3 i)))
    out 1 `shouldBe` [0, 0, 1]
    out 2 `shouldBe` [1, 0, 0]

  -- rowAt is total BY THEOREM. This test re-checks the theorem's conclusion:
  -- a row of a validated matrix is always a valid distribution.
  prop "rowAt always yields a valid Distribution" $
    forAll (genStochasticSq @3) $ \m ->
      case mkTransitionMatrix m of
        Right p ->
          conjoin
            [ case mkDistribution (unDistribution (rowAt p i)) of
                Right _ -> property True
                Left e  -> counterexample ("row " <> show i <> " rejected: " <> show e) False
            | i <- [0 .. 2]
            ]
        Left e -> counterexample ("generator produced a rejected matrix: " <> show e) False