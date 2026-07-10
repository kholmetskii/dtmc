module Dtmc.KernelSpec (spec) where

import Dtmc.Distribution (mkDistribution, unDistribution)
import Dtmc.Generators (genStochasticSq)
import Dtmc.Kernel (rowAt)
import Dtmc.StochasticMatrix (StochasticMatrix, mkStochasticMatrix)
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Static as S
import Test.Hspec
import Test.Hspec.QuickCheck (prop)
import Test.QuickCheck

-- | The 3-cycle: row 0 = e₁, row 1 = e₂, row 2 = e₀.
--
--   under P  : 0 → 1 → 2 → 0
--   under Pᵀ : 0 → 2 → 1 → 0
--
-- Asymmetric, so it catches confused rows/columns. The identity matrix (the old
-- test) is symmetric and catches nothing: row i = column i. docs/TESTING.md T2.
cyclic3 :: StochasticMatrix 3
cyclic3 =
  either (error . show) id $
    mkStochasticMatrix (S.matrix [0, 1, 0, 0, 0, 1, 1, 0, 0])

spec :: Spec
spec = do
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
      case mkStochasticMatrix m of
        Right sm ->
          conjoin
            [ case mkDistribution (unDistribution (rowAt sm i)) of
                Right _ -> property True
                Left e  -> counterexample ("row " <> show i <> " rejected: " <> show e) False
            | i <- [0 .. 2]
            ]
        Left e -> counterexample ("generator produced a rejected matrix: " <> show e) False
