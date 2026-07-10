module Dtmc.SimplexSpec (spec) where

import Dtmc.Generators (bumpSmallest, genSimplexPointList)
import Dtmc.Simplex
import qualified Numeric.LinearAlgebra.Static as S
import Test.Hspec
import Test.Hspec.QuickCheck (prop)
import Test.QuickCheck

-- Verifies: Dtmc.Simplex.validateSimplexPoint
--
-- hspec-discover already wraps spec in describe "Dtmc.Simplex".
spec :: Spec
spec = do
  it "rejects the empty point (n = 0): sum [] = 0 /= 1" $
    validateSimplexPoint (S.vector [] :: S.R 0) `shouldBe` Left (SumOffBy 0.0)

  it "tolerates a -1e-17 rounding sliver" $
    validateSimplexPoint (S.vector [-1e-17, 1.0] :: S.R 2) `shouldBe` Right ()

  it "reports EntryAboveOne before NegativeEntry (check order is contractual)" $
    validateSimplexPoint (S.vector [1.5, -0.5] :: S.R 2)
      `shouldBe` Left (EntryAboveOne 0 1.5)

  -- Regression on the QuickCheck counterexample (seed 180161876). docs/TESTING.md T1.
  it "an entry at 1.0 bumped past the bound reports EntryAboveOne, not SumOffBy" $
    validateSimplexPoint (S.vector [1.000001, 0.0, 0.0] :: S.R 3)
      `shouldBe` Left (EntryAboveOne 0 1.000001)

  prop "accepts a normalised point" $
    forAll (genSimplexPointList 3) $ \xs ->
      validateSimplexPoint (S.vector xs :: S.R 3) === Right ()

  prop "rejects a sum bumped by 1e-6 (>> tolerance) with SumOffBy" $
    forAll (genSimplexPointList 3) $ \xs ->
      case validateSimplexPoint (S.vector (bumpSmallest 1e-6 xs) :: S.R 3) of
        Left (SumOffBy _) -> property True
        other -> counterexample ("expected SumOffBy, got " <> show other) False

  prop "rejects a genuinely negative coordinate with NegativeEntry" $
    forAll (genSimplexPointList 3) $ \xs ->
      let broken = case xs of
            (_ : rest) -> (-1e-6) : rest
            []         -> []
       in case validateSimplexPoint (S.vector broken :: S.R 3) of
            Left (NegativeEntry 0 _) -> property True
            other -> counterexample ("expected NegativeEntry 0, got " <> show other) False
