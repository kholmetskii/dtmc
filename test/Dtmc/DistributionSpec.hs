module Dtmc.DistributionSpec (spec) where

import Dtmc.Distribution
    ( DistributionError(DistributionError),
      SimplexError(SumOffBy, EntryAboveOne, NegativeEntry),
      Distribution(unDistribution),
      validateSimplex,
      mkDistribution,
      approxDistributionEq )
import Dtmc.Generators (bumpSmallest, genSimplexPointList)
import qualified Numeric.LinearAlgebra.Static as S
import Test.Hspec
import Test.Hspec.QuickCheck (prop)
import Test.QuickCheck

-- Verifies: Dtmc.Distribution.validateSimplex, mkDistribution, approxDistributionEq
--
-- The old SimplexSpec merged in here: the predicate is no longer a separate
-- module. hspec-discover already wraps spec in describe "Dtmc.Distribution".
spec :: Spec
spec = do
  ---------------------------------------------------------------------------
  -- validateSimplex: the Δ^{n-1} predicate
  ---------------------------------------------------------------------------
  it "rejects the empty point (n = 0): sum [] = 0 /= 1" $
    validateSimplex (S.vector [] :: S.R 0) `shouldBe` Left (SumOffBy 0.0)

  it "tolerates a -1e-17 rounding sliver" $
    validateSimplex (S.vector [-1e-17, 1.0] :: S.R 2) `shouldBe` Right ()

  it "reports EntryAboveOne before NegativeEntry (check order is contractual)" $
    validateSimplex (S.vector [1.5, -0.5] :: S.R 2)
      `shouldBe` Left (EntryAboveOne 0 1.5)

  -- Regression on the QuickCheck counterexample (seed 180161876). docs/TESTING.md T1.
  it "an entry at 1.0 bumped past the bound reports EntryAboveOne, not SumOffBy" $
    validateSimplex (S.vector [1.000001, 0.0, 0.0] :: S.R 3)
      `shouldBe` Left (EntryAboveOne 0 1.000001)

  prop "accepts a normalised point" $
    forAll (genSimplexPointList 3) $ \xs ->
      validateSimplex (S.vector xs :: S.R 3) === Right ()

  prop "rejects a sum bumped by 1e-6 (>> tolerance) with SumOffBy" $
    forAll (genSimplexPointList 3) $ \xs ->
      case validateSimplex (S.vector (bumpSmallest 1e-6 xs) :: S.R 3) of
        Left (SumOffBy _) -> property True
        other -> counterexample ("expected SumOffBy, got " <> show other) False

  prop "rejects a genuinely negative coordinate with NegativeEntry" $
    forAll (genSimplexPointList 3) $ \xs ->
      let broken = case xs of
            (_ : rest) -> (-1e-6) : rest
            []         -> []
       in case validateSimplex (S.vector broken :: S.R 3) of
            Left (NegativeEntry 0 _) -> property True
            other -> counterexample ("expected NegativeEntry 0, got " <> show other) False

  ---------------------------------------------------------------------------
  -- mkDistribution: the validation boundary
  ---------------------------------------------------------------------------
  -- The constructor validates but does NOT mutate, hence exact equality.
  prop "round-trips exactly" $
    forAll (genSimplexPointList 3) $ \xs ->
      let v = S.vector xs :: S.R 3
       in case mkDistribution v of
            Right d -> property (S.extract (unDistribution d) == S.extract v)
            Left e  -> counterexample ("generator produced a rejected point: " <> show e) False

  prop "wraps the simplex error without inventing a row index" $
    forAll (genSimplexPointList 3) $ \xs ->
      case mkDistribution (S.vector (bumpSmallest 1e-6 xs) :: S.R 3) of
        Left (DistributionError (SumOffBy _)) -> property True
        Left e  -> counterexample ("expected SumOffBy, got " <> show e) False
        Right _ -> counterexample "expected rejection, got Right" False

  prop "approxDistributionEq is reflexive at tolerance 0" $
    forAll (genSimplexPointList 3) $ \xs ->
      case mkDistribution (S.vector xs :: S.R 3) of
        Right d -> property (approxDistributionEq 0 d d)
        Left e  -> counterexample ("generator produced a rejected point: " <> show e) False

  it "rejects a zero-dimensional distribution" $
    case mkDistribution (S.vector [] :: S.R 0) of
      Left (DistributionError (SumOffBy s)) -> s `shouldBe` 0.0
      Left e  -> expectationFailure ("expected SumOffBy 0.0, got " <> show e)
      Right _ -> expectationFailure "expected rejection of R 0"