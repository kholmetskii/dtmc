module Dtmc.DistributionSpec (spec) where

import Dtmc.Distribution
import Dtmc.Generators (bumpSmallest, genSimplexPointList)
import Dtmc.Simplex (SimplexError (..))
import qualified Numeric.LinearAlgebra.Static as S
import Test.Hspec
import Test.Hspec.QuickCheck (prop)
import Test.QuickCheck

-- Verifies: Dtmc.Distribution.mkDistribution, approxDistributionEq
spec :: Spec
spec = do
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

  prop "wraps a negative coordinate as NegativeEntry" $
    forAll (genSimplexPointList 3) $ \xs ->
      let broken = case xs of
            (_ : rest) -> (-1e-6) : rest
            []         -> []
       in case mkDistribution (S.vector broken :: S.R 3) of
            Left (DistributionError (NegativeEntry 0 _)) -> property True
            Left e  -> counterexample ("expected NegativeEntry 0, got " <> show e) False
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
