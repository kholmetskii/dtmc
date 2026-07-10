module Dtmc.DistributionSpec (spec) where

import Dtmc.Distribution
  ( DistributionError (..)
  , SimplexError (..)
  , approxDistributionEq
  , mkDistribution
  , unDistribution
  , validateSimplex
  )
import Dtmc.TestSupport
  ( bumpSmallest
  , genSimplexPoint
  )
import qualified Numeric.LinearAlgebra.Static as S
import Test.Hspec
import Test.Hspec.QuickCheck (prop)
import Test.QuickCheck

spec :: Spec
spec = do
  describe "validateSimplex" $ do
    it "rejects an empty vector" $
      validateSimplex (S.vector [] :: S.R 0)
        `shouldBe` Left (SumOffBy 0)

    it "accepts a tiny negative rounding error" $
      validateSimplex (S.vector [-1e-17, 1] :: S.R 2)
        `shouldBe` Right ()

    it "reports an entry above one" $
      validateSimplex (S.vector [1.5, -0.5] :: S.R 2)
        `shouldBe` Left (EntryAboveOne 0 1.5)

    prop "accepts normalised vectors" $
      forAll (genSimplexPoint 3) $ \entries ->
        validateSimplex (S.vector entries :: S.R 3)
          === Right ()

    prop "rejects vectors whose sum is too large" $
      forAll (genSimplexPoint 3) $ \entries ->
        case validateSimplex
          (S.vector (bumpSmallest 1e-6 entries) :: S.R 3) of
          Left (SumOffBy _) ->
            property True
          result ->
            counterexample
              ("expected SumOffBy, got " <> show result)
              False

    prop "rejects genuinely negative entries" $
      forAll (genSimplexPoint 3) $ \entries ->
        let invalid = case entries of
              _ : rest -> (-1e-6) : rest
              [] -> []
         in case validateSimplex (S.vector invalid :: S.R 3) of
              Left (NegativeEntry 0 _) ->
                property True
              result ->
                counterexample
                  ("expected NegativeEntry 0, got " <> show result)
                  False

  describe "mkDistribution" $ do
    prop "preserves the validated vector" $
      forAll (genSimplexPoint 3) $ \entries ->
        let simplexVector = S.vector entries :: S.R 3
        in case mkDistribution simplexVector of
              Right distribution ->
                S.extract (unDistribution distribution)
                  === S.extract simplexVector
              Left err ->
                counterexample
                  ("generated vector was rejected: " <> show err)
                  False

    prop "wraps simplex validation errors" $
      forAll (genSimplexPoint 3) $ \entries ->
        case mkDistribution
          (S.vector (bumpSmallest 1e-6 entries) :: S.R 3) of
          Left (DistributionError (SumOffBy _)) ->
            property True
          result ->
            counterexample
              ("expected DistributionError SumOffBy, got " <> show result)
              False

    it "rejects a zero-dimensional distribution" $
      case mkDistribution (S.vector [] :: S.R 0) of
        Left (DistributionError (SumOffBy actualSum)) ->
          actualSum `shouldBe` 0
        Left err ->
          expectationFailure
            ("expected SumOffBy 0, got " <> show err)
        Right _ ->
          expectationFailure "expected rejection"

  describe "approxDistributionEq" $ do
    prop "is reflexive at zero tolerance" $
      forAll (genSimplexPoint 3) $ \entries ->
        case mkDistribution (S.vector entries :: S.R 3) of
          Right distribution ->
            property
              (approxDistributionEq 0 distribution distribution)
          Left err ->
            counterexample
              ("generated vector was rejected: " <> show err)
              False
