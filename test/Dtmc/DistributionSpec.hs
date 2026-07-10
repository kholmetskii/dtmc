module Dtmc.DistributionSpec
  ( spec
  ) where

import Dtmc
  ( DistributionError (..)
  , SimplexError (..)
  , approxDistributionEq
  , mkDistribution
  , unDistribution
  )
import Dtmc.TestSupport
  ( bumpSmallest
  , genSimplexPoint
  )
import qualified Numeric.LinearAlgebra.Static as S
import Test.Hspec
  ( Spec
  , describe
  , expectationFailure
  , it
  , shouldBe
  )
import Test.Hspec.QuickCheck
  ( prop
  )
import Test.QuickCheck
  ( counterexample
  , forAll
  , property
  , (===)
  )

spec :: Spec
spec = do
  describe "mkDistribution" $ do
    it "rejects an empty vector" $
      case mkDistribution (S.vector [] :: S.R 0) of
        Left err ->
          err `shouldBe` DistributionError (SumOffBy 0)
        Right _ ->
          expectationFailure "expected rejection"

    it "accepts a tiny negative rounding error" $
      case mkDistribution (S.vector [-1e-17, 1] :: S.R 2) of
        Right _ ->
          pure ()
        Left err ->
          expectationFailure
            ("expected acceptance, got " <> show err)

    it "reports an entry above one" $
      case mkDistribution (S.vector [1.5, -0.5] :: S.R 2) of
        Left err ->
          err `shouldBe`
            DistributionError (EntryAboveOne 0 1.5)
        Right _ ->
          expectationFailure "expected rejection"

    prop "accepts normalised vectors" $
      forAll (genSimplexPoint 3) $ \entries ->
        case mkDistribution (S.vector entries :: S.R 3) of
          Right _ ->
            property True
          Left err ->
            counterexample
              ("generated vector was rejected: " <> show err)
              False

    prop "rejects vectors whose sum is too large" $
      forAll (genSimplexPoint 3) $ \entries ->
        case mkDistribution
          (S.vector (bumpSmallest 1e-6 entries) :: S.R 3) of
          Left (DistributionError (SumOffBy _)) ->
            property True
          result ->
            counterexample
              ("expected DistributionError SumOffBy, got " <> show result)
              False

    prop "rejects genuinely negative entries" $
      forAll (genSimplexPoint 3) $ \entries ->
        let invalid =
              case entries of
                _ : rest -> (-1e-6) : rest
                [] -> []
         in case mkDistribution (S.vector invalid :: S.R 3) of
              Left (DistributionError (NegativeEntry 0 _)) ->
                property True
              result ->
                counterexample
                  ("expected DistributionError NegativeEntry 0, got " <> show result)
                  False

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

  describe "approxDistributionEq" $
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