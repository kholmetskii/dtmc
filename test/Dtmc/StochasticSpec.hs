module Dtmc.StochasticSpec where

import Dtmc.Stochastic
import Numeric.LinearAlgebra (fromLists)
import Test.Hspec
import Test.QuickCheck

spec :: Spec
spec = do
  describe "mkStochastic" $ do
    it "accepts a valid row-stochastic matrix" $ do
      let matrix =
            fromLists
              [ [0.5, 0.5]
              , [0.2, 0.8]
              ]

      mkStochastic matrix `shouldSatisfy` isJust

    it "rejects a matrix with a negative entry" $ do
      let matrix =
            fromLists
              [ [1.0, 0.0]
              , [-0.1, 1.1]
              ]

      mkStochastic matrix `shouldBe` Nothing

    it "rejects a matrix whose rows do not sum to 1" $ do
      let matrix =
            fromLists
              [ [0.5, 0.4]
              , [0.2, 0.8]
              ]

      mkStochastic matrix `shouldBe` Nothing

    it "rejects a non-square matrix" $ do
      let matrix =
            fromLists
              [ [0.5, 0.5, 0.0]
              , [0.2, 0.3, 0.5]
              ]

      mkStochastic matrix `shouldBe` Nothing

  describe "mulStochastic" $ do
    it "the product of two stochastic matrices is stochastic" $
      property prop_productRowStochastic

prop_productRowStochastic :: SameSizedStochastic -> Bool
prop_productRowStochastic (SameSizedStochastic a b) =
  isJust (mkStochastic (unStochastic (mulStochastic a b)))

data SameSizedStochastic =
  SameSizedStochastic Stochastic Stochastic
  deriving (Show)

instance Arbitrary SameSizedStochastic where
  arbitrary = do
    n <- chooseInt (1, 6)
    a <- genStochastic n
    b <- genStochastic n
    pure (SameSizedStochastic a b)

genStochastic :: Int -> Gen Stochastic
genStochastic n = do
  rawRows <- vectorOf n (vectorOf n (choose (0.0, 1000.0)))
  let matrix = fromLists (map normalise rawRows)

  case mkStochastic matrix of
    Just stochastic -> pure stochastic
    Nothing -> error "genStochastic produced a non-stochastic matrix"

normalise :: [Double] -> [Double]
normalise row =
  let positiveRow = map (+ 1.0) row
      rowTotal = sum positiveRow
   in map (/ rowTotal) positiveRow

isJust :: Maybe a -> Bool
isJust (Just _) = True
isJust Nothing = False
