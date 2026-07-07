{-# LANGUAGE InstanceSigs #-}
module Dtmc.StochasticSpec where

import Data.Either ( isRight )
import Dtmc.StochErr ( StochErr (..) )
import Dtmc.Stochastic
  ( Stochastic
  , mkStochastic
  , mulStochastic
  , unStochastic
  )
import Numeric.LinearAlgebra ( fromLists )
import Test.Hspec
  ( Spec
  , describe
  , expectationFailure
  , it
  , shouldBe
  , shouldSatisfy
  )
import Test.QuickCheck
  ( Arbitrary (arbitrary)
  , Gen
  , choose
  , chooseInt
  , frequency
  , property
  , vectorOf
  )

spec :: Spec
spec = do
  describe "mkStochastic" $ do
    it "accepts a valid row-stochastic matrix" $ do
      let matrix =
            fromLists
              [ [0.5, 0.5]
              , [0.2, 0.8]
              ]

      mkStochastic matrix `shouldSatisfy` isRight

    it "accepts a valid row-stochastic matrix with zero entries" $ do
      let matrix =
            fromLists
              [ [1.0, 0.0]
              , [0.0, 1.0]
              ]

      mkStochastic matrix `shouldSatisfy` isRight

    it "rejects a matrix with a negative entry" $ do
      let matrix =
            fromLists
              [ [1.0, 0.0]
              , [-0.1, 1.1]
              ]

      mkStochastic matrix
        `shouldFailWith` NegativeEntry
          { row = 1
          , col = 0
          , val = -0.1
          }

    it "rejects a matrix with an entry above 1" $ do
      let matrix =
            fromLists
              [ [1.1, 0.0]
              , [0.0, 1.0]
              ]

      mkStochastic matrix
        `shouldFailWith` EntryAboveOne
          { row = 0
          , col = 0
          , val = 1.1
          }

    it "rejects a matrix whose rows do not sum to 1" $ do
      let matrix =
            fromLists
              [ [0.5, 0.4]
              , [0.2, 0.8]
              ]

      mkStochastic matrix
        `shouldFailWith` RowSumOffBy
          { row = 0
          , rowSum = 0.9
          }

    it "rejects a non-square matrix" $ do
      let matrix =
            fromLists
              [ [0.5, 0.5, 0.0]
              , [0.2, 0.3, 0.5]
              ]

      mkStochastic matrix
        `shouldFailWith` NonSquareMatrix
          { rowCount = 2
          , colCount = 3
          }

  describe "mulStochastic" $ do
    it "the product of two stochastic matrices is stochastic" $
      property prop_productRowStochastic

prop_productRowStochastic :: SameSizedStochastic -> Bool
prop_productRowStochastic (SameSizedStochastic a b) =
  isRight (mkStochastic (unStochastic (mulStochastic a b)))

data SameSizedStochastic =
  SameSizedStochastic Stochastic Stochastic
  deriving (Show)

instance Arbitrary SameSizedStochastic where
  arbitrary :: Gen SameSizedStochastic
  arbitrary = do
    n <- chooseInt (1, 6)
    a <- genStochastic n
    b <- genStochastic n
    pure (SameSizedStochastic a b)

genStochastic :: Int -> Gen Stochastic
genStochastic n = do
  rawRows <- vectorOf n (genNonZeroRow n)
  let matrix = fromLists (map normalise rawRows)

  case mkStochastic matrix of
    Right stochastic -> pure stochastic
    Left err -> error ("genStochastic produced a non-stochastic matrix: " <> show err)

genNonZeroRow :: Int -> Gen [Double]
genNonZeroRow n = do
  row <- vectorOf n genEntry
  if sum row == 0.0
    then genNonZeroRow n
    else pure row

genEntry :: Gen Double
genEntry =
  frequency
    [ (3, pure 0.0)
    , (7, choose (0.0, 1000.0))
    ]

normalise :: [Double] -> [Double]
normalise row =
  let rowTotal = sum row
   in map (/ rowTotal) row

shouldFailWith :: Either StochErr a -> StochErr -> IO ()
shouldFailWith result expectedErr =
  case result of
    Left actualErr ->
      actualErr `shouldBe` expectedErr
    Right _ ->
      expectationFailure "Expected validation to fail, but it succeeded."