{-# LANGUAGE DataKinds #-}
{-# LANGUAGE InstanceSigs #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.StochasticMatrixSpec where

import Data.Either
  ( isRight
  )
import Dtmc.Generators
  ( genDenseStochasticMatrix2
  , genDenseStochasticRawMatrix2
  )
import Dtmc.StochasticMatrix
  ( StochasticMatrix
  , mkStochasticMatrix
  , mulStochasticMatrix
  , unStochasticMatrix
  )
import Dtmc.ValidationError
  ( ValidationError ( .. )
  )
import Numeric.LinearAlgebra
  ( Matrix
  , fromLists
  , toLists
  )
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
  , property
  )

spec :: Spec
spec = do
  describe "mkStochasticMatrix" $ do
    it "accepts a valid row-stochastic matrix" $ do
      let matrix =
            fromLists
              [ [0.5, 0.5]
              , [0.2, 0.8]
              ]

      mkStochasticMatrix @2 matrix `shouldSatisfy` isRight

    it "accepts a valid row-stochastic matrix with zero entries" $ do
      let matrix =
            fromLists
              [ [1.0, 0.0]
              , [0.0, 1.0]
              ]

      mkStochasticMatrix @2 matrix `shouldSatisfy` isRight

    it "rejects a matrix with a negative entry" $ do
      let matrix =
            fromLists
              [ [1.0, 0.0]
              , [-0.1, 1.1]
              ]

      mkStochasticMatrix @2 matrix
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

      mkStochasticMatrix @2 matrix
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

      mkStochasticMatrix @2 matrix
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

      mkStochasticMatrix @2 matrix
        `shouldFailWith` NonSquareMatrix
          { rowCount = 2
          , colCount = 3
          }

    it "rejects a square matrix with the wrong type-level dimension" $ do
      let matrix =
            fromLists
              [ [0.5, 0.5]
              , [0.2, 0.8]
              ]

      mkStochasticMatrix @3 matrix
        `shouldFailWith` MatrixDimensionMismatch
          { expectedSize = 3
          , actualRows = 2
          , actualCols = 2
          }

    it "round-trips generated stochastic matrices without mutating them" $
      property prop_generatedMatrixRoundTrips

    it "rejects row-sum perturbations with RowSumOffBy" $
      property prop_rejectsRowSumBump

    it "rejects sign flips with NegativeEntry" $
      property prop_rejectsNegativeEntry

  describe "mulStochasticMatrix" $ do
    it "the product of two stochastic matrices is stochastic" $
      property prop_productRowStochastic

prop_productRowStochastic :: SameSizedStochastic -> Bool
prop_productRowStochastic (SameSizedStochastic a b) =
  isRight
    ( mkStochasticMatrix @2
        (unStochasticMatrix (mulStochasticMatrix a b))
    )

prop_generatedMatrixRoundTrips :: DenseRawMatrix2 -> Bool
prop_generatedMatrixRoundTrips (DenseRawMatrix2 matrix) =
  case mkStochasticMatrix @2 matrix of
    Right stochasticMatrix ->
      unStochasticMatrix stochasticMatrix == matrix
    Left _ ->
      False

prop_rejectsRowSumBump :: DenseRawMatrix2 -> Bool
prop_rejectsRowSumBump (DenseRawMatrix2 matrix) =
  case mkStochasticMatrix @2 (bumpSafeEntryInRow0 delta matrix) of
    Left RowSumOffBy { row = 0 } ->
      True
    _ ->
      False

prop_rejectsNegativeEntry :: DenseRawMatrix2 -> Bool
prop_rejectsNegativeEntry (DenseRawMatrix2 matrix) =
  case mkStochasticMatrix @2 (setEntry00 (-delta) matrix) of
    Left NegativeEntry { row = 0, col = 0 } ->
      True
    _ ->
      False

delta :: Double
delta = 1e-6

data SameSizedStochastic =
  SameSizedStochastic (StochasticMatrix 2) (StochasticMatrix 2)
  deriving (Show)

instance Arbitrary SameSizedStochastic where
  arbitrary :: Gen SameSizedStochastic
  arbitrary = do
    a <- genDenseStochasticMatrix2
    b <- genDenseStochasticMatrix2
    pure (SameSizedStochastic a b)

newtype DenseRawMatrix2 =
  DenseRawMatrix2 (Matrix Double)
  deriving (Show)

instance Arbitrary DenseRawMatrix2 where
  arbitrary :: Gen DenseRawMatrix2
  arbitrary =
    DenseRawMatrix2 <$> genDenseStochasticRawMatrix2

bumpSafeEntryInRow0 :: Double -> Matrix Double -> Matrix Double
bumpSafeEntryInRow0 amount matrix =
  case toLists matrix of
    ([x, y] : rowsRest)
      | x <= y ->
          fromLists ([x + amount, y] : rowsRest)
      | otherwise ->
          fromLists ([x, y + amount] : rowsRest)
    _ ->
      matrix

setEntry00 :: Double -> Matrix Double -> Matrix Double
setEntry00 newValue matrix =
  case toLists matrix of
    (row0 : rowsRest) ->
      case row0 of
        (_ : xs) ->
          fromLists ((newValue : xs) : rowsRest)
        [] ->
          matrix
    [] ->
      matrix

shouldFailWith :: Either ValidationError a -> ValidationError -> IO ()
shouldFailWith result expectedErr =
  case result of
    Left actualErr ->
      actualErr `shouldBe` expectedErr
    Right _ ->
      expectationFailure "Expected validation to fail, but it succeeded."