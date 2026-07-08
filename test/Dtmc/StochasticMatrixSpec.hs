{-# LANGUAGE DataKinds #-}
{-# LANGUAGE InstanceSigs #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.StochasticMatrixSpec where

import Data.Either (isRight)
import Dtmc.StochasticMatrix
  ( StochasticMatrix
  , mkStochasticMatrix
  , mulStochasticMatrix
  , unStochasticMatrix
  )
import Dtmc.ValidationError ( ValidationError( .. ) )
import Numeric.LinearAlgebra (fromLists)
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
import Dtmc.Generators (genDenseStochasticMatrix2)

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

  describe "mulStochasticMatrix" $ do
    it "the product of two stochastic matrices is stochastic" $
      property prop_productRowStochastic

prop_productRowStochastic :: SameSizedStochastic -> Bool
prop_productRowStochastic (SameSizedStochastic a b) =
  isRight
    ( mkStochasticMatrix @2
        (unStochasticMatrix (mulStochasticMatrix a b))
    )

data SameSizedStochastic =
  SameSizedStochastic (StochasticMatrix 2) (StochasticMatrix 2)
  deriving (Show)

instance Arbitrary SameSizedStochastic where
  arbitrary :: Gen SameSizedStochastic
  arbitrary = do
    a <- genDenseStochasticMatrix2
    b <- genDenseStochasticMatrix2
    pure (SameSizedStochastic a b)

shouldFailWith :: Either ValidationError a -> ValidationError -> IO ()
shouldFailWith result expectedErr =
  case result of
    Left actualErr ->
      actualErr `shouldBe` expectedErr
    Right _ ->
      expectationFailure "Expected validation to fail, but it succeeded."