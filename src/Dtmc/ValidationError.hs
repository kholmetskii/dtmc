module Dtmc.ValidationError
  ( ValidationError (..)
  ) where

data ValidationError
  = VectorDimensionMismatch
      { expectedLength :: Int
      , actualLength :: Int
      }
  | NonSquareMatrix
      { rowCount :: Int
      , colCount :: Int
      }
  | MatrixDimensionMismatch
      { expectedSize :: Int
      , actualRows :: Int
      , actualCols :: Int
      }
  | NegativeEntry
      { row :: Int
      , col :: Int
      , val :: Double
      }
  | EntryAboveOne
      { row :: Int
      , col :: Int
      , val :: Double
      }
  | RowSumOffBy
      { row :: Int
      , rowSum :: Double
      }
  deriving (Eq, Show)