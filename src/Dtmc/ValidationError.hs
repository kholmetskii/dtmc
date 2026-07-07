module Dtmc.ValidationError
  ( ValidationError (..)
  ) where

data ValidationError
  = NonSquareMatrix
      { rowCount :: Int
      , colCount :: Int
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