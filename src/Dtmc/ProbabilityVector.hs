module Dtmc.ProbabilityVector
  ( ProbVector
  , unProbVector
  , mkProbVector
  , mkProbVectorAt
  , epsilon
  ) where

import Dtmc.StochErr
  ( StochErr
      ( EntryAboveOne
      , NegativeEntry
      , RowSumOffBy
      )
  )
import Numeric.LinearAlgebra (Vector, toList)

-- | A numerically valid probability vector.
--
-- A vector is accepted as a probability vector when:
--
-- * every entry is at least @-epsilon@;
-- * every entry is at most @1 + epsilon@;
-- * the entries sum to @1@ up to @epsilon@.
--
-- Therefore, a value of type 'ProbVector' does not mean that the entries sum
-- to exactly @1@. It means the vector passed numerical probability validation
-- at a constructor boundary.
newtype ProbVector = ProbVector
  { unProbVector :: Vector Double
  }
  deriving (Show)

-- | Numerical tolerance used when checking probability vectors.
--
-- The tolerance should be larger than ordinary floating-point normalisation
-- error, but much smaller than meaningful modelling error.
epsilon :: Double
epsilon = 1e-9

-- | Smart constructor for probability vectors.
--
-- This version uses row index @0@ in error messages. Matrix validation should
-- usually use 'mkProbVectorAt' so the offending matrix row can be reported.
mkProbVector :: Vector Double -> Either StochErr ProbVector
mkProbVector =
  mkProbVectorAt 0

-- | Smart constructor for probability vectors with an explicit row index.
--
-- The row index is used only for informative error reporting when a probability
-- vector is being checked as a row of a stochastic matrix.
mkProbVectorAt :: Int -> Vector Double -> Either StochErr ProbVector
mkProbVectorAt rowIndex vector =
  case firstBadEntry rowIndex (toList vector) of
    Just err -> Left err
    Nothing ->
      let total = sum (toList vector)
       in if abs (total - 1.0) <= epsilon
            then Right (ProbVector vector)
            else Left (RowSumOffBy rowIndex total)

firstBadEntry :: Int -> [Double] -> Maybe StochErr
firstBadEntry rowIndex entries =
  firstBadEntryFrom rowIndex 0 entries

firstBadEntryFrom :: Int -> Int -> [Double] -> Maybe StochErr
firstBadEntryFrom _ _ [] =
  Nothing
firstBadEntryFrom rowIndex colIndex (x : xs)
  | x < -epsilon =
      Just (NegativeEntry rowIndex colIndex x)
  | x > 1.0 + epsilon =
      Just (EntryAboveOne rowIndex colIndex x)
  | otherwise =
      firstBadEntryFrom rowIndex (colIndex + 1) xs