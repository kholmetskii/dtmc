module Dtmc.Stochastic
    ( Stochastic
    , unStochastic
    , mkStochastic
    , mulStochastic
    ) where

import Numeric.LinearAlgebra ( Matrix, cols, rows, toLists )
import qualified Numeric.LinearAlgebra as LA

-- | A row-stochastic matrix.
--
-- A matrix is row-stochastic when:
--
-- * it is square;
-- * every entry is non-negative;
-- * every row sums to 1.
--
-- For a discrete-time Markov chain, a stochastic matrix is a valid transition
-- matrix.
newtype Stochastic = Stochastic
  { 
    unStochastic :: Matrix Double
  }
  deriving (Eq, Show)

-- | Numerical tolerance used when checking stochastic matrices.
epsilon :: Double
epsilon = 1e-9

-- | Smart constructor for stochastic matrices.
--
-- This is the only safe way to construct a 'Stochastic' value from a raw
-- matrix. It checks that the matrix is square, non-negative, and row-stochastic.
mkStochastic :: Matrix Double -> Maybe Stochastic
mkStochastic matrix
  | not (isSquare matrix) = Nothing
  | not (allEntriesNonNegative matrix) = Nothing
  | not (allRowsSumToOne matrix) = Nothing
  | otherwise = Just (Stochastic matrix)

-- | Multiply two stochastic matrices.
--
-- Proof:
--
-- Let A and B be row-stochastic matrices of the same size. Since every entry of
-- A and B is non-negative, every entry of AB is a sum of products of
-- non-negative numbers, so every entry of AB is non-negative.
--
-- Also, for each row i,
--
-- @
-- sum_j (AB)_ij
--   = sum_j sum_k A_ik B_kj
--   = sum_k A_ik sum_j B_kj
--   = sum_k A_ik
--   = 1.
-- @
--
-- Therefore AB is row-stochastic.
--
-- Warning: this function assumes the two matrices have compatible dimensions.
mulStochastic :: Stochastic -> Stochastic -> Stochastic
mulStochastic a b =
  Stochastic (unStochastic a LA.<> unStochastic b)

isSquare :: Matrix Double -> Bool
isSquare matrix =
  rows matrix == cols matrix

allEntriesNonNegative :: Matrix Double -> Bool
allEntriesNonNegative matrix =
  all (all approximatelyNonNegative) (toLists matrix)

approximatelyNonNegative :: Double -> Bool
approximatelyNonNegative x =
  x >= -epsilon

allRowsSumToOne :: Matrix Double -> Bool
allRowsSumToOne matrix =
  all rowSumsToOne (toLists matrix)

rowSumsToOne :: [Double] -> Bool
rowSumsToOne row =
  abs (sum row - 1.0) <= epsilon