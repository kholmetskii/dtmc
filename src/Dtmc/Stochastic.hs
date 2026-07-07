module Dtmc.Stochastic
    ( Stochastic
    , unStochastic
    , mkStochastic
    , mulStochastic
    ) where

import Numeric.LinearAlgebra (Matrix, cols, rows, toLists)
import qualified Numeric.LinearAlgebra as LA

-- | A numerically row-stochastic matrix.
--
-- A matrix is accepted as row-stochastic when:
--
-- * it is square;
-- * every entry is at least @-epsilon@;
-- * every entry is at most @1 + epsilon@;
-- * every row sums to @1@ up to @epsilon@.
--
-- Therefore, a value of type 'Stochastic' does not mean that the rows sum to
-- exactly @1@. It means the matrix passed numerical stochastic validation at a
-- constructor boundary.
--
-- This representation uses 'Double' because the spectral layer of the library
-- will compute eigenvectors, stationary distributions, and
-- Perron--Frobenius-style quantities, which are generally non-rational and
-- therefore numerical.
--
-- For a discrete-time Markov chain, a stochastic matrix is a valid transition
-- matrix.
newtype Stochastic = Stochastic
  { 
    unStochastic :: Matrix Double
  }
  deriving (Eq, Show)

-- | Numerical tolerance used when checking stochastic matrices.
--
-- The tolerance should be larger than ordinary floating-point normalisation
-- error, but much smaller than meaningful modelling error.
epsilon :: Double
epsilon = 1e-9

-- | Smart constructor for stochastic matrices.
--
-- This is the only safe way to construct a 'Stochastic' value from a raw
-- matrix. It checks that the matrix is square and numerically row-stochastic.
mkStochastic :: Matrix Double -> Maybe Stochastic
mkStochastic matrix
  | not (isSquare matrix) = Nothing
  | not (allEntriesInProbabilityRange matrix) = Nothing
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
-- Numerically, floating-point multiplication may introduce tiny drift. This
-- library represents stochasticity up to 'epsilon', so the invariant is
-- understood numerically rather than exactly.
--
-- Warning: this function assumes the two matrices have compatible dimensions.
mulStochastic :: Stochastic -> Stochastic -> Stochastic
mulStochastic a b =
  Stochastic (unStochastic a LA.<> unStochastic b)

isSquare :: Matrix Double -> Bool
isSquare matrix =
  rows matrix == cols matrix

allEntriesInProbabilityRange :: Matrix Double -> Bool
allEntriesInProbabilityRange matrix =
  all (all approximatelyProbability) (toLists matrix)

approximatelyProbability :: Double -> Bool
approximatelyProbability x =
  x >= -epsilon && x <= 1.0 + epsilon

allRowsSumToOne :: Matrix Double -> Bool
allRowsSumToOne matrix =
  all rowSumsToOne (toLists matrix)

rowSumsToOne :: [Double] -> Bool
rowSumsToOne row =
  abs (sum row - 1.0) <= epsilon