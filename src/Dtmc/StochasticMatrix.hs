{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE KindSignatures #-}


module Dtmc.StochasticMatrix
  ( StochasticMatrix
  , unStochasticMatrix
  , mkStochasticMatrix
  , mulStochasticMatrix
  , approxStochasticMatrixEq
  ) where

import Data.Proxy (Proxy (Proxy))
import Dtmc.ProbabilityVector (mkProbabilityVectorAt)
import Dtmc.ValidationError (ValidationError (..))
import GHC.TypeNats (KnownNat, Nat, natVal)
import Numeric.LinearAlgebra (Matrix, cols, rows, toRows)
import qualified Numeric.LinearAlgebra as LA

-- | A numerically row-stochastic @n x n@ matrix.
--
-- A matrix is accepted as a stochastic matrix when:
--
-- * its runtime shape is @n x n@;
-- * each row is accepted as a probability vector of length @n@.
--
-- Therefore, a value of type @StochasticMatrix n@ does not mean that the rows
-- sum to exactly @1@. It means the matrix passed numerical stochastic
-- validation at a constructor boundary.
--
-- This representation uses 'Double' because the spectral layer of the library
-- will compute eigenvectors, stationary distributions, and
-- Perron--Frobenius-style quantities, which are generally non-rational and
-- therefore numerical.
--
-- For a discrete-time Markov chain, a stochastic matrix is a valid transition
-- matrix.
newtype StochasticMatrix (n :: Nat) = StochasticMatrix
  { unStochasticMatrix :: Matrix Double
  }
  deriving (Show)

-- | Smart constructor for stochastic matrices.
--
-- This is the only safe way to construct a @StochasticMatrix n@ value from a
-- raw matrix. It checks that the matrix has shape @n x n@ and that every row is
-- a probability vector of length @n@.
mkStochasticMatrix
  :: forall n
   . KnownNat n
  => Matrix Double
  -> Either ValidationError (StochasticMatrix n)
mkStochasticMatrix matrix
  | not (isSquare matrix) =
      Left
        NonSquareMatrix
          { rowCount = rows matrix
          , colCount = cols matrix
          }
  | actualRows /= expectedSize || actualCols /= expectedSize =
      Left
        MatrixDimensionMismatch
          { expectedSize = expectedSize
          , actualRows = actualRows
          , actualCols = actualCols
          }
  | otherwise =
      validateRows (Proxy @n) 0 (toRows matrix) *> Right (StochasticMatrix matrix)
 where
  expectedSize :: Int
  expectedSize =
    fromIntegral (natVal (Proxy @n))

  actualRows :: Int
  actualRows =
    rows matrix

  actualCols :: Int
  actualCols =
    cols matrix

-- | Multiply two stochastic matrices of the same type-level size.
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
-- library represents stochasticity up to the probability-vector tolerance, so
-- the invariant is understood numerically rather than exactly.
mulStochasticMatrix :: StochasticMatrix n -> StochasticMatrix n -> StochasticMatrix n
mulStochasticMatrix a b =
  StochasticMatrix (unStochasticMatrix a LA.<> unStochasticMatrix b)

-- | Compare two stochastic matrices approximately.
--
-- This is intentionally explicit instead of deriving 'Eq', because exact
-- structural equality on 'Double' probability data is misleading.
approxStochasticMatrixEq :: Double -> StochasticMatrix n -> StochasticMatrix n -> Bool
approxStochasticMatrixEq tolerance a b =
  and
    ( zipWith
        (\x y -> abs (x - y) <= tolerance)
        (LA.toList (LA.flatten (unStochasticMatrix a)))
        (LA.toList (LA.flatten (unStochasticMatrix b)))
    )

isSquare :: Matrix Double -> Bool
isSquare matrix =
  rows matrix == cols matrix

validateRows :: forall n . KnownNat n => Proxy n -> Int -> [LA.Vector Double] -> Either ValidationError ()
validateRows _ _ [] =
  Right ()
validateRows proxy rowIndex (rowVector : rowVectors) =
  mkProbabilityVectorAt @n rowIndex rowVector
    *> validateRows proxy (rowIndex + 1) rowVectors