{-# LANGUAGE DataKinds #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.ProbabilityVector
  ( ProbabilityVector
  , unProbabilityVector
  , mkProbabilityVector
  , mkProbabilityVectorAt
  , probabilityTolerance
  ) where

import Data.Proxy (Proxy (Proxy))
import Dtmc.ValidationError (ValidationError (..))
import GHC.TypeNats (KnownNat, Nat, natVal)
import Numeric.LinearAlgebra (Vector, size, toList)

-- | A numerically valid probability vector of length @n@.
--
-- A vector is accepted as a probability vector when:
--
-- * its runtime length is equal to the type-level dimension @n@;
-- * every entry is at least @-probabilityTolerance@;
-- * every entry is at most @1 + probabilityTolerance@;
-- * the entries sum to @1@ up to @probabilityTolerance@.
--
-- Therefore, a value of type @ProbabilityVector n@ does not mean that the
-- entries sum to exactly @1@. It means the vector passed numerical probability
-- validation at a constructor boundary.
newtype ProbabilityVector (n :: Nat) = ProbabilityVector
  { unProbabilityVector :: Vector Double
  }
  deriving (Show)

-- | Numerical tolerance used when checking probability vectors.
--
-- The tolerance should be larger than ordinary floating-point normalisation
-- error, but much smaller than meaningful modelling error.
probabilityTolerance :: Double
probabilityTolerance = 1e-9

-- | Smart constructor for probability vectors.
--
-- This version uses row index @0@ in error messages. Matrix validation should
-- usually use 'mkProbabilityVectorAt' so the offending matrix row can be
-- reported.
mkProbabilityVector :: forall n . KnownNat n => Vector Double -> Either ValidationError (ProbabilityVector n)
mkProbabilityVector =
  mkProbabilityVectorAt 0

-- | Smart constructor for probability vectors with an explicit row index.
--
-- The row index is used only for informative error reporting when a probability
-- vector is being checked as a row of a stochastic matrix.
mkProbabilityVectorAt :: forall n . KnownNat n => Int -> Vector Double -> Either ValidationError (ProbabilityVector n)
mkProbabilityVectorAt rowIndex vector
  | actualLength /= expectedLength =
      Left
        VectorDimensionMismatch
          { expectedLength = expectedLength
          , actualLength = actualLength
          }
  | otherwise =
      case firstBadEntry rowIndex (toList vector) of
        Just err ->
          Left err
        Nothing ->
          let total = sum (toList vector)
           in if abs (total - 1.0) <= probabilityTolerance
                then Right (ProbabilityVector vector)
                else Left (RowSumOffBy rowIndex total)
 where
  expectedLength :: Int
  expectedLength =
    fromIntegral (natVal (Proxy @n))

  actualLength :: Int
  actualLength =
    size vector

firstBadEntry :: Int -> [Double] -> Maybe ValidationError
firstBadEntry rowIndex entries =
  firstBadEntryFrom rowIndex 0 entries

firstBadEntryFrom :: Int -> Int -> [Double] -> Maybe ValidationError
firstBadEntryFrom _ _ [] =
  Nothing
firstBadEntryFrom rowIndex colIndex (x : xs)
  | x < -probabilityTolerance =
      Just (NegativeEntry rowIndex colIndex x)
  | x > 1.0 + probabilityTolerance =
      Just (EntryAboveOne rowIndex colIndex x)
  | otherwise =
      firstBadEntryFrom rowIndex (colIndex + 1) xs