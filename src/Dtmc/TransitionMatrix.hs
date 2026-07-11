module Dtmc.TransitionMatrix (
    TransitionMatrix,
    TransitionError (..),
    mkTransitionMatrix,
    unTransitionMatrix,
    mulTransitionMatrix,
    rowAt,
    approxTransitionMatrixEq,
) where

import Data.Bifunctor (
    first,
 )
import Data.Finite (
    Finite,
    getFinite,
 )
import Data.Foldable (
    traverse_,
 )
import Dtmc.Internal.Simplex (
    SimplexError,
    validateSimplex,
 )
import Dtmc.Internal.Types (
    Distribution (..),
    TransitionMatrix (..),
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

data TransitionError = InRow Int SimplexError
    deriving (Eq, Show)

mkTransitionMatrix :: (KnownNat n) => S.Sq n -> Either TransitionError (TransitionMatrix n)
mkTransitionMatrix matrix =
    TransitionMatrix{unTransitionMatrix = matrix} <$ traverse_ validateRow (zip [0 ..] (S.toRows matrix))
  where
    validateRow (index, row) =
        first (InRow index) (validateSimplex row)

mulTransitionMatrix :: (KnownNat n) => TransitionMatrix n -> TransitionMatrix n -> TransitionMatrix n
mulTransitionMatrix = (<>)

rowAt :: (KnownNat n) => TransitionMatrix n -> Finite n -> Distribution n
rowAt (TransitionMatrix matrix) index =
    Distribution{unDistribution = S.toRows matrix !! fromIntegral (getFinite index)}

approxTransitionMatrixEq :: (KnownNat n) => Double -> TransitionMatrix n -> TransitionMatrix n -> Bool
approxTransitionMatrixEq tolerance (TransitionMatrix left) (TransitionMatrix right) =
    and (zipWith close (entries left) (entries right))
  where
    entries = LA.toList . LA.flatten . S.extract
    close x y = abs (x - y) <= tolerance
