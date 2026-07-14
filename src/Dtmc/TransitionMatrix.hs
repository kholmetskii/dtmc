-- |
-- Module      : Dtmc.TransitionMatrix
-- Description : Row-stochastic one-step transition matrices.
--
-- Public interface for 'TransitionMatrix', the one-step dynamics of a
-- discrete-time Markov chain on @n@ states. 'mkTransitionMatrix' enforces the
-- row-stochastic invariant (each row is a 'Distribution' over next states);
-- composition of steps is inherited from the 'Semigroup' instance in
-- "Dtmc.Internal.Types" and exposed here as 'mulTransitionMatrix'.
module Dtmc.TransitionMatrix (
    TransitionMatrix,
    TransitionError (..),
    mkTransitionMatrix,
    unTransitionMatrix,
    mulTransitionMatrix,
    rowAt,
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
import Numeric.LinearAlgebra.Static qualified as S

-- | A row of the matrix was not a valid distribution: carries the zero-based
-- row index together with the underlying simplex failure.
data TransitionError = InRow Int SimplexError
    deriving (Eq, Show)

-- | Smart constructor: accept a raw matrix only if every row passes
-- 'validateSimplex', reporting the first offending row otherwise. The sole
-- sanctioned way to build a 'TransitionMatrix'.
mkTransitionMatrix :: (KnownNat n) => S.Sq n -> Either TransitionError (TransitionMatrix n)
mkTransitionMatrix matrix =
    TransitionMatrix{unTransitionMatrix = matrix} <$ traverse_ validateRow (zip [0 ..] (S.toRows matrix))
  where
    validateRow (index, row) =
        first (InRow index) (validateSimplex row)

-- | Compose two steps by matrix multiplication (the 'Semigroup' '<>'). No
-- re-validation is needed: the product of two row-stochastic matrices is
-- row-stochastic, so the result is valid by construction.
mulTransitionMatrix :: (KnownNat n) => TransitionMatrix n -> TransitionMatrix n -> TransitionMatrix n
mulTransitionMatrix = (<>)

-- | The @i@-th row as a 'Distribution': the conditional law of the next state
-- given the chain is currently in state @i@. The 'Finite' index keeps @i@
-- statically in range, so the lookup is total.
rowAt :: (KnownNat n) => TransitionMatrix n -> Finite n -> Distribution n
rowAt (TransitionMatrix matrix) index =
    Distribution{unDistribution = S.toRows matrix !! fromIntegral (getFinite index)}
