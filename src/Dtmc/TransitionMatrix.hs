-- |
-- Module      : Dtmc.TransitionMatrix
-- Description : Row-stochastic transition matrices and their monoid.
--
-- Public interface for t'TransitionMatrix', the one-step law of a discrete-time
-- Markov chain on @n@ states. 'mkTransitionMatrix' enforces the row-stochastic
-- invariant (each row is a t'Distribution' over next states). The one-step
-- matrix generates a monoid of @k@-step transitions, all backed by the 'Monoid'
-- instance in "Dtmc.Internal.Types": composition ('mulTransitionMatrix', the
-- '<>' product), the zero-step identity ('identityMatrix', 'mempty'), and
-- @k@-step powers ('matrixPower').
module Dtmc.TransitionMatrix (
    TransitionMatrix,
    TransitionMatrixError (..),
    mkTransitionMatrix,
    unTransitionMatrix,
    mulTransitionMatrix,
    identityMatrix,
    matrixPower,
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
import Data.Semigroup (
    mtimesDefault,
 )
import Dtmc.Internal.Simplex (
    validateSimplex,
 )
import Dtmc.Internal.Types (
    Distribution (..),
    TransitionMatrix (..),
    unsafeTransitionMatrix,
 )
import Dtmc.Simplex (
    SimplexError,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )

-- | A row of the matrix was not a valid distribution: carries the zero-based
-- row index together with the underlying simplex failure.
data TransitionMatrixError = InRow Int SimplexError
    deriving (Eq, Show)

-- | Smart constructor: accept a raw matrix only if every row passes
-- @validateSimplex@, reporting the first offending row otherwise. The sole
-- sanctioned way to build a t'TransitionMatrix'.
mkTransitionMatrix :: (KnownNat n) => S.Sq n -> Either TransitionMatrixError (TransitionMatrix n)
mkTransitionMatrix matrix =
    unsafeTransitionMatrix matrix <$ traverse_ validateRow (zip [0 ..] (S.toRows matrix))
  where
    validateRow (index, row) =
        first (InRow index) (validateSimplex row)

-- | Compose two steps by matrix multiplication (the 'Semigroup' '<>'). No
-- re-validation is needed: the product of two row-stochastic matrices is
-- row-stochastic, so the result is valid by construction.
mulTransitionMatrix :: (KnownNat n) => TransitionMatrix n -> TransitionMatrix n -> TransitionMatrix n
mulTransitionMatrix = (<>)

-- | The @n*n@ identity as a transition matrix: the zero-step transition
-- (@mempty@), which leaves any distribution unchanged.
identityMatrix :: (KnownNat n) => TransitionMatrix n
identityMatrix = mempty

-- | The @k@-step transition matrix @p^k@, the @k@-fold monoidal product
-- (@matrixPower 0 = identityMatrix@). This is Chapman--Kolmogorov in matrix
-- form: @matrixPower (m + n) p == matrixPower m p <> matrixPower n p@.
matrixPower :: (KnownNat n) => Natural -> TransitionMatrix n -> TransitionMatrix n
matrixPower = mtimesDefault

-- | The @i@-th row as a t'Distribution': the conditional law of the next state
-- given the chain is currently in state @i@. The 'Finite' index keeps @i@
-- statically in range, so the lookup is total.
rowAt :: (KnownNat n) => TransitionMatrix n -> Finite n -> Distribution n
rowAt p index =
    Distribution{unDistribution = S.toRows (unTransitionMatrix p) !! fromIntegral (getFinite index)}
