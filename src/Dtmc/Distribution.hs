{- |
Module      : Dtmc.Distribution
Description : Probability distributions over a finite state space.

Probability distributions over @n@ states, used as initial and marginal laws
of a DTMC. 'mkDistribution' checks the simplex invariant with the @1e-9@
tolerance documented by 'SimplexError'.
-}
module Dtmc.Distribution (
    Distribution,
    DistributionError (..),
    mkDistribution,
    unDistribution,
) where

import Data.Bifunctor (
    first,
 )
import Dtmc.Distribution.Internal (
    Distribution (Distribution),
    unDistribution,
 )
import Dtmc.Simplex (
    SimplexError,
 )
import Dtmc.Simplex.Internal (
    validateSimplex,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra.Static qualified as S

{- | A simplex failure while constructing a distribution. The wrapper keeps it
distinct from a transition-matrix row failure.
-}
newtype DistributionError
    = -- | Wrap the underlying simplex failure.
      DistributionError SimplexError
    deriving (Eq, Show)

{- | Validate a raw state distribution. On success the vector is preserved
exactly; tolerated coordinates are not clamped and the total is not
renormalised. For @n = 0@, returns
@Left (DistributionError (SumOffBy 0))@.

Time: @O(n)@. Space: @O(n)@ for validation.
-}
mkDistribution :: (KnownNat n) => S.R n -> Either DistributionError (Distribution n)
mkDistribution vector =
    Distribution vector <$ first DistributionError (validateSimplex vector)
