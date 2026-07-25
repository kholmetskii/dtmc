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
    probabilityAt,
) where

import Data.Bifunctor (
    first,
 )
import Data.Finite (
    Finite,
    getFinite,
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
import Numeric.LinearAlgebra qualified as LA
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

{- | Read the stored probability @mu(j)@ of a single state.

Mathematically this is the coordinate of the law at index @j@: for an initial
or marginal distribution @mu@, @probabilityAt mu j = mu(j) = P(X = j)@. The
bounded 'Finite' @n@ index makes the lookup total, so there is no
out-of-range or partial case.

The stored coordinate is returned exactly as held. No clamping to @[0, 1]@,
renormalisation, or revalidation is performed, so a value carried through
tolerated construction or floating-point arithmetic is reported unchanged and
may lie slightly outside @[0, 1]@.

Time: @O(n)@ to materialise the underlying vector; the index itself is @O(1)@.
-}
probabilityAt :: (KnownNat n) => Distribution n -> Finite n -> Double
probabilityAt distribution index =
    S.extract (unDistribution distribution)
        `LA.atIndex` fromIntegral (getFinite index)
