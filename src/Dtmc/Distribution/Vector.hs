{- |
Module      : Dtmc.Distribution.Vector
Description : Dense probability vectors over type-indexed finite state spaces.

t'DistributionVector' stores a probability law over @{0 .. n-1}@ in a static
vector. 'mkDistributionVector' checks the simplex invariant with the @1e-9@
tolerance documented by 'Dtmc.Simplex.SimplexError'.
-}
module Dtmc.Distribution.Vector (
    DistributionVector,
    mkDistributionVector,
    unDistributionVector,
) where

import Data.Bifunctor (
    first,
 )
import Dtmc.Distribution (
    DistributionError (DistributionError),
 )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector (DistributionVector),
    unDistributionVector,
 )
import Dtmc.Simplex.Internal (
    validateSimplex,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra.Static qualified as S

{- | Validate a raw state distribution vector. On success it is preserved
exactly; tolerated coordinates are not clamped and the total is not
renormalised. For @n = 0@, returns
@Left (DistributionError (SumOffBy 0))@.

Time: @O(n)@. Space: @O(n)@ for validation.
-}
mkDistributionVector ::
    (KnownNat n) =>
    S.R n ->
    Either DistributionError (DistributionVector n)
mkDistributionVector vector =
    DistributionVector vector <$ first DistributionError (validateSimplex vector)
