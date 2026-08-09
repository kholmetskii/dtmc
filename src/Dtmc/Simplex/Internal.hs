{- |
Module      : Dtmc.Simplex.Internal
Description : Validation of the probability-simplex invariant.

Shared simplex validation for distribution and transition-matrix smart
constructors. It also provides the small-negative repair used before
categorical sampling. Successful validation preserves the input; it does not
clamp or renormalise coordinates.
-}
module Dtmc.Simplex.Internal (
    validateSimplex,
    validateSimplexEntries,
    snapToSimplex,
) where

import Dtmc.Simplex (
    SimplexError (..),
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

-- Keep validation and sampling repair on one private absolute threshold.
tolerance :: Double
tolerance = 1e-9

{- | Accept a vector iff every coordinate is in
@[-tolerance, 1 + tolerance]@ and its total is in
@[1 - tolerance, 1 + tolerance]@. Reports the first coordinate error before
checking the total and does not modify accepted values.

An empty vector yields @Left (SumOffBy 0)@. If the coordinate scan finds no
bound violation, @NaN@ yields @Left (SumOffBy NaN)@; infinities fail their
coordinate bound.

Time: @O(n)@. Space: @O(n)@ for dynamic-vector conversion.
-}
validateSimplex :: (KnownNat n) => S.R n -> Either SimplexError ()
validateSimplex =
    validateSimplexEntries . LA.toList . S.extract

{- | Validate a finite list with the same tolerance and error ordering as
'validateSimplex'. Entry indices refer to the supplied list order. The input
is preserved by callers; this function performs no clamping or
renormalisation.

An empty list yields @Left (SumOffBy 0)@. Time: @O(n)@; space: @O(1)@ beyond
the supplied list.
-}
validateSimplexEntries :: [Double] -> Either SimplexError ()
validateSimplexEntries entries =
    case firstInvalidEntry 0 entries of
        Just err -> Left err
        Nothing
            | abs (total - 1.0) <= tolerance -> Right ()
            | otherwise -> Left (SumOffBy total)
  where
    total = foldl' (+) 0 entries

-- Scan separately so a coordinate error reports its index before the total.
firstInvalidEntry :: Int -> [Double] -> Maybe SimplexError
firstInvalidEntry _ [] = Nothing
firstInvalidEntry index (entry : rest)
    | entry < negate tolerance =
        Just (NegativeEntry index entry)
    | entry > 1.0 + tolerance =
        Just (EntryAboveOne index entry)
    | otherwise =
        firstInvalidEntry (index + 1) rest

{- | Replace coordinates in @[-tolerance, 0)@ with zero before categorical
sampling. Non-negative values are unchanged; the result is not renormalised
or clamped above one.

Raises an error for a coordinate below @-tolerance@ or a @NaN@ coordinate.
An empty vector remains empty.

Time and result space: @O(n)@.
-}
snapToSimplex :: LA.Vector Double -> LA.Vector Double
snapToSimplex =
    LA.cmap snap
  where
    snap value
        | value >= 0 = value
        | value >= negate tolerance = 0
        | otherwise =
            error
                ( "Dtmc.Simplex.Internal.snapToSimplex: probability coordinate "
                    <> show value
                    <> " is below -tolerance"
                )
