{- |
Module      : Dtmc.Simplex.Internal
Description : Construction and repair of probability-simplex values.

Shared simplex construction for distribution and transition-matrix smart
constructors. Accepted values are made canonical by clamping tolerated bound
error and normalising the repaired total.
-}
module Dtmc.Simplex.Internal (
    simplexTolerance,
    canonicaliseSimplex,
    canonicaliseSimplexEntries,
) where

import Dtmc.Simplex (
    SimplexError (..),
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

-- | Absolute tolerance shared by simplex construction and sampling repair.
simplexTolerance :: Double
simplexTolerance = 1e-9

{- | Construct a canonical simplex vector when every coordinate is in
@[-simplexTolerance, 1 + simplexTolerance]@ and its total is in
@[1 - simplexTolerance, 1 + simplexTolerance]@. Tolerated negative coordinates
are clamped to zero, tolerated coordinates above one are clamped to one, and
the repaired coordinates are divided by their computed total.

Reports the first non-finite or bound error before checking the total. An
empty vector yields @Left (SumOffBy 0)@.

Time: @O(n)@. Space: @O(n)@ for dynamic-vector conversion.
-}
canonicaliseSimplex :: (KnownNat n) => S.R n -> Either SimplexError (S.R n)
canonicaliseSimplex vector =
    S.vector <$> canonicaliseSimplexEntries (LA.toList (S.extract vector))

{- | Construct a canonical finite list with the same tolerance, repair, and
error ordering as 'canonicaliseSimplex'. Entry indices refer to the supplied
list order.

An empty list yields @Left (SumOffBy 0)@. Time and result space: @O(n)@.
-}
canonicaliseSimplexEntries :: [Double] -> Either SimplexError [Double]
canonicaliseSimplexEntries entries =
    case firstInvalidEntry 0 entries of
        Just err -> Left err
        Nothing
            | abs (total - 1.0) <= simplexTolerance ->
                Right (map (/ repairedTotal) repaired)
            | otherwise -> Left (SumOffBy total)
  where
    total = foldl' (+) 0 entries
    repaired = map repair entries
    repairedTotal = foldl' (+) 0 repaired

    repair entry
        | entry < 0 = 0
        | entry > 1 = 1
        | otherwise = entry

-- Scan separately so a coordinate error reports its index before the total.
firstInvalidEntry :: Int -> [Double] -> Maybe SimplexError
firstInvalidEntry _ [] = Nothing
firstInvalidEntry index (entry : rest)
    | isNaN entry || isInfinite entry =
        Just (NonFiniteEntry index)
    | entry < negate simplexTolerance =
        Just (NegativeEntry index entry)
    | entry > 1.0 + simplexTolerance =
        Just (EntryAboveOne index entry)
    | otherwise =
        firstInvalidEntry (index + 1) rest
