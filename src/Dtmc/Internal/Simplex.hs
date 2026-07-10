module Dtmc.Internal.Simplex (
    SimplexError (..),
    simplexTolerance,
    validateSimplex,
) where

import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

data SimplexError
    = NegativeEntry Int Double
    | EntryAboveOne Int Double
    | SumOffBy Double
    deriving (Eq, Show)

simplexTolerance :: Double
simplexTolerance = 1e-9

validateSimplex :: (KnownNat n) => S.R n -> Either SimplexError ()
validateSimplex vector =
    case firstInvalidEntry 0 entries of
        Just err -> Left err
        Nothing
            | abs (total - 1.0) <= simplexTolerance -> Right ()
            | otherwise -> Left (SumOffBy total)
  where
    entries = LA.toList (S.extract vector)
    total = sum entries

firstInvalidEntry :: Int -> [Double] -> Maybe SimplexError
firstInvalidEntry _ [] = Nothing
firstInvalidEntry index (entry : rest)
    | entry < negate simplexTolerance =
        Just (NegativeEntry index entry)
    | entry > 1.0 + simplexTolerance =
        Just (EntryAboveOne index entry)
    | otherwise =
        firstInvalidEntry (index + 1) rest
