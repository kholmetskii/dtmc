-- |
-- Module      : Dtmc.Internal.Simplex
-- Description : Validation of the probability-simplex invariant.
--
-- Shared numeric check used by both the distribution and transition-matrix smart
-- constructors: it decides whether a raw vector is a valid probability
-- distribution (a point on the standard simplex). All comparisons allow
-- 'simplexTolerance' of slack so floating-point values are not spuriously
-- rejected.
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

-- | Why a vector failed to be a probability distribution. Indices are
-- zero-based and each 'Double' echoes the offending value for diagnostics.
data SimplexError
    -- | Entry at this index is negative beyond tolerance.
    = NegativeEntry Int Double
    -- | Entry at this index exceeds one beyond tolerance.
    | EntryAboveOne Int Double
    -- | Entries were individually in range but summed to this value, not one.
    | SumOffBy Double
    deriving (Eq, Show)

-- | Absolute slack allowed on every simplex comparison (non-negativity, the
-- upper bound of one, and the unit-sum condition). Small enough to catch real
-- modelling errors while tolerating rounding noise.
simplexTolerance :: Double
simplexTolerance = 1e-9

-- | Check that @vector@ is a probability distribution: scan for the first
-- out-of-range coordinate, and if none is found require the total to be one
-- within 'simplexTolerance'. The per-entry scan runs first so a malformed
-- coordinate is reported at its index rather than masked by the sum check.
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

-- | Walk the coordinates left to right, returning the first that is negative
-- or greater than one (each beyond 'simplexTolerance'), tagged with its index.
-- 'Nothing' means every coordinate is individually within @[0,1]@.
firstInvalidEntry :: Int -> [Double] -> Maybe SimplexError
firstInvalidEntry _ [] = Nothing
firstInvalidEntry index (entry : rest)
    | entry < negate simplexTolerance =
        Just (NegativeEntry index entry)
    | entry > 1.0 + simplexTolerance =
        Just (EntryAboveOne index entry)
    | otherwise =
        firstInvalidEntry (index + 1) rest
