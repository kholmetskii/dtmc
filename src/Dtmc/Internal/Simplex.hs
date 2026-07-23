-- |
-- Module      : Dtmc.Internal.Simplex
-- Description : Validation of the probability-simplex invariant.
--
-- Shared numeric check used by both the distribution and transition-matrix smart
-- constructors: it decides whether a raw vector is a valid probability
-- distribution (a point on the standard simplex). All comparisons allow the
-- shared 'tolerance' of slack (from "Dtmc.Approx") so floating-point values are
-- not spuriously rejected. It also hosts 'snapToSimplex', the floating-point
-- repair applied before categorical sampling.
module Dtmc.Internal.Simplex (
    validateSimplex,
    snapToSimplex,
) where

import Dtmc.Approx (
    tolerance,
 )
import Dtmc.Simplex (
    SimplexError (..),
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

-- | Check that @vector@ is a probability distribution: scan for the first
-- out-of-range coordinate, and if none is found require the total to be one
-- within 'tolerance'. The per-entry scan runs first so a malformed coordinate is
-- reported at its index rather than masked by the sum check.
validateSimplex :: (KnownNat n) => S.R n -> Either SimplexError ()
validateSimplex vector =
    case firstInvalidEntry 0 entries of
        Just err -> Left err
        Nothing
            | abs (total - 1.0) <= tolerance -> Right ()
            | otherwise -> Left (SumOffBy total)
  where
    entries = LA.toList (S.extract vector)
    total = sum entries

-- | Walk the coordinates left to right, returning the first that is negative
-- or greater than one (each beyond 'tolerance'), tagged with its index.
-- 'Nothing' means every coordinate is individually within @[0,1]@.
firstInvalidEntry :: Int -> [Double] -> Maybe SimplexError
firstInvalidEntry _ [] = Nothing
firstInvalidEntry index (entry : rest)
    | entry < negate tolerance =
        Just (NegativeEntry index entry)
    | entry > 1.0 + tolerance =
        Just (EntryAboveOne index entry)
    | otherwise =
        firstInvalidEntry (index + 1) rest

-- | Snap small negative coordinates -- those within 'tolerance' of zero, an
-- artefact of floating-point arithmetic -- to exactly zero, so a probability
-- vector is accepted by a categorical sampler. A coordinate more negative than
-- that signals a real invariant violation (a programmer error) and fails loudly.
snapToSimplex :: LA.Vector Double -> LA.Vector Double
snapToSimplex =
    LA.cmap snap
  where
    snap value
        | value >= 0 = value
        | value >= negate tolerance = 0
        | otherwise =
            error
                ( "Dtmc.Internal.Simplex.snapToSimplex: probability coordinate "
                    <> show value
                    <> " is below -tolerance"
                )
