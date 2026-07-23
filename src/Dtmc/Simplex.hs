-- |
-- Module      : Dtmc.Simplex
-- Description : The probability-simplex validation error.
--
-- The public error describing how a raw vector fails to be a probability
-- distribution (a point on the standard simplex). It is embedded in both
-- 'Dtmc.Distribution.DistributionError' and
-- 'Dtmc.TransitionMatrix.TransitionMatrixError', so it is part of the public
-- API even though the validation that produces it lives in the internal
-- "Dtmc.Internal.Simplex".
module Dtmc.Simplex (
    SimplexError (..),
) where

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
