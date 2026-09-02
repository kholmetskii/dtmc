{- |
Module      : Dtmc.Simplex
Description : Probability-simplex construction errors.

Errors shared by the distribution and transition-matrix smart constructors.
Construction uses an absolute tolerance of @1e-9@ for coordinates and the
total.
-}
module Dtmc.Simplex (
    SimplexError (..),
) where

{- | Why a vector failed to be a probability distribution. Bound errors carry
a zero-based index and the offending value; 'NonFiniteEntry' identifies a
@NaN@ or infinite coordinate; 'SumOffBy' carries the computed total.

The first coordinate error takes precedence over the total. Coordinate bounds
accept @[-1e-9, 1 + 1e-9]@; the total succeeds when
@abs (total - 1) <= 1e-9@. Smart constructors clamp accepted coordinates to
@[0, 1]@ and normalise before storage. An empty vector yields @SumOffBy 0@.
-}
data SimplexError
    = -- | Coordinate is @NaN@ or infinite.
      NonFiniteEntry Int
    | -- | Coordinate less than @-1e-9@.
      NegativeEntry Int Double
    | -- | Coordinate greater than @1 + 1e-9@.
      EntryAboveOne Int Double
    | -- | No coordinate error, but the total is outside tolerance.
      SumOffBy Double
    deriving (Eq, Show)
