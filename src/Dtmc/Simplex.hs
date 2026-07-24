{- |
Module      : Dtmc.Simplex
Description : The probability-simplex validation error.

Errors shared by the distribution and transition-matrix smart constructors.
Validation uses an absolute tolerance of @1e-9@ for coordinates and the total.
-}
module Dtmc.Simplex (
    SimplexError (..),
) where

{- | Why a vector failed to be a probability distribution. Coordinate errors
carry a zero-based index and the offending value; 'SumOffBy' carries the
computed total.

The first coordinate error takes precedence over the total. Coordinate bounds
are inclusive; the total succeeds when @abs (total - 1) <= 1e-9@. An empty
vector yields @SumOffBy 0@; with no coordinate bound error, a @NaN@ coordinate
yields @SumOffBy NaN@.
-}
data SimplexError
    = -- | Coordinate less than @-1e-9@.
      NegativeEntry Int Double
    | -- | Coordinate greater than @1 + 1e-9@.
      EntryAboveOne Int Double
    | -- | No coordinate error, but the total is outside tolerance.
      SumOffBy Double
    deriving (Eq, Show)
