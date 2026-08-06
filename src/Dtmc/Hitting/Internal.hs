{- |
Module      : Dtmc.Hitting.Internal
Description : Raw carrier for expected-time results (unsafe underbelly).

The 'MeanTime' result type behind "Dtmc.Hitting". Its constructors are exposed
so internal solvers can build results directly. 'FiniteMean' performs no
validation: callers can construct negative, non-finite, or @NaN@ values, and
library results reach 'InfiniteMean' from support-graph reachability rather
than from floating-point overflow.
-}
module Dtmc.Hitting.Internal (
    MeanTime (..),
) where

{- | An expected number of transitions, represented either by a 'Double' or an
exact infinite case. Library results use 'InfiniteMean' based on support-graph
reachability rather than floating-point overflow.

'FiniteMean' performs no validation: callers can construct negative,
non-finite, or @NaN@ values. Derived ordering places every 'FiniteMean'
constructor before 'InfiniteMean'; comparisons between finite constructors
inherit the behavior of 'Double', including @NaN@.
-}
data MeanTime
    = -- | A mathematically non-negative finite mean, subject to solver rounding.
      FiniteMean Double
    | -- | The target or return is not reached with probability one.
      InfiniteMean
    deriving (Eq, Ord, Show)
