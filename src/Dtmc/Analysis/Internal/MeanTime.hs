{- |
Module      : Dtmc.Analysis.Internal.MeanTime
Description : Internal carrier for expected-time analysis results.

The shared 'MeanTime' result type for hitting- and return-time expectations. Its constructors are exposed
so internal solvers can build results directly. 'FiniteMean' performs no
validation: callers can construct negative, non-finite, or @NaN@ values, and
library results reach 'InfiniteMean' from support-graph reachability rather
than from floating-point overflow.
-}
module Dtmc.Analysis.Internal.MeanTime (
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
