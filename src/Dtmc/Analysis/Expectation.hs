{- |
Module      : Dtmc.Analysis.Expectation
Description : Finite and infinite expectations of non-negative quantities.

A shared result type for expectations that may be mathematically infinite.
It is used by hitting-time, return-time, and total visit-count analysis.
-}
module Dtmc.Analysis.Expectation (
    Expectation (..),
) where

{- | An expectation of a non-negative random quantity.

'FiniteExpectation' performs no validation: callers can construct negative,
non-finite, or @NaN@ values. Library functions use 'InfiniteExpectation' for a
structural mathematical infinity, not floating-point overflow. Derived
ordering places every 'FiniteExpectation' before 'InfiniteExpectation'; finite
comparisons inherit the behavior of 'Double', including @NaN@.
-}
data Expectation
    = -- | A mathematically non-negative finite expectation, subject to rounding.
      FiniteExpectation Double
    | -- | A mathematically infinite expectation.
      InfiniteExpectation
    deriving (Eq, Ord, Show)
