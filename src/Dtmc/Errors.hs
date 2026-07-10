-- |
-- Module      : Dtmc.Errors
--
-- The shared error vocabulary.
--
-- 'SimplexError' is the common core: a violation of Δ^{n-1}, naming the
-- coordinate and the value. The two carriers wrap it with their own context —
-- a standalone distribution has no rows, a transition matrix does. Hence
-- 'DistributionError' carries no row index and 'TransitionError' does.
--
-- Kept in @other-modules@: users obtain these constructors by re-export from
-- the module that produces them ('Dtmc.Distribution', 'Dtmc.TransitionMatrix'),
-- so 'Dtmc.Errors' is a vocabulary, not a public entry point. See
-- docs/DECISIONS.md D10.
module Dtmc.Errors
  ( SimplexError (..)
  , DistributionError (..)
  , TransitionError (..)
  ) where

-- | Why a point fails to lie in Δ^{n-1}. 'Int' is the coordinate index.
data SimplexError
  = NegativeEntry Int Double
  | EntryAboveOne Int Double
  | SumOffBy Double
  deriving (Eq, Show)

-- | No row index: a standalone distribution has no rows.
newtype DistributionError = DistributionError SimplexError
  deriving (Eq, Show)

-- | @InRow 2 (SumOffBy 1.03)@ reads: "row 2 violates the distribution invariant".
--
-- @data@, not @newtype@: two fields. There are no dimension constructors here —
-- @NonSquareMatrix@ and @MatrixDimensionMismatch@ became unrepresentable once
-- the size moved into the type (docs/DECISIONS.md D1).
data TransitionError = InRow Int SimplexError
  deriving (Eq, Show)
