{- |
Module      : Dtmc.Analysis
Description : Curated facade for DTMC analysis.

Re-exports fixed-time probability queries, hitting and return times,
classification, and stationary distributions. Focused modules remain
available when a component should depend on only one mathematical subject.
-}
module Dtmc.Analysis (
    module FixedTime,
    module HittingTime,
    module ReturnTime,
    module Classification,
    module Stationary,
) where

import Dtmc.Analysis.Classification as Classification
import Dtmc.Analysis.FixedTime as FixedTime
import Dtmc.Analysis.HittingTime as HittingTime
import Dtmc.Analysis.ReturnTime as ReturnTime hiding (
    LinearSystemError (..),
    MeanTime (..),
 )
import Dtmc.Analysis.Stationary as Stationary hiding (
    LinearSystemError (..),
 )
