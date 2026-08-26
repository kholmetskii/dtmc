{- |
Module      : Dtmc.Analysis.ReturnTime
Description : Exact, bounded, eventual, and expected first-return times.

First-return quantities for DTMCs. For state @i@,
@T_i = inf { t >= 1 | X_t = i }@, so time zero is never a return. Scalar
exact-time and bounded queries work through any locally finite 'Transition'.
All-state, eventual, and expected queries use a finite 'TransitionMatrix'.

Exact-time and strictly bounded queries use finite recurrences. Eventual and
expected queries use support classification and checked 'Double' linear
solves. Results are not clamped or renormalised.
-}
module Dtmc.Analysis.ReturnTime (
    LinearSystemError (..),
    MeanTime (..),
    returnTimeProbabilitiesAt,
    returnTimeProbabilityAt,
    returnTimeProbabilitiesBefore,
    returnTimeProbabilityBefore,
    returnProbabilities,
    returnProbability,
    expectedReturnTimes,
    expectedReturnTime,
) where

import Dtmc.Analysis.ReturnTime.Internal (
    LinearSystemError (..),
    MeanTime (..),
    expectedReturnTime,
    expectedReturnTimes,
    returnProbabilities,
    returnProbability,
    returnTimeProbabilitiesAt,
    returnTimeProbabilitiesBefore,
    returnTimeProbabilityAt,
    returnTimeProbabilityBefore,
 )
