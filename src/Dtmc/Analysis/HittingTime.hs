{- |
Module      : Dtmc.Analysis.HittingTime
Description : Exact, bounded, eventual, competing, and expected hitting times.

Hitting-time quantities for DTMCs. Scalar exact-time and bounded queries work
through any locally finite 'Transition', including kernels on infinite state
spaces. All-state, eventual, competing, and expected queries use a finite
'TransitionMatrix'. For a target set @A@,
@H_A = inf { t >= 0 | X_t in A }@.

Exact-time and strictly bounded queries use finite recurrences. Eventual,
competing, and expected queries use support reachability and checked 'Double'
linear solves. Results are not clamped or renormalised.
-}
module Dtmc.Analysis.HittingTime (
    LinearSystemError (..),
    MeanTime (..),
    hittingTimeProbabilitiesAt,
    hittingTimeProbabilityAt,
    hittingTimeProbabilitiesBefore,
    hittingTimeProbabilityBefore,
    hittingProbabilities,
    hittingProbability,
    hittingBeforeProbabilities,
    hittingBeforeProbability,
    expectedHittingTimes,
    expectedHittingTime,
) where

import Dtmc.Analysis.HittingTime.Internal (
    LinearSystemError (..),
    MeanTime (..),
    expectedHittingTime,
    expectedHittingTimes,
    hittingBeforeProbabilities,
    hittingBeforeProbability,
    hittingProbabilities,
    hittingProbability,
    hittingTimeProbabilitiesAt,
    hittingTimeProbabilitiesBefore,
    hittingTimeProbabilityAt,
    hittingTimeProbabilityBefore,
 )
