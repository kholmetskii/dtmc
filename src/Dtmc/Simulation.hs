{- |
Module      : Dtmc.Simulation
Description : Sampling states and running the chain forward.

Random sampling from dense or sparse state distributions, plus shared
simulation through any locally finite 'Transition'. All operations mutate
the supplied MWC generator in any 'PrimMonad'.
-}
module Dtmc.Simulation (
    sample,
    step,
    simulateN,
) where

import Control.Monad.Primitive (
    PrimMonad,
    PrimState,
 )
import Dtmc.Distribution (
    Distribution (..),
 )
import Dtmc.Kernel (
    Transition (..),
 )
import Dtmc.Simplex.Internal (
    snapToSimplex,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.Natural (
    Natural,
 )
import System.Random.MWC qualified as MWC
import System.Random.MWC.Distributions qualified as MWCD

{- | Draw a state from any finite-support 'Distribution'. Before sampling,
stored weights in @[-1e-9, 0)@ are replaced by zero; the categorical sampler
scales by the resulting total, so no explicit renormalisation is stored.

Requires a non-empty distribution. A weight below @-1e-9@ or a @NaN@ raises
an error during repair; the backend may also reject invalid unchecked weights.

Time and temporary space: @O(s)@ for stored support size @s@.
-}
sample ::
    (Distribution distribution, PrimMonad m) =>
    distribution ->
    MWC.Gen (PrimState m) ->
    m (DistributionState distribution)
sample distribution generator = do
    index <- MWCD.categorical weights generator
    pure (states !! index)
  where
    entries = distributionWeights distribution
    states = map fst entries
    weights = snapToSimplex (LA.fromList (map snd entries))

{- | Sample one transition from a state through any 'Transition'. Passing
each result back with the same generator advances one trajectory. The finite
transition law inherits the repair and error behavior of 'sample'.
-}
step ::
    (PrimMonad m, Transition kernel) =>
    kernel ->
    TransitionState kernel ->
    MWC.Gen (PrimState m) ->
    m (TransitionState kernel)
step kernel state =
    sample (transitionLaw kernel state)

{- | Simulate exactly @k@ transitions through any 'Transition' and return
the trajectory including its initial state. The result has length @k + 1@.
-}
simulateN ::
    (PrimMonad m, Transition kernel) =>
    Natural ->
    kernel ->
    TransitionState kernel ->
    MWC.Gen (PrimState m) ->
    m [TransitionState kernel]
simulateN transitions kernel initial generator =
    go transitions initial [initial]
  where
    go 0 _ reversed = pure (reverse reversed)
    go remaining current reversed = do
        next <- step kernel current generator
        go (remaining - 1) next (next : reversed)
