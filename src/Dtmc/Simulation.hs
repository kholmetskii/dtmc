{- |
Module      : Dtmc.Simulation
Description : Sampling states and running the chain forward.

Random sampling from dense or sparse state distributions, plus shared
simulation through any locally finite 'MarkovKernel'. All operations mutate
the supplied MWC generator in any 'PrimMonad'.
-}
module Dtmc.Simulation (
    sampleFrom,
    sampleSparseFrom,
    step,
    simulateN,
) where

import Control.Monad.Primitive (
    PrimMonad,
    PrimState,
 )
import Data.Finite (
    Finite,
    finite,
 )
import Dtmc.Distribution (
    SparseDistribution,
    sparseEntries,
 )
import Dtmc.Distribution.Internal (
    Distribution,
    unDistribution,
 )
import Dtmc.Kernel (
    MarkovKernel (..),
 )
import Dtmc.Simplex.Internal (
    snapToSimplex,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )
import System.Random.MWC qualified as MWC
import System.Random.MWC.Distributions qualified as MWCD

{- | Draw a state in @{0 .. n-1}@ with probability proportional to the stored
weights. Before sampling, coordinates in @[-1e-9, 0)@ are replaced by zero;
the categorical sampler scales by the resulting total, so no explicit
renormalisation is stored.

Requires a non-empty distribution. A coordinate below @-1e-9@ or a @NaN@
coordinate raises an error during repair; the backend may also reject invalid
unchecked weights. The returned 'Finite' index is in range.

Time and temporary space: @O(n)@ per draw.
-}
sampleFrom ::
    (KnownNat n, PrimMonad m) =>
    Distribution n ->
    MWC.Gen (PrimState m) ->
    m (Finite n)
sampleFrom distribution generator = do
    index <- MWCD.categorical weights generator
    pure (finite (fromIntegral index))
  where
    weights =
        snapToSimplex
            (S.extract (unDistribution distribution))

{- | Sample one transition from a state through any 'MarkovKernel'. Passing
each result back with the same generator advances one trajectory. The finite
transition law inherits the repair and error behavior of 'sampleSparseFrom'.
-}
step ::
    (PrimMonad m, MarkovKernel kernel) =>
    kernel ->
    KernelState kernel ->
    MWC.Gen (PrimState m) ->
    m (KernelState kernel)
step kernel state =
    sampleSparseFrom (transitionLaw kernel state)

{- | Draw one state from a validated sparse law. Tolerated negative weights
are snapped to zero with the same rule as 'sampleFrom'.

Time and temporary space: @O(s)@ for stored support size @s@.
-}
sampleSparseFrom ::
    (PrimMonad m) =>
    SparseDistribution state ->
    MWC.Gen (PrimState m) ->
    m state
sampleSparseFrom distribution generator = do
    index <- MWCD.categorical weights generator
    pure (states !! index)
  where
    entries = sparseEntries distribution
    states = map fst entries
    weights = snapToSimplex (LA.fromList (map snd entries))

{- | Simulate exactly @k@ transitions through any 'MarkovKernel' and return
the trajectory including its initial state. The result has length @k + 1@.
-}
simulateN ::
    (PrimMonad m, MarkovKernel kernel) =>
    Natural ->
    kernel ->
    KernelState kernel ->
    MWC.Gen (PrimState m) ->
    m [KernelState kernel]
simulateN transitions kernel initial generator =
    go transitions initial [initial]
  where
    go 0 _ reversed = pure (reverse reversed)
    go remaining current reversed = do
        next <- step kernel current generator
        go (remaining - 1) next (next : reversed)
