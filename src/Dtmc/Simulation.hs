{- |
Module      : Dtmc.Simulation
Description : Sampling states and running the chain forward.

Random sampling from state distributions and transition probabilities. Both
operations mutate the supplied MWC generator in any 'PrimMonad'. Passing each
sampled state back to 'step' while reusing the generator produces a trajectory.
-}
module Dtmc.Simulation (
    sampleFrom,
    step,
) where

import Control.Monad.Primitive (
    PrimMonad,
    PrimState,
 )
import Data.Finite (
    Finite,
    finite,
 )
import Dtmc.Distribution.Internal (
    Distribution,
    unDistribution,
 )
import Dtmc.Simplex.Internal (
    snapToSimplex,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    rowAt,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra.Static qualified as S
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

{- | Sample one transition from the supplied state using its stored matrix row.
The 'Finite' state index is always in range. Passing each result back with the
same generator advances one trajectory.

Rows produced by unchecked or floating-point matrix arithmetic inherit the
repair and error behavior of 'sampleFrom'. Sampling itself takes @O(n)@ time
and temporary space.
-}
step ::
    (KnownNat n, PrimMonad m) =>
    TransitionMatrix n ->
    Finite n ->
    MWC.Gen (PrimState m) ->
    m (Finite n)
step matrix state =
    sampleFrom (rowAt matrix state)
