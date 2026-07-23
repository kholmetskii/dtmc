-- |
-- Module      : Dtmc.Simulation
-- Description : Sampling states and running the chain forward.
--
-- Stochastic-simulation layer over the pure types. 'sampleFrom' draws a single
-- state from a t'Distribution', and 'step' advances the chain one transition by
-- sampling from the current state's row. Both run in any 'PrimMonad' with an MWC
-- generator, so threading a single 'MWC.Gen' realises a trajectory.
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
import Dtmc.Internal.Simplex (
    snapToSimplex,
 )
import Dtmc.Internal.Types (
    Distribution,
    unDistribution,
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

-- | Draw one state index from a distribution, using its entries as categorical
-- weights. The result is a valid 'Finite' @n@ because the weight vector has
-- length @n@; coordinates first pass through @snapToSimplex@ to absorb tiny
-- negative rounding artefacts the sampler would otherwise reject.
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

-- | One step of the chain: from the current @state@, sample the next state from
-- the corresponding row of the transition matrix. Iterating 'step' with the
-- same generator produces a trajectory.
step ::
    (KnownNat n, PrimMonad m) =>
    TransitionMatrix n ->
    Finite n ->
    MWC.Gen (PrimState m) ->
    m (Finite n)
step matrix state =
    sampleFrom (rowAt matrix state)
