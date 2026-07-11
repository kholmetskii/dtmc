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
    simplexTolerance,
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
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import System.Random.MWC qualified as MWC
import System.Random.MWC.Distributions qualified as MWCD

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
        sanitizeWeights
            (S.extract (unDistribution distribution))

step ::
    (KnownNat n, PrimMonad m) =>
    TransitionMatrix n ->
    Finite n ->
    MWC.Gen (PrimState m) ->
    m (Finite n)
step matrix state =
    sampleFrom (rowAt matrix state)

sanitizeWeights :: LA.Vector Double -> LA.Vector Double
sanitizeWeights =
    LA.cmap sanitize
  where
    sanitize value
        | value >= 0 = value
        | value >= negate simplexTolerance = 0
        | otherwise =
            error
                ( "Dtmc.Simulation: internal invariant violation: "
                    <> "probability coordinate "
                    <> show value
                    <> " is below -simplexTolerance"
                )
