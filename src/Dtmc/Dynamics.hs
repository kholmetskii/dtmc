module Dtmc.Dynamics (
    evolve,
    evolveN,
    identityMatrix,
    matrixPower,
    chapmanKolmogorov,
) where

import Data.Semigroup (
    mtimesDefault,
 )
import Dtmc.Internal.Types (
    Distribution (..),
    TransitionMatrix (..),
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.Natural (Natural)
import Numeric.LinearAlgebra.Static qualified as S


evolve :: (KnownNat n) => Distribution n -> TransitionMatrix n -> Distribution n
evolve (Distribution v) (TransitionMatrix m) =
    Distribution
        { unDistribution = S.tr m S.#> v
        }

identityMatrix :: (KnownNat n) => TransitionMatrix n
identityMatrix = mempty

matrixPower :: (KnownNat n) => Natural -> TransitionMatrix n -> TransitionMatrix n
matrixPower = mtimesDefault

evolveN ::
    (KnownNat n) =>
    Natural ->
    Distribution n ->
    TransitionMatrix n ->
    Distribution n
evolveN k mu p =
    evolve mu (matrixPower k p)

chapmanKolmogorov ::
    (KnownNat n) =>
    Natural ->
    Natural ->
    TransitionMatrix n ->
    TransitionMatrix n
chapmanKolmogorov m n p =
    matrixPower m p <> matrixPower n p
