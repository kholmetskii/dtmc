module Dtmc.Distribution (
    Distribution,
    DistributionError (..),
    mkDistribution,
    unDistribution,
    approxDistributionEq,
) where

import Data.Bifunctor (
    first,
 )
import Dtmc.Internal.Simplex (
    SimplexError,
    validateSimplex,
 )
import Dtmc.Internal.Types (
    Distribution (..),
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

newtype DistributionError = DistributionError SimplexError
    deriving (Eq, Show)

mkDistribution :: (KnownNat n) => S.R n -> Either DistributionError (Distribution n)
mkDistribution vector =
    Distribution vector <$ first DistributionError (validateSimplex vector)

approxDistributionEq :: (KnownNat n) => Double -> Distribution n -> Distribution n -> Bool
approxDistributionEq tolerance (Distribution left) (Distribution right) =
    and (zipWith close (entries left) (entries right))
  where
    entries = LA.toList . S.extract
    close x y = abs (x - y) <= tolerance
