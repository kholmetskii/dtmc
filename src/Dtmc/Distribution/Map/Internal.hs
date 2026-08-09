{- |
Module      : Dtmc.Distribution.Map.Internal
Description : Unsafe carrier for map-backed distributions.

The public smart constructor validates the simplex invariant. Internal callers
may construct values only when their operation preserves that invariant up to
floating-point error.
-}
module Dtmc.Distribution.Map.Internal (
    DistributionMap (DistributionMap),
    unDistributionMap,
) where

import Data.Map.Strict (
    Map,
 )
import Data.Map.Strict qualified as Map
import Dtmc.Distribution (
    Distribution (..),
 )

{- | A finite-support probability distribution backed by a strict map. The
internal constructor performs no validation.
-}
newtype DistributionMap state
    = DistributionMap (Map state Double)

type role DistributionMap nominal

deriving instance (Eq state) => Eq (DistributionMap state)
deriving instance (Show state) => Show (DistributionMap state)

-- | Return the canonical state-to-weight map without copying or validation.
unDistributionMap :: DistributionMap state -> Map state Double
unDistributionMap (DistributionMap weights) = weights

instance Distribution (DistributionMap state) where
    type DistributionState (DistributionMap state) = state

    probabilityAt distribution state =
        Map.findWithDefault 0 state (unDistributionMap distribution)

    distributionWeights = Map.toAscList . unDistributionMap
