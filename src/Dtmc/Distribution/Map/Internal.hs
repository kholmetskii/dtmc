{- |
Module      : Dtmc.Distribution.Map.Internal
Description : Unsafe carrier for map-backed distributions.

The public smart constructor validates and canonicalises the simplex
invariant. Internal callers may construct values only when their operation
preserves that invariant up to floating-point error.
-}
module Dtmc.Distribution.Map.Internal (
    DistributionMap (DistributionMap),
    unDistributionMap,
    denseWeights,
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
    = -- | Wrap an unchecked state-to-weight map.
      DistributionMap (Map state Double)

type role DistributionMap nominal

deriving instance (Eq state) => Eq (DistributionMap state)
deriving instance (Show state) => Show (DistributionMap state)

{- | Return the canonical state-to-weight map without copying or validation.

Complexity: @O(1)@ time and @O(1)@ space.
-}
unDistributionMap :: DistributionMap state -> Map state Double
unDistributionMap (DistributionMap weights) = weights

{- | Return the weights of a map-backed distribution over a supplied ascending
state list, inserting exact zeros for absent states. Stored states absent from
the supplied list are ignored; lawful 'Dtmc.State.FiniteState' enumerations
contain every value of their state type.

Complexity: @O(n + s)@ time for @n@ requested states and stored support size
@s@, with @O(s)@ temporary space and @O(n)@ result space.
-}
denseWeights :: (Ord state) => [state] -> DistributionMap state -> [Double]
denseWeights states = align states . Map.toAscList . unDistributionMap
  where
    align [] _ = []
    align remaining [] = replicate (length remaining) 0
    align allStates@(state : rest) allWeights@((storedState, weight) : weights) =
        case compare storedState state of
            LT -> align allStates weights
            EQ -> weight : align rest weights
            GT -> 0 : align rest allWeights

instance Distribution (DistributionMap state) where
    type DistributionState (DistributionMap state) = state

    probabilityAt distribution state =
        Map.findWithDefault 0 state (unDistributionMap distribution)

    distributionWeights = Map.toAscList . unDistributionMap
