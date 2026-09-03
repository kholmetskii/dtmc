{- |
Module      : Dtmc.Distribution.Map
Description : Map-backed finite-support probability distributions.

t'DistributionMap' stores the nonzero coordinates of a probability law in a
'Data.Map.Strict.Map'. The represented state type is otherwise unrestricted.
-}
module Dtmc.Distribution.Map (
    DistributionMap,
    mkDistributionMap,
    unDistributionMap,
    pointMass,
    toDistributionMap,
) where

import Data.Bifunctor (
    first,
 )
import Data.Map.Strict qualified as Map
import Dtmc.Distribution (
    Distribution (..),
    DistributionError (DistributionError),
 )
import Dtmc.Distribution.Map.Internal (
    DistributionMap (DistributionMap),
    unDistributionMap,
 )
import Dtmc.Simplex.Internal (
    canonicaliseSimplexEntries,
 )

{- | Combine duplicate states, remove entries whose combined weight is exactly
zero, and construct a canonical finite-support probability law. Tolerated
coordinate error is clamped to @[0, 1]@, the repaired weights are normalised,
and any weights repaired to zero are omitted. Input order is irrelevant.

Time: @O(m log m)@ for @m@ supplied entries. Space: @O(s)@ for @s@ distinct
states.
-}
mkDistributionMap ::
    (Ord state) =>
    [(state, Double)] ->
    Either DistributionError (DistributionMap state)
mkDistributionMap entries =
    DistributionMap
        . Map.fromDistinctAscList
        . filter ((/= 0) . snd)
        . zip (Map.keys combined)
        <$> first
            DistributionError
            (canonicaliseSimplexEntries (Map.elems combined))
  where
    combined = Map.filter (/= 0) (Map.fromListWith (+) entries)

-- | The point mass concentrated on one state. This is valid by construction.
pointMass :: state -> DistributionMap state
pointMass state = DistributionMap (Map.singleton state 1)

{- | Convert any distribution representation to a map without revalidation or
renormalisation. Exact-zero coordinates are already omitted by the
'distributionWeights' contract.

Time and result space: @O(s)@ for stored support size @s@.
-}
toDistributionMap ::
    (Distribution distribution) =>
    distribution ->
    DistributionMap (DistributionState distribution)
toDistributionMap =
    DistributionMap . Map.fromDistinctAscList . distributionWeights
