{- |
Module      : Dtmc.Distribution.Map
Description : Map-backed finite-support probability distributions.

t'DistributionMap' stores the nonzero coordinates of a probability law in a
'Data.Map.Strict.Map'. The represented state type is otherwise unrestricted.
-}
module Dtmc.Distribution.Map (
    DistributionMap,
    fromList,
    fromDistribution,
    pointMass,
    toMap,
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

{- | Construct a canonical finite-support probability law. Duplicate states
are combined, entries whose combined weight is exactly zero are removed, and
input order is ignored. Tolerated coordinate error is clamped to @[0, 1]@;
the repaired weights are normalised, and weights repaired to zero are omitted.

Complexity: @O(m log m)@ time for @m@ supplied entries, with @O(s)@ temporary
and result space for @s@ distinct states.
-}
fromList ::
    (Ord state) =>
    [(state, Double)] ->
    Either DistributionError (DistributionMap state)
fromList entries =
    DistributionMap
        . Map.fromDistinctAscList
        . filter ((/= 0) . snd)
        . zip (Map.keys combined)
        <$> first
            DistributionError
            (canonicaliseSimplexEntries (Map.elems combined))
  where
    combined = Map.filter (/= 0) (Map.fromListWith (+) entries)

{- | Construct the point mass concentrated on one state.

Complexity: @O(1)@ time and @O(1)@ result space.
-}
pointMass :: state -> DistributionMap state
pointMass state = DistributionMap (Map.singleton state 1)

{- | Convert any distribution representation to a map without revalidation or
renormalisation. Exact-zero coordinates are already omitted by the
'distributionWeights' contract.

Complexity: the cost of 'distributionWeights', plus @O(s)@ time and @O(s)@
temporary and result space for @s@ returned weights.
-}
fromDistribution ::
    (Distribution distribution) =>
    distribution ->
    DistributionMap (DistributionState distribution)
fromDistribution =
    DistributionMap . Map.fromDistinctAscList . distributionWeights

{- | Project the stored coordinates as a strict map. Exact-zero coordinates
are already omitted, so the result carries the mathematical support with its
weights.

Complexity: @O(1)@ time and space; the stored map is shared, not copied.
-}
toMap :: DistributionMap state -> Map.Map state Double
toMap = unDistributionMap
