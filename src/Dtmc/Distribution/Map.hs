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
    mapStates,
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

{- | Push a distribution through a deterministic state mapping. Weights whose
states map to the same target are added, and an exact-zero combined weight is
removed.

No validation, clamping, or renormalisation is performed. A valid input
therefore remains a probability distribution up to ordinary floating-point
summation error.

Complexity: @O(s log(s + 1))@ time, @O(s)@ temporary space, and @O(r)@
result space for @s@ stored source states and @r@ distinct target states.
-}
mapStates ::
    (Ord target) =>
    (source -> target) ->
    DistributionMap source ->
    DistributionMap target
mapStates transform =
    DistributionMap
        . Map.filter (/= 0)
        . Map.mapKeysWith (+) transform
        . unDistributionMap

{- | Convert any distribution representation to a map without revalidation or
renormalisation. Weights reported for the same state are added and an
exact-zero combined weight is removed, so an instance that reports states out
of order, or reports one twice, still yields a structurally sound map. Whether
the reported weights form a probability law remains the obligation of the
'Distribution' instance.

Complexity: the cost of 'distributionWeights', plus @O(s log s)@ time and
@O(s)@ temporary and result space for @s@ returned weights.
-}
fromDistribution ::
    (Distribution distribution, Ord (DistributionState distribution)) =>
    distribution ->
    DistributionMap (DistributionState distribution)
fromDistribution =
    DistributionMap
        . Map.filter (/= 0)
        . Map.fromListWith (+)
        . distributionWeights

{- | Project the stored coordinates as a strict map. Exact-zero coordinates
are already omitted, so the result carries the mathematical support with its
weights.

Complexity: @O(1)@ time and space; the stored map is shared, not copied.
-}
toMap :: DistributionMap state -> Map.Map state Double
toMap = unDistributionMap
