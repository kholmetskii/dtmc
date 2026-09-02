{- |
Module      : Dtmc.Distribution
Description : Shared abstraction for finite-support probability distributions.

The 'Distribution' class captures the read-only probability operations shared
by concrete distribution representations. Implementations live in
"Dtmc.Distribution.Vector" and "Dtmc.Distribution.Map".
-}
module Dtmc.Distribution (
    Distribution (..),
    DistributionError (..),
) where

import Dtmc.Simplex (
    SimplexError,
 )

{- | A simplex failure while constructing a distribution representation.
For a map-backed law, coordinate indices refer to ascending state order after
duplicate states have been combined and exact-zero weights omitted. The
wrapper keeps distribution failures distinct from transition-matrix row
failures.
-}
newtype DistributionError
    = -- | Wrap the underlying simplex failure.
      DistributionError SimplexError
    deriving (Eq, Show)

{- | A discrete probability distribution with finite stored support.

The class exposes observations common to every representation. Conversions
belong to the target representation module, so this abstraction does not
depend on a particular carrier.
-}
class Distribution distribution where
    -- | State type carried by the distribution representation.
    type DistributionState distribution

    {- | Read the stored probability of one state, returning exactly zero when
    the representation does not store that state. The value is returned
    without clamping or revalidation.
    -}
    probabilityAt ::
        (Ord (DistributionState distribution)) =>
        distribution ->
        DistributionState distribution ->
        Double

    {- | Return canonical ascending state weights. Exact-zero weights are
    omitted. Custom instances and numerically derived values are returned
    without revalidation.
    -}
    distributionWeights ::
        distribution ->
        [(DistributionState distribution, Double)]

    {- | States with strictly positive stored weight, in ascending order.
    Non-positive coordinates from custom instances or unchecked numerical
    operations are not mathematical support.
    -}
    support :: distribution -> [DistributionState distribution]
    support distribution =
        [ state
        | (state, weight) <- distributionWeights distribution
        , weight > 0
        ]
