{- |
Module      : Dtmc.Distribution.Vector
Description : Dense probability vectors over finite state types.

t'DistributionVector' stores a probability law over a 'FiniteState' type in a
statically sized vector. Coordinates follow its canonical state order.
'mkDistributionVector' checks and canonicalises the simplex invariant with the
@1e-9@ tolerance documented by 'Dtmc.Simplex.SimplexError'.
-}
module Dtmc.Distribution.Vector (
    DistributionVector,
    mkDistributionVector,
    unDistributionVector,
) where

import Data.Bifunctor (
    first,
 )
import Dtmc.Distribution (
    DistributionError (DistributionError),
 )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector (DistributionVector),
    unDistributionVector,
 )
import Dtmc.Simplex.Internal (
    canonicaliseSimplex,
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
 )
import Numeric.LinearAlgebra.Static qualified as S

{- | Construct a state distribution vector. Tolerated coordinate error is
clamped to @[0, 1]@ and the repaired vector is normalised before storage. For
a state type of cardinality zero, returns
@Left (DistributionError (SumOffBy 0))@.

Time and result space: @O(k)@, where @k@ is the state
cardinality.
-}
mkDistributionVector ::
    (FiniteState state) =>
    S.R (Cardinality state) ->
    Either DistributionError (DistributionVector state)
mkDistributionVector vector = do
    canonical <- first DistributionError (canonicaliseSimplex vector)
    pure (DistributionVector canonical)
