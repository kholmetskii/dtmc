{- |
Module      : Dtmc.Distribution.Vector
Description : Dense probability vectors over finite state types.

t'DistributionVector' stores a probability law over a 'FiniteState' type in a
statically sized vector. Coordinates follow its canonical state order.
'fromList' checks and canonicalises the simplex invariant with the @1e-9@
tolerance documented by 'Dtmc.Simplex.SimplexError'. Explicit @hmatrix@
interoperability lives in "Dtmc.Distribution.Vector.HMatrix".
-}
module Dtmc.Distribution.Vector (
    DistributionVector,
    fromList,
    toList,
) where

import Dtmc.Distribution (
    DistributionError,
 )
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.Distribution.Map.Internal (
    denseWeights,
 )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector (DistributionVector),
    unDistributionVector,
 )
import Dtmc.State (
    FiniteState,
    finiteStates,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

{- | Construct a dense state distribution from labelled weights. Duplicate
states are combined, missing states receive weight zero, tolerated coordinate
error is clamped to @[0, 1]@, and the repaired weights are normalised before
storage. Input order is irrelevant. For a state type of cardinality zero,
returns @Left (DistributionError (SumOffBy 0))@.

Complexity: excluding evaluation of 'finiteStates', @O(m log m + n + s)@
time for @m@ supplied entries, @s@ distinct states, and state cardinality
@n@, with @O(n + s)@ temporary space and @O(n)@ result space.
-}
fromList ::
    (FiniteState state) =>
    [(state, Double)] ->
    Either DistributionError (DistributionVector state)
fromList entries = do
    distribution <- DistributionMap.fromList entries
    pure $
        DistributionVector $
            S.vector (denseWeights finiteStates distribution)

{- | Return every stored coordinate in canonical state order, including exact
zeros. This is a representation-neutral copy of the dense vector.

Complexity: @O(n)@ time and @O(n)@ temporary and result space for state
cardinality @n@.
-}
toList :: (FiniteState state) => DistributionVector state -> [Double]
toList = LA.toList . S.extract . unDistributionVector
