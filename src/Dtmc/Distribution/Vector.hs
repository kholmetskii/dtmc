{- |
Module      : Dtmc.Distribution.Vector
Description : Dense probability vectors over finite state types.

t'DistributionVector' stores a probability law over a 'FiniteState' type in a
statically sized vector. Coordinates follow its canonical state order, so
'fromList' and 'toList' are a positional pair: both speak the same list of
weights, one coordinate per state. 'fromList' checks and canonicalises the
simplex invariant with the @1e-9@ tolerance documented by
'Dtmc.Simplex.SimplexError'.

To build a vector from /labelled/ weights, where duplicates should combine
and missing states should default to zero, use
'Dtmc.Distribution.Map.fromList' and read the coordinates off the result:

> Vector.fromList [probabilityAt m s | s <- finiteStates]
-}
module Dtmc.Distribution.Vector (
    DistributionVector,
    DistributionVectorError (..),
    fromList,
    toList,
) where

import Data.Bifunctor (
    bimap,
 )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector (DistributionVector),
    unDistributionVector,
 )
import Dtmc.Simplex (
    SimplexError,
 )
import Dtmc.Simplex.Internal (
    canonicaliseSimplexEntries,
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.State.Internal (
    stateCardinalityInt,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

{- | Why a list of weights was rejected as a state distribution.
-}
data DistributionVectorError
    = -- | The state cardinality and the supplied number of weights.
      WrongLength Int Int
    | -- | The weights failed simplex validation. The coordinate index the
      -- 'SimplexError' carries is zero-based in canonical state order.
      InWeights SimplexError
    deriving (Eq, Show)

{- | Construct a dense state distribution from one weight per state, in
canonical state order. The list must have exactly as many entries as the
state type has inhabitants; tolerated coordinate error is clamped to @[0, 1]@
and the repaired weights are normalised before storage. This is the exact
inverse of 'toList' up to that repair.

For a state type of cardinality zero the only accepted input is @[]@, which
is rejected as @Left (InWeights (SumOffBy 0))@: the empty simplex has no
points.

Complexity: @O(n)@ time and @O(n)@ temporary and result space for state
cardinality @n@.
-}
fromList ::
    forall state.
    (FiniteState state) =>
    [Double] ->
    Either DistributionVectorError (DistributionVector state)
fromList weights
    | supplied /= dimension = Left (WrongLength dimension supplied)
    | otherwise =
        bimap InWeights (DistributionVector . S.vector) canonicalised
  where
    dimension = stateCardinalityInt @state
    supplied = length weights
    canonicalised = canonicaliseSimplexEntries weights

{- | Return every stored coordinate in canonical state order, including exact
zeros. This is a representation-neutral copy of the dense vector.

Complexity: @O(n)@ time and @O(n)@ temporary and result space for state
cardinality @n@.
-}
toList :: (FiniteState state) => DistributionVector state -> [Double]
toList = LA.toList . S.extract . unDistributionVector
