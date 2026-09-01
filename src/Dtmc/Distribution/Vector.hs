{- |
Module      : Dtmc.Distribution.Vector
Description : Dense probability vectors over finite state types.

t'DistributionVector' stores a probability law over a 'FiniteState' type in a
statically sized vector. Coordinates follow its canonical state order.
'mkDistributionVector' checks the simplex invariant with the @1e-9@ tolerance
documented by 'Dtmc.Simplex.SimplexError'.
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
    validateSimplex,
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
 )
import Numeric.LinearAlgebra.Static qualified as S

{- | Validate a raw state distribution vector. On success it is preserved
exactly; tolerated coordinates are not clamped and the total is not
renormalised. For a state type of cardinality zero, returns
@Left (DistributionError (SumOffBy 0))@.

Time: @O(k)@. Space: @O(k)@ for validation, where @k@ is the state
cardinality.
-}
mkDistributionVector ::
    (FiniteState state) =>
    S.R (Cardinality state) ->
    Either DistributionError (DistributionVector state)
mkDistributionVector vector =
    DistributionVector vector <$ first DistributionError (validateSimplex vector)
