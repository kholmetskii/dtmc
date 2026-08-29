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
    expectation,
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
import Numeric.LinearAlgebra qualified as LA
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

{- | The expectation of a state-indexed quantity under this distribution,
@E_mu(f) = sum_i mu(i) f(i)@.

The second argument holds one value per state in the canonical order of the
'FiniteState' instance -- exactly the shape returned by the @...ByState@
analysis queries. Applying this to such a result turns a query conditioned on
one initial state into the same query under an arbitrary initial law, because
probability is linear in the initial distribution:

@P(A) = sum_j mu(j) P(A | X_0 = j)@

For example, with @v = eventualProbabilityByState p targets@ the value
@expectation mu v@ is @P_mu(H_A < infinity)@.

Ordinary floating-point summation. The result is not clamped, so a quantity
that lies outside @[0, 1]@ stays outside it. The result is a plain 'Double';
the @Expectation@ type used by hitting- and return-time means is unrelated and
carries their separate infinite case.

Time: @O(k)@ for state cardinality @k@. Space: @O(k)@ for extraction.
-}
expectation ::
    (FiniteState state) =>
    DistributionVector state ->
    S.R (Cardinality state) ->
    Double
expectation (DistributionVector weights) values =
    S.extract weights LA.<.> S.extract values
