{- |
Module      : Dtmc.Dynamics
Description : Deterministic forward evolution of distributions.

Deterministic push-forward of a state distribution through a DTMC. For a row
transition matrix @P@ and column distribution @mu@,
@mu'(j) = sum_i mu(i) P(i,j)@.
-}
module Dtmc.Dynamics (
    evolve,
    evolveN,
    probabilityAtTime,
) where

import Data.Finite (
    Finite,
 )
import Dtmc.Distribution (
    probabilityAt,
 )
import Dtmc.Distribution.Internal (
    Distribution (Distribution),
 )
import Dtmc.TransitionMatrix (
    matrixPower,
 )
import Dtmc.TransitionMatrix.Internal (
    TransitionMatrix,
    unTransitionMatrix,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (Natural)

{- | The next-state distribution @mu' = transpose(P) mu@.

Exact probability inputs produce a probability distribution. The result is
wrapped without validation, clamping, or renormalisation, so tolerated input
error and floating-point rounding are preserved and may make a subsequent
validation fail.

Time: @O(n^2)@. Result space: @O(n)@.
-}
evolve :: (KnownNat n) => Distribution n -> TransitionMatrix n -> Distribution n
evolve (Distribution v) p =
    Distribution (S.tr (unTransitionMatrix p) S.#> v)

{- | The distribution after @k@ transitions, computed as
@evolve mu (matrixPower k p)@. Exponent zero is the original distribution
mathematically.

This powers the matrix rather than iterating 'evolve', so the two calculations
may differ by floating-point rounding. The result is not revalidated.

Time: @O(n^2 + n^3 log(k + 1))@.
-}
evolveN ::
    (KnownNat n) =>
    Natural ->
    Distribution n ->
    TransitionMatrix n ->
    Distribution n
evolveN k mu p =
    evolve mu (matrixPower k p)

{- | The marginal probability @P(X_k = j)@ that the chain occupies state @j@
after @k@ steps, started from the initial law @initial@ and driven by @p@. It
reads the target coordinate of the evolved distribution:

@
probabilityAtTime k initial p j == 'probabilityAt' ('evolveN' k initial p) j
@

At @k = 0@ no step has been taken, so the result is
@'probabilityAt' initial j@, the initial probability of @j@.

This is distinct from 'probabilityAt': 'probabilityAt' reads a coordinate from
an already-computed distribution, whereas 'probabilityAtTime' first evolves
@initial@ through @k@ steps of @p@ and then reads coordinate @j@. The 'Double'
result is not clamped, renormalised, or revalidated.

Time: @O(n^2 + n^3 log(k + 1))@, dominated by the matrix power inside
'evolveN'.
-}
probabilityAtTime ::
    (KnownNat n) =>
    Natural ->
    Distribution n ->
    TransitionMatrix n ->
    Finite n ->
    Double
probabilityAtTime k initial p = probabilityAt (evolveN k initial p) 
