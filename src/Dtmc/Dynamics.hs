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
) where

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
