{- |
Module      : Dtmc.Dynamics
Description : Deterministic forward evolution of distributions.

Deterministic push-forward of a state distribution through a DTMC. Dense
finite laws use transition matrices; sparse finite-support laws use any
locally finite 'Transition'. In both cases,
@mu'(j) = sum_i mu(i) P(i,j)@.
-}
module Dtmc.Dynamics (
    evolve,
    evolveN,
    evolveVector,
    evolveVectorN,
) where

import Dtmc.Distribution (
    Distribution (..),
 )
import Dtmc.Distribution.Map (
    toDistributionMap,
 )
import Dtmc.Distribution.Map.Internal (
    DistributionMap (DistributionMap),
    unDistributionMap,
 )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector (DistributionVector),
 )
import Dtmc.Dynamics.Internal (
    pushSparseWeights,
 )
import Dtmc.Transition (
    Transition (..),
 )
import Dtmc.Transition.Matrix (
    matrixPower,
 )
import Dtmc.Transition.Matrix.Internal (
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
evolveVector ::
    (KnownNat n) =>
    DistributionVector n ->
    TransitionMatrix n ->
    DistributionVector n
evolveVector (DistributionVector v) p =
    DistributionVector (S.tr (unTransitionMatrix p) S.#> v)

{- | The distribution after @k@ transitions, computed as
@evolveVector mu (matrixPower k p)@. Exponent zero is the original distribution
mathematically.

This powers the matrix rather than iterating 'evolveVector', so the two calculations
may differ by floating-point rounding. The result is not revalidated.

Time: @O(n^2 + n^3 log(k + 1))@.
-}
evolveVectorN ::
    (KnownNat n) =>
    Natural ->
    DistributionVector n ->
    TransitionMatrix n ->
    DistributionVector n
evolveVectorN k mu p =
    evolveVector mu (matrixPower k p)

{- | Push any finite-support 'Distribution' through one locally finite kernel
step. The result uses t'DistributionMap' because a general kernel does not
provide a finite global state enumeration. It is not revalidated, clamped, or
renormalised.

Time is @O(e log r)@ for @e@ traversed support edges and @r@ result states.
-}
evolve ::
    ( Distribution distribution
    , Transition kernel
    , DistributionState distribution ~ TransitionState kernel
    , Ord (TransitionState kernel)
    ) =>
    distribution ->
    kernel ->
    DistributionMap (TransitionState kernel)
evolve distribution kernel =
    DistributionMap
        ( pushSparseWeights
            (unDistributionMap (toDistributionMap distribution))
            kernel
        )

{- | Apply 'evolve' exactly @k@ times. At @k = 0@ the initial law is converted
to an equivalent t'DistributionMap' without revalidation. No state-space
enumeration or truncation is performed.
-}
evolveN ::
    ( Distribution distribution
    , Transition kernel
    , DistributionState distribution ~ TransitionState kernel
    , Ord (TransitionState kernel)
    ) =>
    Natural ->
    distribution ->
    kernel ->
    DistributionMap (TransitionState kernel)
evolveN steps initial kernel = go steps (toDistributionMap initial)
  where
    go 0 distribution = distribution
    go remaining distribution =
        let next = evolve distribution kernel
         in next `seq` go (remaining - 1) next
