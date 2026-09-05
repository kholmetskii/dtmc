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
    fromDistribution,
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
import Dtmc.State (
    FiniteState,
 )
import Dtmc.Transition (
    Transition (..),
 )
import Dtmc.Transition.Matrix (
    power,
 )
import Dtmc.Transition.Matrix.Internal (
    TransitionMatrix,
    unTransitionMatrix,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (Natural)

{- | Compute the next-state distribution @mu' = transpose(P) mu@.

Exact probability inputs produce a probability distribution. The result is
wrapped without validation, clamping, or renormalisation, so error from custom
or numerically derived inputs and floating-point rounding is preserved and may
make a subsequent validation fail.

Complexity: @O(n^2)@ time, @O(n^2)@ temporary space in the worst case, and
@O(n)@ result space for state cardinality @n@.
-}
evolveVector ::
    (FiniteState state) =>
    DistributionVector state ->
    TransitionMatrix state ->
    DistributionVector state
evolveVector (DistributionVector v) p =
    DistributionVector (S.tr (unTransitionMatrix p) S.#> v)

{- | Compute the distribution after @k@ transitions as
@evolveVector mu (power k p)@. Exponent zero is the original distribution
mathematically.

This powers the matrix rather than iterating 'evolveVector', so the two
calculations may differ by floating-point rounding. The result is not
revalidated.

Complexity: @O(n^2 + n^3 log(k + 1))@ time, @O(n^2)@ temporary space, and
@O(n)@ result space.
-}
evolveVectorN ::
    (FiniteState state) =>
    Natural ->
    DistributionVector state ->
    TransitionMatrix state ->
    DistributionVector state
evolveVectorN k mu p =
    evolveVector mu (power k p)

{- | Push any finite-support 'Distribution' through one locally finite kernel
step. The result uses t'DistributionMap' because a general kernel does not
provide a finite global state enumeration. It is not revalidated, clamped, or
renormalised.

For the complexity bounds, @s@ is the number of source states, @e@ the number
of traversed support edges, @u@ the number of distinct destinations
encountered, and @r@ the number retained after exact-zero removal.

Complexity: excluding 'distributionWeights' and 'transitionLaw' evaluation,
@O(s + e log(u + 1) + u)@ time, @O(s + u)@ temporary space, and @O(r)@ result
space.
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
            (unDistributionMap (fromDistribution distribution))
            kernel
        )

{- | Apply 'evolve' exactly @k@ times. At @k = 0@ the initial law is converted
to an equivalent t'DistributionMap' without revalidation. No state-space
enumeration or truncation is performed.

For a positive step count, let @s@, @e@, and @u@ be upper bounds per step on
the source states, traversed support edges, and distinct destinations
encountered; let @r@ be the final support size.

Complexity: excluding the initial 'distributionWeights' call and all
'transitionLaw' evaluations, @O(k (s + e log(u + 1) + u))@ time,
@O(s + u)@ temporary space, and @O(r)@ result space. At @k = 0@, the cost is
that of 'Dtmc.Distribution.Map.fromDistribution'.
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
evolveN steps initial kernel = go steps (fromDistribution initial)
  where
    go 0 distribution = distribution
    go remaining distribution =
        let next = evolve distribution kernel
         in next `seq` go (remaining - 1) next
