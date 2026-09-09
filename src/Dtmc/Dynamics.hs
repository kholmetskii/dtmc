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
import Numeric.LinearAlgebra qualified as LA
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
    DistributionVector state ->
    TransitionMatrix state ->
    DistributionVector state
evolveVector (DistributionVector v) p =
    DistributionVector (LA.tr (unTransitionMatrix p) LA.#> v)

{- | Compute the distribution after @k@ transitions. A cost model chooses
between repeated matrix-vector multiplication and powering the matrix, so a
moderate number of steps does not construct a full matrix power unnecessarily.
Exponent zero returns the original distribution.

The two strategies are mathematically equivalent but may differ by ordinary
floating-point rounding. The result is not revalidated.

Complexity: @O(k n^2)@ time when @k <= n@ and
@O(n^2 + n^3 log(k + 1))@ otherwise. Temporary space is @O(n)@ in the
iterated case and @O(n^2)@ in the powered case; result space is @O(n)@.
-}
evolveVectorN ::
    (FiniteState state) =>
    Natural ->
    DistributionVector state ->
    TransitionMatrix state ->
    DistributionVector state
evolveVectorN k mu@(DistributionVector initial) p
    | k == 0 = mu
    | useIteration = DistributionVector (iterateVector k initial)
    | otherwise = evolveVector mu (power k p)
  where
    matrix = unTransitionMatrix p
    transposed = LA.tr matrix
    dimension = LA.rows matrix

    -- A deliberately conservative threshold retains the highly tuned matrix
    -- power path for long runs on small matrices. Focused benchmarks show the
    -- repeated matrix-vector path winning once the matrix dimension reaches
    -- the requested step count.
    useIteration = toInteger k <= toInteger dimension

    iterateVector 0 vector = vector
    iterateVector remaining vector =
        let next = transposed LA.#> vector
         in next `seq` iterateVector (remaining - 1) next

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
