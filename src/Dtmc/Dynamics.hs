{- |
Module      : Dtmc.Dynamics
Description : Deterministic forward evolution of distributions.

Deterministic push-forward of a state distribution through a DTMC. Dense
finite laws use transition matrices; sparse finite-support laws use any
locally finite 'MarkovKernel'. In both cases,
@mu'(j) = sum_i mu(i) P(i,j)@.
-}
module Dtmc.Dynamics (
    evolve,
    evolveN,
    evolveSparse,
    evolveSparseN,
) where

import Dtmc.Distribution.Internal (
    Distribution (Distribution),
    SparseDistribution (SparseDistribution),
    unSparseDistribution,
 )
import Dtmc.Dynamics.Internal (
    pushSparseWeights,
 )
import Dtmc.Kernel (
    MarkovKernel (..),
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

{- | Push a sparse distribution through one locally finite kernel step. The
result is sparse and is not revalidated, clamped, or renormalised.

Time is @O(e log r)@ for @e@ traversed support edges and @r@ result states.
-}
evolveSparse ::
    (MarkovKernel kernel, Ord (KernelState kernel)) =>
    SparseDistribution (KernelState kernel) ->
    kernel ->
    SparseDistribution (KernelState kernel)
evolveSparse distribution kernel =
    SparseDistribution
        (pushSparseWeights (unSparseDistribution distribution) kernel)

{- | Apply 'evolveSparse' exactly @k@ times. At @k = 0@ the original sparse
distribution is returned unchanged. No state-space enumeration or truncation
is performed.
-}
evolveSparseN ::
    (MarkovKernel kernel, Ord (KernelState kernel)) =>
    Natural ->
    SparseDistribution (KernelState kernel) ->
    kernel ->
    SparseDistribution (KernelState kernel)
evolveSparseN steps initial kernel = go steps initial
  where
    go 0 distribution = distribution
    go remaining distribution =
        let next = evolveSparse distribution kernel
         in next `seq` go (remaining - 1) next
