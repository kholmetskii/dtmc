{- |
Module      : Dtmc.Analysis.Absorption
Description : Canonical decomposition, fundamental matrix, and absorption.

Absorption analysis for a finite chain. Ordering the states so that the
transient set @T@ comes first and the recurrent set @R@ last puts the
transition matrix in block form

@
P = [ Q  R' ]
    [ 0  S  ]
@

and the fundamental matrix of the transient block is
@G = sum_(n >= 0) Q^n = (I - Q)^-1@, whose entry @G(i,j)@ is the expected
number of visits to @j@ starting from @i@ for transient @i@ and @j@.

Every finite chain has this decomposition, including the degenerate cases
@T = empty@ (no transient states) and a single recurrent class. No witness
type is required: unlike stationarity, the analysis is defined for every
finite transition matrix.

Absorption probabilities are @B = G R'@, where @B(i,k)@ is the probability
that the /first/ recurrent state the chain visits is @k@. This is not the
same as the probability of ever visiting @k@: a recurrent class with more
than one state can be entered at one member and later reach another, which
'Dtmc.Analysis.HittingTime.eventualProbability' counts and @B@ does not. The
two agree only after summing over a whole recurrent class.

Unless stated otherwise, complexity bounds exclude 'FiniteState' method
costs. For those bounds, @n@ is the state count, @E@ the support-edge count,
@t@ the transient-state count, and @s@ an initial distribution's stored
support size.
-}
module Dtmc.Analysis.Absorption (
    LinearSystemError (..),
    Expectation (..),
    canonicalOrder,
    fundamentalMatrix,
    probability,
    probabilityGivenInitialState,
    expectation,
    expectationGivenInitialState,
) where

import Data.Array.Unboxed qualified as Unboxed
import Dtmc.Analysis.Classification (
    recurrentState,
    recurrentStates,
    transientStates,
 )
import Dtmc.Analysis.Expectation (
    Expectation (..),
 )
import Dtmc.Analysis.HittingTime qualified as Hitting
import Dtmc.Analysis.Initial.Internal (
    expectationUnderEither,
    probabilityUnderEither,
 )
import Dtmc.Analysis.LinearSystem (
    LinearSystemError (..),
 )
import Dtmc.Analysis.LinearSystem.Internal (
    fundamental,
    subMatrix,
 )
import Dtmc.Distribution (
    Distribution (..),
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
 )
import Dtmc.State.Internal (
    stateCardinalityInt,
    stateIndexInt,
 )
import Dtmc.Transition.Matrix.Internal (
    TransitionMatrix,
    unTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

toIndex :: (FiniteState state) => state -> Int
toIndex = stateIndexInt

{- | Return the transient and recurrent states, each in the canonical order
of the 'FiniteState' instance. This ordering puts the matrix in block form
and indexes the rows and columns of 'fundamentalMatrix'.

Membership is decided from the support graph, so the split is exact and
involves no floating-point comparison.

Complexity: with shared graph facts cached, @O(n)@ time and @O(n)@ temporary
and result space. On an unforced matrix, the first full evaluation takes
@O(n^2 + (n + E) log(n + 1))@ time, @O(n^2 + n + E)@ temporary space, and
retains @O(n + E)@ graph-cache space.
-}
canonicalOrder ::
    (FiniteState state) =>
    TransitionMatrix state ->
    ([state], [state])
canonicalOrder p =
    (transientStates p, recurrentStates p)

{- | Compute the fundamental matrix @G = (I - Q)^-1@ of the transient block,
together with the transient states that index its rows and columns.

Entry @(i,j)@ is @E(V_j | X_0 = i)@, the expected number of visits to
transient @j@ from transient @i@. Every entry is finite, so the result uses
'Double' rather than 'Expectation'.
'Dtmc.Analysis.VisitCount.totalExpectation' gives the same entries one at a
time and extends to recurrent targets, where the value is infinite.

A chain with no transient states returns @([], [])@ without a solve.
Otherwise the numerical behaviour and errors of the shared @(I - Q)@ solver
apply; @rho(Q) < 1@ holds in exact arithmetic because every transient state
reaches a recurrent one.

Complexity: including first-time graph classification,
@O(n^2 + (n + E) log(n + 1) + t^3)@ time,
@O(n^2 + n + E + t^2)@ temporary space, @O(n + E)@ retained graph-cache
space, and @O(t^2)@ result space.
-}
fundamentalMatrix ::
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError ([state], [[Double]])
fundamentalMatrix p
    | null transient = Right ([], [])
    | otherwise = do
        g <- fundamental (subMatrix transientIdx transientIdx matrix)
        pure (transient, LA.toLists g)
  where
    transient = transientStates p
    transientIdx = map toIndex transient
    matrix = S.extract (unTransitionMatrix p)

{- | Compute absorption probabilities into one recurrent state in canonical
state order. Coordinate @i@ is the probability that the supplied target is
the first recurrent state visited when starting from @i@.

Boundary values are exact and taken without a solve:

* a target that is not recurrent gives an all-zero vector;
* a recurrent starting state has already arrived, so its coordinate is @1@
  when it is the target and @0@ otherwise.

Transient coordinates are the corresponding column of @B = G R'@ and inherit
the numerical behaviour and errors of 'fundamentalMatrix'.

Complexity: including first-time graph classification,
@O(n^2 + (n + E) log(n + 1) + t^3)@ worst-case time,
@O(n^2 + n + E + t^2)@ temporary space, @O(n + E)@ retained graph-cache
space, and @O(n)@ result space. A non-recurrent target avoids the numerical
solve.
-}
probabilityByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError (S.R (Cardinality state))
probabilityByState p target
    | not (recurrentState p target) =
        Right (S.vector (replicate dim 0))
    | null transientIdx =
        Right (S.vector [arrived i | i <- [0 .. dim - 1]])
    | otherwise = do
        g <- fundamental (subMatrix transientIdx transientIdx matrix)
        let exits = LA.flatten (subMatrix transientIdx [targetIdx] matrix)
            solved = LA.toList (g LA.#> exits)
            interior :: Unboxed.UArray Int Double
            interior =
                Unboxed.accumArray
                    (\_ x -> x)
                    0
                    (0, dim - 1)
                    (zip transientIdx solved)
            valueAt i
                | transientMask Unboxed.! i = interior Unboxed.! i
                | otherwise = arrived i
        pure (S.vector [valueAt i | i <- [0 .. dim - 1]])
  where
    dim = stateCardinalityInt @state
    matrix = S.extract (unTransitionMatrix p)
    targetIdx = toIndex target
    transientIdx = map toIndex (transientStates p)
    transientMask :: Unboxed.UArray Int Bool
    transientMask =
        Unboxed.accumArray
            (\_ x -> x)
            False
            (0, dim - 1)
            [(i, True) | i <- transientIdx]
    arrived i = if i == targetIdx then 1 else 0

{- | Compute, under an arbitrary initial distribution, the probability that
the supplied target is the first recurrent state visited. The result is the
initial-law mixture of the state-conditioned absorption probabilities. A
non-recurrent target gives exactly zero without a numerical solve.

Complexity: excluding 'distributionWeights', @O(n^3 + s)@ worst-case time,
@O(n^2 + s)@ temporary space, and @O(1)@ result space. The matrix may retain
@O(n + E)@ graph-cache space.
-}
probability ::
    ( FiniteState state
    , Distribution distribution
    , DistributionState distribution ~ state
    ) =>
    TransitionMatrix state ->
    state ->
    distribution ->
    Either LinearSystemError Double
probability p target initial =
    probabilityUnderEither initial (probabilityGivenInitialState p target)

{- | Compute the probability that the supplied target is the first recurrent
state visited, conditioned on @X_0 = i@. A recurrent initial state is already
absorbed at time zero.

Partial application shares one lazy all-state table. A non-recurrent target
produces exact zeros without a numerical solve.

Complexity: the first forced query takes @O(n^3)@ worst-case time and
@O(n^2)@ temporary space and may retain an @O(n)@ all-state result and
@O(n + E)@ graph cache. Subsequent shared lookups take @O(1)@ time and
space; the scalar result occupies @O(1)@ space.
-}
probabilityGivenInitialState ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    state ->
    Either LinearSystemError Double
probabilityGivenInitialState p target =
    \i -> (`LA.atIndex` toIndex i) <$> values
  where
    values = S.extract <$> probabilityByState p target

{- | Compute the expected number of transitions until the chain first enters
the recurrent states under an arbitrary initial distribution. The value is
mathematically finite for every valid finite chain; numerical failures come
from the checked transient-state solve.

Complexity: excluding 'distributionWeights', @O(n^3 + s)@ worst-case time,
@O(n^2 + s)@ temporary space, and @O(1)@ result space. The matrix may retain
@O(n + E)@ graph-cache space.
-}
expectation ::
    ( FiniteState state
    , Distribution distribution
    , DistributionState distribution ~ state
    ) =>
    TransitionMatrix state ->
    distribution ->
    Either LinearSystemError Expectation
expectation p initial =
    expectationUnderEither initial (expectationGivenInitialState p)

{- | Compute the expected number of transitions until the chain first enters
the recurrent states, conditioned on @X_0 = i@. A recurrent initial state has
expectation zero; every transient state has a mathematically finite value.

Partial application shares the lazy all-state result of the checked
transient-state solve.

Complexity: the first forced query takes @O(n^3)@ worst-case time and
@O(n^2)@ temporary space and may retain an @O(n)@ all-state result and
@O(n + E)@ graph cache. Subsequent shared lookups take @O(1)@ time and
space; the scalar result occupies @O(1)@ space.
-}
expectationGivenInitialState ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError Expectation
expectationGivenInitialState p =
    Hitting.expectationGivenInitialState p (recurrentStates p)
