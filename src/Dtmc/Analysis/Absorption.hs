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

{- | The transient and recurrent states, each in the canonical order of the
'FiniteState' instance. This is the ordering that puts the matrix in block
form, and the row and column ordering of 'fundamentalMatrix'.

Membership is decided from the support graph, so the split is exact and
involves no floating-point comparison.

Time: @O(n)@ after graph facts are cached.
-}
canonicalOrder ::
    (FiniteState state) =>
    TransitionMatrix state ->
    ([state], [state])
canonicalOrder p =
    (transientStates p, recurrentStates p)

{- | The fundamental matrix @G = (I - Q)^-1@ of the transient block, together
with the transient states that index its rows and columns.

Entry @(i,j)@ is @E(V_j | X_0 = i)@, the expected number of visits to
transient @j@ from transient @i@. Every entry is finite, so the result uses
'Double' rather than 'Expectation'.
'Dtmc.Analysis.VisitCount.totalExpectation' gives the same entries one at a
time and extends to recurrent targets, where the value is infinite.

A chain with no transient states returns @([], [])@ without a solve.
Otherwise the numerical behaviour and errors of the shared @(I - Q)@ solver
apply; @rho(Q) < 1@ holds in exact arithmetic because every transient state
reaches a recurrent one.

Time: @O(n^2 + t^3)@ for @t@ transient states. Result space: @O(t^2)@.
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

{- | Absorption probabilities into one recurrent state, in canonical state
order: coordinate @i@ is the probability that the first recurrent state
visited is the supplied target, started from @i@.

Boundary values are exact and taken without a solve:

* a target that is not recurrent gives an all-zero vector;
* a recurrent starting state has already arrived, so its coordinate is @1@
  when it is the target and @0@ otherwise.

Transient coordinates are the corresponding column of @B = G R'@ and inherit
the numerical behaviour and errors of 'fundamentalMatrix'.

Time: @O(n^2 + t^3)@. Result space: @O(n)@.
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

{- | The probability under an arbitrary initial distribution that the first
recurrent state visited is the supplied target. The result is the
initial-distribution mixture of the state-conditioned absorption
probabilities.
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

{- | Probability that the first recurrent state visited is the supplied
target, conditioned on @X_0 = i@.
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

{- | Expected number of transitions until absorption under an arbitrary
initial distribution.
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

{- | Expected number of transitions until absorption conditioned on
@X_0 = i@.
-}
expectationGivenInitialState ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError Expectation
expectationGivenInitialState p =
    Hitting.expectationGivenInitialState p (recurrentStates p)
