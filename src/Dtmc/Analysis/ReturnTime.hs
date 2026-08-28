{- |
Module      : Dtmc.Analysis.ReturnTime
Description : Exact, bounded, eventual, and expected first-return times.

First-return quantities for DTMCs. For state @i@,
@T_i = inf { t >= 1 | X_t = i }@, so time zero is never a return. Scalar
exact-time and bounded queries work through any locally finite 'Transition'.
All-state, eventual, and expected queries use a finite 'TransitionMatrix'.

Exact-time and strictly bounded queries use finite recurrences. Eventual and
expected queries use support classification and checked 'Double' linear
solves. Results are not clamped or renormalised.
-}
module Dtmc.Analysis.ReturnTime (
    LinearSystemError (..),
    Expectation (..),
    probability,
    probabilityByState,
    eventualProbability,
    eventualProbabilityByState,
    expectation,
    expectationByState,
) where

import Control.Monad (
    foldM,
 )
import Data.Array.Unboxed qualified as Unboxed
import Data.Finite (
    getFinite,
 )
import Data.Map.Strict (
    Map,
 )
import Data.Map.Strict qualified as Map
import Data.Proxy (
    Proxy (..),
 )
import Dtmc.Analysis.Classification (
    recurrentState,
    transientStates,
 )
import Dtmc.Analysis.Expectation (
    Expectation (..),
 )
import Dtmc.Analysis.Event (
    DiscreteEvent (..),
 )
import Dtmc.Analysis.HittingTime qualified as Hit
import Dtmc.Analysis.Internal.LinearSystem (
    fundamental,
    subMatrix,
 )
import Dtmc.Analysis.LinearSystem (
    LinearSystemError (..),
 )
import Dtmc.Distribution.Vector.Internal (
    unDistributionVector,
 )
import Dtmc.Dynamics.Internal (
    pushSparseWeights,
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
    finiteStates,
    stateIndex,
 )
import Dtmc.Transition (
    Transition (..),
 )
import Dtmc.Transition.Matrix (
    rowAt,
 )
import Dtmc.Transition.Matrix.Internal (
    TransitionMatrix,
    unTransitionMatrix,
 )
import GHC.TypeNats (
    natVal,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )

toIndex :: (FiniteState state) => state -> Int
toIndex = fromIntegral . getFinite . stateIndex

iterateNatural :: Natural -> (a -> a) -> a -> a
iterateNatural steps f = go steps
  where
    go 0 value = value
    go remaining value =
        let next = f value
         in next `seq` go (remaining - 1) next

advanceUntilTarget ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    kernel ->
    (TransitionState kernel -> Bool) ->
    Map (TransitionState kernel) Double ->
    (Map (TransitionState kernel) Double, Double)
advanceUntilTarget kernel isTarget survivors =
    (remaining, hitMass)
  where
    advanced = pushSparseWeights survivors kernel
    (hits, remaining) = Map.partitionWithKey (\state _ -> isTarget state) advanced
    hitMass = sum (Map.elems hits)

-- that origin. Row @i@ holds the current-state mass of paths started at @i@
-- that avoided @i@ at every positive time so far. The diagonal of the
-- advanced matrix is therefore the first-return mass at the new time, and
-- clearing it removes those completed paths from later steps.
firstReturnStep ::
    LA.Matrix Double ->
    LA.Matrix Double ->
    (LA.Matrix Double, LA.Vector Double)
firstReturnStep matrix survivors =
    (advanced - LA.diag masses, masses)
  where
    advanced = survivors LA.<> matrix
    masses = LA.takeDiag advanced

{- | Exact first-return probabilities
@f_i(t) = P(T_i = t | X_0 = i)@ in state order, where
@T_i = inf { r >= 1 | X_r = i }@.

At @t = 0@ every entry is exactly @0@. At @t = 1@ the result is the diagonal
of the transition matrix, so an absorbing state has mass exactly @1@ at time
one and exactly @0@ at every later first-return time.

The implementation evolves an @n x n@ survivor matrix whose row @i@ contains
the mass of paths started at @i@ that have not returned to @i@. At each step
the new diagonal is the first-return mass and is cleared before continuing.
This computes every starting state's law together without a linear solve.
Results use ordinary 'Double' arithmetic without clamping or renormalisation.

Time: @O(t n^3)@; temporary space: @O(n^2)@; result space: @O(n)@.
-}
returnTimeProbabilities ::
    forall state.
    (FiniteState state) =>
    Natural ->
    TransitionMatrix state ->
    S.R (Cardinality state)
returnTimeProbabilities time p =
    S.vector (LA.toList latest)
  where
    dim = fromIntegral (natVal (Proxy @(Cardinality state)))
    matrix = S.extract (unTransitionMatrix p)
    initial = (LA.ident dim, LA.konst 0 dim)
    advance (survivors, _) = firstReturnStep matrix survivors
    (_, latest) = iterateNatural time advance initial

{- | Exact first-return probability @P(T_i^+ = t | X_0 = i)@ through any
'Transition'. Time zero is exactly zero; a self-loop returns at time one.
-}
returnTimeProbability ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    Natural ->
    kernel ->
    TransitionState kernel ->
    Double
returnTimeProbability 0 _ _ = 0
returnTimeProbability time kernel initialState =
    go time (Map.singleton initialState 1)
  where
    isInitial state = state == initialState

    go 0 _ = 0
    go _ survivors | Map.null survivors = 0
    go remaining survivors =
        let (next, returnMass) = advanceUntilTarget kernel isInitial survivors
         in if remaining == 1
                then returnMass
                else go (remaining - 1) next

{- | Strictly bounded first-return probabilities
@r_i(c) = P(T_i < c | X_0 = i)@ in state order.

Because first returns occur only at positive times, bounds @c = 0@ and @c = 1@
both give an all-zero vector. At @c = 2@ the result is the transition-matrix
diagonal. An absorbing state therefore has value exactly @1@ for every
@c >= 2@.

The implementation accumulates the exact first-return masses at times
@1, ..., c - 1@ while maintaining the shared survivor matrix described by
'returnTimeProbabilities'. In exact arithmetic,
@r_i(c + 1) - r_i(c) = P(T_i = c | X_0 = i)@. Results use ordinary 'Double'
arithmetic and are not clamped or renormalised.

Time: @O(c n^3)@; temporary space: @O(n^2)@; result space: @O(n)@.
-}
returnTimeBeforeProbabilities ::
    forall state.
    (FiniteState state) =>
    Natural ->
    TransitionMatrix state ->
    S.R (Cardinality state)
returnTimeBeforeProbabilities bound p =
    S.vector (LA.toList cumulative)
  where
    dim = fromIntegral (natVal (Proxy @(Cardinality state)))
    matrix = S.extract (unTransitionMatrix p)
    steps
        | bound == 0 = 0
        | otherwise = bound - 1
    initial = (LA.ident dim, LA.konst 0 dim)
    advance (survivors, total) =
        let (remaining, mass) = firstReturnStep matrix survivors
         in (remaining, total + mass)
    (_, cumulative) = iterateNatural steps advance initial

{- | Strict bounded first-return probability @P(T_i^+ < c | X_0 = i)@ through
any 'Transition'. Bounds @0@ and @1@ are exactly zero.
-}
returnTimeBeforeProbability ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    Natural ->
    kernel ->
    TransitionState kernel ->
    Double
returnTimeBeforeProbability bound kernel initialState =
    go bound (Map.singleton initialState 1) 0
  where
    isInitial state = state == initialState

    go remaining _ total | remaining <= 1 = total
    go _ survivors total | Map.null survivors = total
    go remaining survivors total =
        let (next, returnMass) = advanceUntilTarget kernel isInitial survivors
            cumulative = total + returnMass
         in cumulative `seq` go (remaining - 1) next cumulative

{- | First-return probabilities
@f_i = P(T_i < infinity | X_0 = i)@ in state order. Recurrent states are
exactly @1@ from support classification.

For all transient states, one fundamental-matrix solve computes
@N = (I - Q)^-1@ and @f_i = 1 - 1/N(i,i)@. Transient results inherit solver
rounding and are not clamped to @[0,1]@. Returns 'Left' if the transient system
fails the numerical contract.

Worst-case time: @O(n^3)@; temporary space: @O(n^2)@; result space: @O(n)@.
-}
returnProbabilities ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError (S.R (Cardinality state))
returnProbabilities p = do
    transientReturns <-
        if null transient
            then Right []
            else do
                nMatrix <- fundamental (subMatrix transientIdx transientIdx matrix)
                pure
                    [ 1 - 1 / (nMatrix `LA.atIndex` (k, k))
                    | k <- [0 .. length transient - 1]
                    ]
    let transientValues :: Unboxed.UArray Int Double
        transientValues =
            Unboxed.accumArray
                (\_ x -> x)
                0
                (0, dim - 1)
                (zip transientIdx transientReturns)
        valueAt i
            | recurrentState p i = 1
            | otherwise = transientValues Unboxed.! toIndex i
    pure (S.vector [valueAt i | i <- finiteStates])
  where
    dim = fromIntegral (natVal (Proxy @(Cardinality state)))
    transient = transientStates p
    transientIdx = map toIndex transient
    matrix = S.extract (unTransitionMatrix p)

{- | The probability of returning to one state after at least one transition.
A recurrent-state query returns exactly @1@ without forcing the
fundamental-matrix solve. The first transient query costs @O(n^3)@ worst case;
partial application shares that solve, making later lookups @O(1)@.

Transient queries inherit the numerical behavior and errors of
'returnProbabilities'.
-}
returnProbability ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError Double
returnProbability p =
    \i ->
        if recurrentState p i
            then Right 1
            else (`LA.atIndex` toIndex i) <$> probabilities
  where
    probabilities = S.extract <$> returnProbabilities p

{- | Expected first-return times @E(T_i | X_0 = i)@ in state order. Transient
states are exactly 'InfiniteExpectation'. Each recurrent state uses the
first-step identity
@m_i = 1 + sum_j P(i,j) eta_j({i})@, where @eta_j({i})@ is the expected
hitting time of the singleton target @{i}@ from @j@.

The implementation performs a separate singleton hitting-time solve for each
recurrent state. For @r@ recurrent states, time is @O(n^2 + r n^3)@, at most
@O(n^4)@; each solve uses @O(n^2)@ temporary space and the result uses
@O(n)@ space.
-}
returnTimeExpectations ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError [Expectation]
-- Reuse the hitting-time path; an all-state Kac calculation would require
-- stationary-distribution machinery not otherwise present in this module.
returnTimeExpectations p =
    traverse (returnTimeExpectation p) finiteStates

{- | The expected first-return time for one state. A transient state returns
'InfiniteExpectation' without a linear solve. A recurrent state performs one
singleton hitting-time solve and uses only stored transition probabilities
greater than zero; zero and tolerated negative entries contribute nothing.

A recurrent query takes @O(n^3)@ time and @O(n^2)@ temporary space. After
classification is cached, a transient query takes @O(1)@. Recurrent queries
inherit the numerical behavior and errors of 'Hit.expectation'.
-}
returnTimeExpectation ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError Expectation
returnTimeExpectation p i
    | recurrentState p i = returnTimeExpectationFrom p i
    | otherwise = Right InfiniteExpectation

returnTimeExpectationFrom ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError Expectation
returnTimeExpectationFrom p i =
    foldM addTerm (FiniteExpectation 1) (zip finiteStates row)
  where
    eta = Hit.expectation p [i]
    row = LA.toList (S.extract (unDistributionVector (rowAt p i)))
    addTerm acc (j, pij)
        | pij <= 0 = Right acc
        | otherwise = do
            hitting <- eta j
            pure $
                case (acc, hitting) of
                    (FiniteExpectation total, FiniteExpectation hittingTime) ->
                        FiniteExpectation (total + pij * hittingTime)
                    _ -> InfiniteExpectation

-- Direct survivor mass @P(T_i > t)@ through a locally finite transition.
-- Unlike hitting time, the initial state is not removed at time zero: a first
-- return can occur only after at least one transition.
returnTimeAfterProbability ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    Natural ->
    kernel ->
    TransitionState kernel ->
    Double
returnTimeAfterProbability time kernel initialState =
    go time (Map.singleton initialState 1)
  where
    isInitial state = state == initialState
    go 0 survivors = sum (Map.elems survivors)
    go _ survivors | Map.null survivors = 0
    go remaining survivors =
        let (next, _) = advanceUntilTarget kernel isInitial survivors
         in next `seq` go (remaining - 1) next

-- All-state first-return survivor mass. Row @i@ of the survivor matrix tracks
-- paths started at @i@ that have not revisited @i@ at a positive time.
returnTimeAfterProbabilities ::
    forall state.
    (FiniteState state) =>
    Natural ->
    TransitionMatrix state ->
    S.R (Cardinality state)
returnTimeAfterProbabilities time p =
    S.vector (LA.toList (survivors LA.#> LA.konst 1 dim))
  where
    dim = fromIntegral (natVal (Proxy @(Cardinality state)))
    matrix = S.extract (unTransitionMatrix p)
    advance current = fst (firstReturnStep matrix current)
    survivors = iterateNatural time advance (LA.ident dim)

{- | Probability of a finite-threshold event in the first-return time
@T_i = inf { t >= 1 | X_t = i }@.

'EqualTo' and the lower tails reuse the direct exact/bounded recurrences.
'GreaterThan' and 'AtLeast' use surviving mass directly, include the atom at
infinity, and avoid complement cancellation. Because time zero is excluded,
'EqualTo' @0@, 'LessThan' @1@, and 'AtMost' @0@ are exactly zero, while
'GreaterThan' @0@, 'AtLeast' @0@, and 'AtLeast' @1@ are exactly one.

The scalar query works through any locally finite 'Transition'. Results use
ordinary 'Double' arithmetic without clamping or renormalisation.
-}
probability ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    DiscreteEvent ->
    kernel ->
    TransitionState kernel ->
    Double
probability event kernel initialState =
    case event of
        EqualTo time -> returnTimeProbability time kernel initialState
        LessThan bound -> returnTimeBeforeProbability bound kernel initialState
        AtMost time -> returnTimeBeforeProbability (time + 1) kernel initialState
        GreaterThan time -> returnTimeAfterProbability time kernel initialState
        AtLeast 0 -> 1
        AtLeast time -> returnTimeAfterProbability (time - 1) kernel initialState

{- | First-return event probabilities in canonical initial-state order.
Coordinate @i@ is @P_i(T_i in E)@ for the supplied 'DiscreteEvent' @E@.
Upper tails use the direct survivor matrix and include mass at infinity.

Time is @O(t n^3)@ for threshold @t@; temporary space is @O(n^2)@ and result
space is @O(n)@.
-}
probabilityByState ::
    forall state.
    (FiniteState state) =>
    DiscreteEvent ->
    TransitionMatrix state ->
    S.R (Cardinality state)
probabilityByState event matrix =
    case event of
        EqualTo time -> returnTimeProbabilities time matrix
        LessThan bound -> returnTimeBeforeProbabilities bound matrix
        AtMost time -> returnTimeBeforeProbabilities (time + 1) matrix
        GreaterThan time -> returnTimeAfterProbabilities time matrix
        AtLeast 0 -> S.vector [1 | _ <- finiteStates @state]
        AtLeast time -> returnTimeAfterProbabilities (time - 1) matrix

{- | Probability of eventually returning to one state after at least one
transition.
-}
eventualProbability ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError Double
eventualProbability = returnProbability

{- | Eventual first-return probabilities in canonical state order.
-}
eventualProbabilityByState ::
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError (S.R (Cardinality state))
eventualProbabilityByState = returnProbabilities

{- | Expected first-return time for one state. Transient states have
'InfiniteExpectation'.
-}
expectation ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError Expectation
expectation = returnTimeExpectation

{- | Expected first-return times in canonical state order.
-}
expectationByState ::
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError [Expectation]
expectationByState = returnTimeExpectations
