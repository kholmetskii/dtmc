{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

{- |
Module      : Dtmc.Analysis.VisitCount
Description : Finite- and infinite-horizon visit-count analysis.

Exact analysis of visit counts. The bounded functions count visits to a state
predicate before a strict time bound and work through any locally finite
'Transition'. The total-count functions analyze visits to one state over the
entire path of a finite 'TransitionMatrix'. Both notions include the initial
state at time zero.

For a target state @i@, the total visit count is
@V_i = sum_(t = 0)^infinity 1_{X_t = i}@. Its law is determined by the
probability of ever hitting @i@ and the probability of returning to @i@. The
implementation uses exact graph classification for zero and infinite cases,
and checked 'Double' linear solves for the remaining probabilities. It does
not simulate, truncate an infinite series, clamp, or renormalise results.
-}
module Dtmc.Analysis.VisitCount (
    LinearSystemError (..),
    Expectation (..),
    totalProbability,
    totalProbabilityGivenInitialState,
    infiniteProbability,
    infiniteProbabilityGivenInitialState,
    totalExpectation,
    totalExpectationGivenInitialState,
    boundedLaw,
    boundedProbability,
    boundedProbabilityGivenInitialState,
    boundedExpectation,
    boundedExpectationGivenInitialState,
    occupationMatrix,
) where

import Data.Finite (
    getFinite,
 )
import Data.Map.Strict qualified as Map
import Dtmc.Analysis.Absorption (
    fundamentalMatrix,
 )
import Dtmc.Analysis.Classification (
    accessible,
    recurrentState,
 )
import Dtmc.Analysis.Event (
    DiscreteEvent (..),
    matchesDiscreteEvent,
 )
import Dtmc.Analysis.Expectation (
    Expectation (..),
 )
import Dtmc.Analysis.HittingTime qualified as Hit
import Dtmc.Analysis.Initial.Internal (
    expectationUnderEither,
    probabilityUnderEither,
 )
import Dtmc.Analysis.LinearSystem (
    LinearSystemError (..),
 )
import Dtmc.Analysis.ReturnTime qualified as Return
import Dtmc.Distribution (
    Distribution (..),
 )
import Dtmc.Distribution.Map (
    pointMass,
 )
import Dtmc.Distribution.Map.Internal (
    DistributionMap (DistributionMap),
    unDistributionMap,
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
    TransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )

toIndex :: (FiniteState state) => state -> Int
toIndex = fromIntegral . getFinite . stateIndex

{- | Probabilities of exactly @n@ total visits in canonical state order.
The count is the first argument and the target is the third; coordinate @j@ is
@P(V_i = n | X_0 = j)@.

Writing @h_ji = P_j(H_i < infinity)@ and
@f_i = P_i(T_i^+ < infinity)@, a transient target has

* @P_j(V_i = 0) = 1 - h_ji@;
* @P_j(V_i = n) = h_ji f_i^(n - 1) (1 - f_i)@ for @n >= 1@;

For a recurrent target, every positive finite count has probability zero.
Recurrence is decided from the support graph, not by comparing a computed
return probability with one.

Structural zero cases avoid a linear solve. Other cases inherit the numerical
behavior and errors of the state-conditioned eventual hitting and return
queries.
-}
exactProbabilitiesByState ::
    forall state.
    (FiniteState state) =>
    Natural ->
    TransitionMatrix state ->
    state ->
    Either LinearSystemError (S.R (Cardinality state))
exactProbabilitiesByState count matrix target
    | count == 0 = mapProbabilities (1 -) <$> hitting
    | recurrentState matrix target = Right zeroProbabilities
    | otherwise = do
        hits <- hitting
        returning <- Return.eventualProbabilityGivenInitialState matrix target
        let finiteMass = returning ^ (count - 1) * (1 - returning)
        pure (mapProbabilities (* finiteMass) hits)
  where
    hitting :: Either LinearSystemError (S.R (Cardinality state))
    hitting =
        S.vector
            <$> traverse
                (Hit.eventualProbabilityGivenInitialState matrix [target])
                finiteStates
    zeroProbabilities = S.vector [0 | _ <- finiteStates @state]
    mapProbabilities transform =
        S.vector . map transform . LA.toList . S.extract

{- | Probabilities of infinitely many visits in canonical initial-state order.
For target @i@, coordinate @j@ is @P(V_i = infinity | X_0 = j)@. A recurrent
target returns its hitting probabilities; a transient target returns an exact
zero vector without a linear solve.

Recurrence is decided from the support graph. The recurrent case inherits the
numerical behavior and errors of
'Dtmc.Analysis.HittingTime.eventualProbabilityGivenInitialState'.
-}
infiniteProbabilitiesByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError (S.R (Cardinality state))
infiniteProbabilitiesByState matrix target
    | recurrentState matrix target =
        S.vector
            <$> traverse
                (Hit.eventualProbabilityGivenInitialState matrix [target])
                finiteStates
    | otherwise = Right (S.vector [0 | _ <- finiteStates @state])

{- | Probability of infinitely many visits from one initial state. Argument
order is matrix, target, then initial state. Partially applying the matrix and
target shares the all-state computation.
-}
infiniteProbabilityGivenInitialState ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    state ->
    Either LinearSystemError Double
infiniteProbabilityGivenInitialState matrix target =
    \initial -> (`LA.atIndex` toIndex initial) <$> probabilities
  where
    probabilities = S.extract <$> infiniteProbabilitiesByState matrix target

{- | Expected total visits to the target in canonical initial-state order.

For a transient target @i@, coordinate @j@ is
@h_ji / (1 - f_i)@. For a recurrent target it is zero when @i@ is
unreachable from @j@ and 'InfiniteExpectation' otherwise. The recurrent case
is decided entirely from the support graph and requires no linear solve.

Transient results inherit the numerical behavior and errors of
the state-conditioned eventual hitting and return queries.
-}
totalExpectationsByState ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError [Expectation]
totalExpectationsByState matrix target
    | recurrentState matrix target =
        Right
            [ if accessible matrix initial target
                then InfiniteExpectation
                else FiniteExpectation 0
            | initial <- finiteStates
            ]
    | otherwise = do
        hits <-
            traverse
                (Hit.eventualProbabilityGivenInitialState matrix [target])
                finiteStates
        returning <- Return.eventualProbabilityGivenInitialState matrix target
        pure
            [ FiniteExpectation (hit / (1 - returning))
            | hit <- hits
            ]

{- | Expected total visits to the target from one initial state. Argument
order is matrix, target, then initial state. Partially applying the matrix and
target shares the all-state computation.
-}
totalExpectationGivenInitialState ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    state ->
    Either LinearSystemError Expectation
totalExpectationGivenInitialState matrix target =
    \initial -> (!! toIndex initial) <$> expectations
  where
    expectations = totalExpectationsByState matrix target

iterateNatural :: Natural -> (value -> value) -> value -> value
iterateNatural steps advance = go steps
  where
    go 0 value = value
    go remaining value =
        let next = advance value
         in next `seq` go (remaining - 1) next

{- | Distribution of
@N_A(c) = sum_(t = 0)^(c - 1) 1_A(X_t)@, the number of visits to the supplied
state predicate strictly before time @c@.

The time bound is the first argument, consistently with the other
finite-horizon APIs. Bound zero returns a point mass at count zero. At a
positive bound, the initial state at time zero is included. Consequently the
result is supported on counts from zero through the bound.

The result is computed from the exact finite reachable support of the joint
process @(X_t, N_A(t + 1))@. Ordinary 'Double' arithmetic is preserved without
clamping or renormalisation.
-}
boundedLaw ::
    ( Distribution distribution
    , Transition transition
    , DistributionState distribution ~ TransitionState transition
    , Ord (TransitionState transition)
    ) =>
    Natural ->
    distribution ->
    transition ->
    (TransitionState transition -> Bool) ->
    DistributionMap Natural
boundedLaw bound initial transition isVisited
    | bound == 0 = DistributionMap (Map.singleton 0 1)
    | otherwise = DistributionMap (countMarginal finalJoint)
  where
    initialJoint =
        Map.fromList
            [ ((state, if isVisited state then 1 else 0), weight)
            | (state, weight) <- distributionWeights initial
            , weight /= 0
            ]
    finalJoint =
        iterateNatural (bound - 1) advanceJoint initialJoint

    advanceJoint joint =
        Map.filter (/= 0) (Map.foldlWithKey' advanceState Map.empty joint)

    advanceState accumulated (state, count) stateWeight =
        Map.foldlWithKey'
            (advanceDestination count stateWeight)
            accumulated
            (unDistributionMap (transitionLaw transition state))

    advanceDestination count stateWeight accumulated nextState transitionWeight =
        Map.insertWith
            (+)
            (nextState, count + if isVisited nextState then 1 else 0)
            (stateWeight * transitionWeight)
            accumulated

    countMarginal =
        Map.filter (/= 0)
            . Map.foldlWithKey'
                (\counts (_, count) weight -> Map.insertWith (+) count weight counts)
                Map.empty

{- | Expected number of visits before the strict time bound. Using
@E(N_A(c)) = sum_(t = 0)^(c - 1) P(X_t in A)@, this evolves only the state
marginal rather than constructing the joint count distribution. The result
lies mathematically between zero and the bound, subject to ordinary
floating-point error.
-}
boundedExpectation ::
    ( Distribution distribution
    , Transition transition
    , DistributionState distribution ~ TransitionState transition
    , Ord (TransitionState transition)
    ) =>
    Natural ->
    distribution ->
    transition ->
    (TransitionState transition -> Bool) ->
    Double
boundedExpectation bound initial transition isVisited =
    go bound (Map.fromList (distributionWeights initial)) 0
  where
    go 0 _ expectation = expectation
    go remaining weights expectation =
        let visitProbability =
                Map.foldlWithKey'
                    ( \total state weight ->
                        if isVisited state then total + weight else total
                    )
                    0
                    weights
            cumulative = expectation + visitProbability
         in if remaining == 1
                then cumulative
                else
                    let next = pushSparseWeights weights transition
                     in cumulative `seq` next `seq` go (remaining - 1) next cumulative

{- | Probability under an arbitrary initial distribution of a finite-threshold
event in the total number of visits
@V_i = sum_(t = 0)^infinity 1_{X_t = i}@.

The initial state at time zero is included. Upper-tail events include the atom
at infinity. For @h_ji = P_j(H_i < infinity)@ and transient-target return
probability @f_i@, the implementation evaluates
@P_j(V_i > n) = h_ji f_i^n@ directly. Recurrent targets use support
classification, so their positive mass is placed structurally at infinity.

The implementation mixes the internally shared state-conditioned results
under the supplied initial distribution.
-}
totalProbability ::
    ( FiniteState state
    , Distribution distribution
    , DistributionState distribution ~ state
    ) =>
    DiscreteEvent ->
    TransitionMatrix state ->
    state ->
    distribution ->
    Either LinearSystemError Double
totalProbability event matrix target initial =
    probabilityUnderEither initial (totalProbabilityGivenInitialState event matrix target)

-- | Probability of a total-visit event conditioned on @X_0 = j@.
totalProbabilityGivenInitialState ::
    (FiniteState state) =>
    DiscreteEvent ->
    TransitionMatrix state ->
    state ->
    state ->
    Either LinearSystemError Double
totalProbabilityGivenInitialState event matrix target =
    \initial -> (`LA.atIndex` toIndex initial) <$> probabilities
  where
    probabilities = S.extract <$> totalProbabilityByState event matrix target

{- | Total-visit event probabilities in canonical initial-state order.
Coordinate @j@ is @P_j(V_i in E)@ for the supplied target @i@ and
'DiscreteEvent' @E@. Upper tails are evaluated directly and include infinitely
many visits; lower tails contain finite counts only.

Graph classification and any required checked linear solves are shared across
all initial states. Structural zero and one boundaries avoid a solve.
-}
totalProbabilityByState ::
    forall state.
    (FiniteState state) =>
    DiscreteEvent ->
    TransitionMatrix state ->
    state ->
    Either LinearSystemError (S.R (Cardinality state))
totalProbabilityByState event matrix target =
    case event of
        EqualTo count -> exactProbabilitiesByState count matrix target
        LessThan 0 -> Right zeros
        LessThan bound -> atMost (bound - 1)
        AtMost count -> atMost count
        GreaterThan count -> after count
        AtLeast 0 -> Right ones
        AtLeast count -> after (count - 1)
  where
    recurrent = recurrentState matrix target
    hits :: Either LinearSystemError (S.R (Cardinality state))
    hits =
        S.vector
            <$> traverse
                (Hit.eventualProbabilityGivenInitialState matrix [target])
                finiteStates
    zeros = S.vector [0 | _ <- finiteStates @state]
    ones = S.vector [1 | _ <- finiteStates @state]
    mapValues transform =
        S.vector . map transform . LA.toList . S.extract

    atMost count
        | recurrent = mapValues (1 -) <$> hits
        | otherwise = do
            hitValues <- hits
            returning <- Return.eventualProbabilityGivenInitialState matrix target
            pure
                ( mapValues
                    (\hit -> 1 - hit * returning ^ count)
                    hitValues
                )

    after count
        | recurrent = hits
        | count == 0 = hits
        | otherwise = do
            hitValues <- hits
            returning <- Return.eventualProbabilityGivenInitialState matrix target
            pure (mapValues (* (returning ^ count)) hitValues)

{- | Probability under an arbitrary initial distribution of infinitely many
visits.
-}
infiniteProbability ::
    ( FiniteState state
    , Distribution distribution
    , DistributionState distribution ~ state
    ) =>
    TransitionMatrix state ->
    state ->
    distribution ->
    Either LinearSystemError Double
infiniteProbability matrix target initial =
    probabilityUnderEither initial (infiniteProbabilityGivenInitialState matrix target)

-- | Expected total visits under an arbitrary initial distribution.
totalExpectation ::
    ( FiniteState state
    , Distribution distribution
    , DistributionState distribution ~ state
    ) =>
    TransitionMatrix state ->
    state ->
    distribution ->
    Either LinearSystemError Expectation
totalExpectation matrix target initial =
    expectationUnderEither initial (totalExpectationGivenInitialState matrix target)

{- | Probability of a 'DiscreteEvent' in the number of visits strictly before
a finite time bound. The bounded law has no atom at infinity, and values
outside its support contribute exactly zero.
-}
boundedProbability ::
    ( Distribution distribution
    , Transition transition
    , DistributionState distribution ~ TransitionState transition
    , Ord (TransitionState transition)
    ) =>
    Natural ->
    DiscreteEvent ->
    distribution ->
    transition ->
    (TransitionState transition -> Bool) ->
    Double
boundedProbability bound event initial transition isVisited =
    sum
        [ weight
        | (count, weight) <- distributionWeights law
        , matchesDiscreteEvent event count
        ]
  where
    law = boundedLaw bound initial transition isVisited

-- | Probability of a bounded visit-count event conditioned on @X_0 = i@.
boundedProbabilityGivenInitialState ::
    ( Transition transition
    , Ord (TransitionState transition)
    ) =>
    Natural ->
    DiscreteEvent ->
    TransitionState transition ->
    transition ->
    (TransitionState transition -> Bool) ->
    Double
boundedProbabilityGivenInitialState bound event initial =
    boundedProbability bound event (pointMass initial)

{- | Expected visits before a strict finite time bound conditioned on
@X_0 = i@.
-}
boundedExpectationGivenInitialState ::
    ( Transition transition
    , Ord (TransitionState transition)
    ) =>
    Natural ->
    TransitionState transition ->
    transition ->
    (TransitionState transition -> Bool) ->
    Double
boundedExpectationGivenInitialState bound initial =
    boundedExpectation bound (pointMass initial)

{- | The occupation matrix of the chain, also known as its Green function:
entry @(i, j)@ is @sum_(n >= 0) (P^n)(i,j) = E(V_j | X_0 = i)@, the expected
total number of visits to @j@ started from @i@. Rows and columns follow the
canonical order of the 'FiniteState' instance.

Unlike 'Dtmc.Analysis.Absorption.fundamentalMatrix', which is the finite
@T x T@ block, this is defined on the whole state space and therefore needs
'Expectation': a recurrent target reachable from @i@ is visited infinitely
often almost surely. The four cases are

* @j@ transient and @i@ transient: the corresponding entry of
  @(I - Q)^-1@;
* @j@ transient and @i@ recurrent: exactly zero, because a recurrent class is
  closed and cannot reach a transient state;
* @j@ recurrent and reachable from @i@: 'InfiniteExpectation';
* @j@ recurrent and unreachable from @i@: exactly zero.

Only the transient block needs arithmetic; the infinite and zero entries come
from the support graph, so they are exact.
'Dtmc.Analysis.VisitCount.totalExpectation' computes single entries by a
different route and agrees with this one.

Time: @O(n^2 + t^3 + n E)@ for @t@ transient states and @E@ support edges.
Result space: @O(n^2)@.
-}
occupationMatrix ::
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError [[Expectation]]
occupationMatrix p = do
    (transient, block) <- fundamentalMatrix p
    let table =
            Map.fromList
                [ ((i, j), value)
                | (i, row) <- zip transient block
                , (j, value) <- zip transient row
                ]
        valueAt i j
            | recurrentState p j =
                if accessible p i j
                    then InfiniteExpectation
                    else FiniteExpectation 0
            | otherwise =
                FiniteExpectation (Map.findWithDefault 0 (i, j) table)
    pure [[valueAt i j | j <- finiteStates] | i <- finiteStates]
