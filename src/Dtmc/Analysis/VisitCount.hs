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
    visitCountProbabilities,
    visitCountProbability,
    infiniteVisitProbabilities,
    infiniteVisitProbability,
    visitCountExpectations,
    visitCountExpectation,
    visitCountDistributionBefore,
    visitCountProbabilityBefore,
    visitCountExpectationBefore,
) where

import Data.Finite (
    getFinite,
 )
import Data.Map.Strict qualified as Map
import Dtmc.Analysis.Classification (
    accessible,
    recurrentState,
 )
import Dtmc.Analysis.Expectation (
    Expectation (..),
 )
import Dtmc.Analysis.HittingTime (
    hittingProbabilities,
 )
import Dtmc.Analysis.LinearSystem (
    LinearSystemError (..),
 )
import Dtmc.Analysis.ReturnTime (
    returnProbability,
 )
import Dtmc.Distribution (
    Distribution (..),
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
behavior and errors of 'hittingProbabilities' and 'returnProbability'.
-}
visitCountProbabilities ::
    forall state.
    (FiniteState state) =>
    Natural ->
    TransitionMatrix state ->
    state ->
    Either LinearSystemError (S.R (Cardinality state))
visitCountProbabilities count matrix target
    | count == 0 = mapProbabilities (1 -) <$> hitting
    | recurrentState matrix target = Right zeroProbabilities
    | otherwise = do
        hits <- hitting
        returning <- returnProbability matrix target
        let finiteMass = returning ^ (count - 1) * (1 - returning)
        pure (mapProbabilities (* finiteMass) hits)
  where
    hitting = hittingProbabilities matrix [target]
    zeroProbabilities = S.vector [0 | _ <- finiteStates @state]
    mapProbabilities transform =
        S.vector . map transform . LA.toList . S.extract

{- | Probability of exactly @n@ total visits from one initial state. Argument
order is count, matrix, target, then initial state, matching the exact-time
hitting and return APIs. Partially applying the first three arguments shares
the all-state computation.
-}
visitCountProbability ::
    (FiniteState state) =>
    Natural ->
    TransitionMatrix state ->
    state ->
    state ->
    Either LinearSystemError Double
visitCountProbability count matrix target =
    \initial -> (`LA.atIndex` toIndex initial) <$> probabilities
  where
    probabilities = S.extract <$> visitCountProbabilities count matrix target

{- | Probabilities of infinitely many visits in canonical initial-state order.
For target @i@, coordinate @j@ is @P(V_i = infinity | X_0 = j)@. A recurrent
target returns its hitting probabilities; a transient target returns an exact
zero vector without a linear solve.

Recurrence is decided from the support graph. The recurrent case inherits the
numerical behavior and errors of 'hittingProbabilities'.
-}
infiniteVisitProbabilities ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError (S.R (Cardinality state))
infiniteVisitProbabilities matrix target
    | recurrentState matrix target = hittingProbabilities matrix [target]
    | otherwise = Right (S.vector [0 | _ <- finiteStates @state])

{- | Probability of infinitely many visits from one initial state. Argument
order is matrix, target, then initial state. Partially applying the matrix and
target shares the all-state computation.
-}
infiniteVisitProbability ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    state ->
    Either LinearSystemError Double
infiniteVisitProbability matrix target =
    \initial -> (`LA.atIndex` toIndex initial) <$> probabilities
  where
    probabilities = S.extract <$> infiniteVisitProbabilities matrix target

{- | Expected total visits to the target in canonical initial-state order.

For a transient target @i@, coordinate @j@ is
@h_ji / (1 - f_i)@. For a recurrent target it is zero when @i@ is
unreachable from @j@ and 'InfiniteExpectation' otherwise. The recurrent case
is decided entirely from the support graph and requires no linear solve.

Transient results inherit the numerical behavior and errors of
'hittingProbabilities' and 'returnProbability'.
-}
visitCountExpectations ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError [Expectation]
visitCountExpectations matrix target
    | recurrentState matrix target =
        Right
            [ if accessible matrix initial target
                then InfiniteExpectation
                else FiniteExpectation 0
            | initial <- finiteStates
            ]
    | otherwise = do
        hits <- hittingProbabilities matrix [target]
        returning <- returnProbability matrix target
        pure
            [ FiniteExpectation (hit / (1 - returning))
            | hit <- LA.toList (S.extract hits)
            ]

{- | Expected total visits to the target from one initial state. Argument
order is matrix, target, then initial state. Partially applying the matrix and
target shares the all-state computation.
-}
visitCountExpectation ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    state ->
    Either LinearSystemError Expectation
visitCountExpectation matrix target =
    \initial -> (!! toIndex initial) <$> expectations
  where
    expectations = visitCountExpectations matrix target

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
visitCountDistributionBefore ::
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
visitCountDistributionBefore bound initial transition isVisited
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

{- | Probability of exactly the requested number of visits before the strict
time bound. The requested count is the final argument, mirroring a lookup in
'visitCountDistributionBefore'. Counts outside the distribution's support
return exactly zero.
-}
visitCountProbabilityBefore ::
    ( Distribution distribution
    , Transition transition
    , DistributionState distribution ~ TransitionState transition
    , Ord (TransitionState transition)
    ) =>
    Natural ->
    distribution ->
    transition ->
    (TransitionState transition -> Bool) ->
    Natural ->
    Double
visitCountProbabilityBefore bound initial transition isVisited =
    probabilityAt
        (visitCountDistributionBefore bound initial transition isVisited)

{- | Expected number of visits before the strict time bound. Using
@E(N_A(c)) = sum_(t = 0)^(c - 1) P(X_t in A)@, this evolves only the state
marginal rather than constructing the joint count distribution. The result
lies mathematically between zero and the bound, subject to ordinary
floating-point error.
-}
visitCountExpectationBefore ::
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
visitCountExpectationBefore bound initial transition isVisited =
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
