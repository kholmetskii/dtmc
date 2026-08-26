{- |
Module      : Dtmc.Analysis.VisitCount
Description : Finite-horizon visit-count distributions, probabilities, and expectations.

Exact analysis of the number of visits to a state predicate before a strict
time bound. The implementation tracks the finite-support joint law of the
current state and accumulated count, so it works with finite matrices and
locally finite kernels without simulation, truncation, or global state-space
enumeration.
-}
module Dtmc.Analysis.VisitCount (
    visitCountDistributionBefore,
    visitCountProbabilityBefore,
    visitCountExpectationBefore,
) where

import Data.Map.Strict qualified as Map
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
import Dtmc.Transition (
    Transition (..),
 )
import Numeric.Natural (
    Natural,
 )

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
visitCountProbabilityBefore bound initial transition isVisited count =
    probabilityAt
        (visitCountDistributionBefore bound initial transition isVisited)
        count

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
