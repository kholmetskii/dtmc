{- |
Module      : Dtmc.Analysis.ReturnTime
Description : Exact, bounded, eventual, and expected first-return times.

First-return quantities for DTMCs. For state @i@,
@T_i^+ = inf { t >= 1 | X_t = i }@, so time zero is never a return. Scalar
exact-time and bounded queries work through any locally finite 'Transition'.
Eventual and expected queries use a finite 'TransitionMatrix'.

Exact-time and strictly bounded queries use finite recurrences. Eventual
queries use support classification and checked 'Double' linear solves.
Expected return times use class stationary distributions and Kac's formula.
Results are not clamped or renormalised.

Unless stated otherwise, complexity bounds exclude 'FiniteState' method
costs. Bounds over abstract distributions or transitions also identify the
excluded typeclass-method costs.

For finite-matrix bounds, @n@ is the state count and @E@ the support-edge
count.
-}
module Dtmc.Analysis.ReturnTime (
    LinearSystemError (..),
    Expectation (..),
    probability,
    probabilityGivenInitialState,
    eventualProbability,
    eventualProbabilityGivenInitialState,
    expectation,
    expectationGivenInitialState,
) where

import Data.Array qualified as Array
import Data.Array.Unboxed qualified as Unboxed
import Data.Map.Strict (
    Map,
 )
import Data.Map.Strict qualified as Map
import Dtmc.Analysis.Classification (
    recurrentState,
    transientStates,
 )
import Dtmc.Analysis.Event (
    DiscreteEvent (..),
 )
import Dtmc.Analysis.Expectation (
    Expectation (..),
 )
import Dtmc.Analysis.Initial.Internal (
    expectationUnderEither,
    probabilityUnder,
    probabilityUnderEither,
 )
import Dtmc.Analysis.LinearSystem (
    LinearSystemError (..),
 )
import Dtmc.Analysis.LinearSystem.Internal (
    fundamental,
    subMatrix,
 )
import Dtmc.Analysis.Stationary (
    stationaryDistributions,
 )
import Dtmc.Distribution (
    Distribution (..),
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
 )
import Dtmc.State.Internal (
    stateCardinalityInt,
    stateIndexInt,
 )
import Dtmc.Transition (
    Transition (..),
 )
import Dtmc.Transition.Matrix.Internal (
    TransitionMatrix,
    unTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )

toIndex :: (FiniteState state) => state -> Int
toIndex = stateIndexInt

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

{- | Compute the exact first-return probability
@P(T_i^+ = t | X_0 = i)@ through any 'Transition'. Time zero is exactly zero;
a self-loop returns at time one.

Complexity: excluding 'transitionLaw',
@O(k (w + e log(u + 1) + u) + 1)@ time, @O(w + u)@ temporary space, and
@O(1)@ result space, where @k = t@ and @w@, @e@, and @u@ bound per-step
survivor states, traversed transition edges, and accumulated destinations.
-}
exactProbabilityAt ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    Natural ->
    kernel ->
    TransitionState kernel ->
    Double
exactProbabilityAt 0 _ _ = 0
exactProbabilityAt time kernel initialState =
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

{- | Compute the strict bounded first-return probability
@P(T_i^+ < c | X_0 = i)@ through any 'Transition'. Bounds @0@ and @1@ are
exactly zero.

Complexity: excluding 'transitionLaw',
@O(k (w + e log(u + 1) + u) + 1)@ time, @O(w + u)@ temporary space, and
@O(1)@ result space, where @k = c@ and @w@, @e@, and @u@ bound per-step
survivor states, traversed transition edges, and accumulated destinations.
-}
lowerTailProbability ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    Natural ->
    kernel ->
    TransitionState kernel ->
    Double
lowerTailProbability bound kernel initialState =
    go bound (Map.singleton initialState 1) 0
  where
    isInitial state = state == initialState

    go remaining _ total | remaining <= 1 = total
    go _ survivors total | Map.null survivors = total
    go remaining survivors total =
        let (next, returnMass) = advanceUntilTarget kernel isInitial survivors
            cumulative = total + returnMass
         in cumulative `seq` go (remaining - 1) next cumulative

{- | Compute first-return probabilities
@f_i = P(T_i^+ < infinity | X_0 = i)@ in state order. Recurrent states are
exactly @1@ from support classification.

For all transient states, one fundamental-matrix solve computes
@N = (I - Q)^-1@ and @f_i = 1 - 1/N(i,i)@. Transient results inherit solver
rounding and are not clamped to @[0,1]@. Returns 'Left' if the transient system
fails the numerical contract.

Complexity: @O(n^3)@ worst-case time, @O(n^2)@ temporary space,
@O(n + E)@ retained graph-cache space, and @O(n)@ result space for @n@
states and @E@ support edges.
-}
eventualProbabilitiesByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError (S.R (Cardinality state))
eventualProbabilitiesByState p = do
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
    dim = stateCardinalityInt @state
    transient = transientStates p
    transientIdx = map toIndex transient
    matrix = S.extract (unTransitionMatrix p)

{- | Compute the probability of returning to one state after at least one
transition. A recurrent-state query returns exactly @1@ without forcing the
fundamental-matrix solve. Partial application shares the all-state transient
solve.

Transient queries inherit the numerical behaviour and errors of
@eventualProbabilitiesByState@.

Complexity: the first transient query takes @O(n^3)@ worst-case time and
@O(n^2)@ temporary space and may retain an @O(n)@ all-state result; later
lookups take @O(1)@ time and space. A recurrent query avoids the solve. The
matrix may retain @O(n + E)@ graph-cache space, and the scalar result occupies
@O(1)@ space.
-}
eventualProbabilityGivenInitialState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError Double
eventualProbabilityGivenInitialState p =
    \i ->
        if recurrentState p i
            then Right 1
            else (`LA.atIndex` toIndex i) <$> probabilities
  where
    probabilities = S.extract <$> eventualProbabilitiesByState p

{- | Compute the expected first-return time for one state. A transient state
returns 'InfiniteExpectation' without a numerical solve. For a recurrent
state @i@, Kac's formula gives @E_i T_i^+ = 1 / pi_i@, where @pi@ is the
stationary distribution of its closed communicating class.

Partial application shares the stationary distributions and resulting
all-state table. A transient query does not force that table. Numerical
failures are those of 'stationaryDistributions' or a non-positive or
non-finite recurrent stationary probability.

Complexity: the first recurrent query takes @O(n^3)@ worst-case time and
@O(n^2)@ temporary space and may retain an @O(n)@ all-state result; later
lookups take @O(1)@ time and space. A transient query avoids the stationary
solves. The matrix may retain @O(n + E)@ graph-cache space, and the scalar
result occupies @O(1)@ space.
-}
expectationGivenInitialState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError Expectation
expectationGivenInitialState p =
    \i ->
        if recurrentState p i
            then (Array.! toIndex i) <$> recurrentExpectations
            else Right InfiniteExpectation
  where
    recurrentExpectations = recurrentReturnExpectations p

{- | Compute expected return times for all recurrent states from one set of
class-stationary solves. The array also contains infinity at transient
coordinates, although callers decide transience structurally before lookup.

Complexity: @O(n^3)@ worst-case time, @O(n^2)@ temporary space,
@O(n + E)@ retained graph-cache space, and @O(n)@ result space.
-}
recurrentReturnExpectations ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError (Array.Array Int Expectation)
recurrentReturnExpectations p = do
    classes <- stationaryDistributions p
    recurrentEntries <- concat <$> traverse entriesForClass classes
    pure
        ( Array.accumArray
            (\_ value -> value)
            InfiniteExpectation
            (0, stateCardinalityInt @state - 1)
            recurrentEntries
        )
  where
    entriesForClass (members, distribution) =
        traverse (entry vector) members
      where
        vector = S.extract (unDistributionVector distribution)

    entry vector member
        | stationaryProbability <= 0 || not (finite reciprocal) = Left NonFiniteSolution
        | otherwise = Right (index, FiniteExpectation reciprocal)
      where
        index = toIndex member
        stationaryProbability = vector `LA.atIndex` index
        reciprocal = 1 / stationaryProbability

    finite value = not (isNaN value || isInfinite value)

-- Direct survivor mass @P(T_i > t)@ through a locally finite transition.
-- Unlike hitting time, the initial state is not removed at time zero: a first
-- return can occur only after at least one transition.
upperTailProbability ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    Natural ->
    kernel ->
    TransitionState kernel ->
    Double
upperTailProbability time kernel initialState =
    go time (Map.singleton initialState 1)
  where
    isInitial state = state == initialState
    go 0 survivors = sum (Map.elems survivors)
    go _ survivors | Map.null survivors = 0
    go remaining survivors =
        let (next, _) = advanceUntilTarget kernel isInitial survivors
         in next `seq` go (remaining - 1) next

{- | Compute the probability of a finite-threshold event in the first-return
time @T_i^+ = inf { t >= 1 | X_t = i }@.

'EqualTo' and the lower tails reuse the direct exact/bounded recurrences.
'GreaterThan' and 'AtLeast' use surviving mass directly, include the atom at
infinity, and avoid complement cancellation. Because time zero is excluded,
'EqualTo' @0@, 'LessThan' @1@, and 'AtMost' @0@ are exactly zero, while
'GreaterThan' @0@, 'AtLeast' @0@, and 'AtLeast' @1@ are exactly one.

The initial state is sampled from the supplied distribution; each path then
measures return to its own sampled state. The query works through any locally
finite 'Transition'. Results use ordinary 'Double' arithmetic without
clamping or renormalisation.

For the complexity bounds, @s@ is the initial stored support size, @k@ the
event threshold, and @w@, @e@, and @u@ are per-step upper bounds on survivor
states, traversed transition edges, and accumulated destinations.

Complexity: excluding 'distributionWeights' and 'transitionLaw',
@O(s (k (w + e log(u + 1) + u) + 1))@ time, @O(s + w + u)@ temporary
space, and @O(1)@ result space.
-}
probability ::
    ( Distribution distribution
    , Transition kernel
    , DistributionState distribution ~ TransitionState kernel
    , Ord (TransitionState kernel)
    ) =>
    DiscreteEvent ->
    kernel ->
    distribution ->
    Double
probability event kernel initial =
    probabilityUnder initial (probabilityGivenInitialState event kernel)

{- | Compute the probability of a finite-threshold first-return event
conditioned on @X_0 = i@. The return target is that same initial state, and
time zero is not a return.

For the complexity bounds, @k@ is the event threshold and @w@, @e@, and @u@
are per-step upper bounds on survivor states, traversed transition edges, and
accumulated destinations.

Complexity: excluding 'transitionLaw',
@O(k (w + e log(u + 1) + u) + 1)@ time, @O(w + u)@ temporary space, and
@O(1)@ result space.
-}
probabilityGivenInitialState ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    DiscreteEvent ->
    kernel ->
    TransitionState kernel ->
    Double
probabilityGivenInitialState event kernel initialState =
    case event of
        EqualTo time -> exactProbabilityAt time kernel initialState
        LessThan bound -> lowerTailProbability bound kernel initialState
        AtMost time -> lowerTailProbability (time + 1) kernel initialState
        GreaterThan time -> upperTailProbability time kernel initialState
        AtLeast 0 -> 1
        AtLeast time -> upperTailProbability (time - 1) kernel initialState

{- | Compute, under an arbitrary initial distribution, the probability of
eventually returning to the sampled initial state after at least one
transition. Recurrent states contribute exactly one; transient-state values
come from one shared checked fundamental-matrix solve.

Complexity: excluding 'distributionWeights', @O(n^3 + s)@ worst-case time,
@O(n^2 + s)@ temporary space, and @O(1)@ result space for @n@ states and
initial stored support size @s@. The matrix may retain @O(n + E)@ graph-cache
space.
-}
eventualProbability ::
    ( FiniteState state
    , Distribution distribution
    , DistributionState distribution ~ state
    ) =>
    TransitionMatrix state ->
    distribution ->
    Either LinearSystemError Double
eventualProbability matrix initial =
    probabilityUnderEither initial (eventualProbabilityGivenInitialState matrix)

{- | Compute the expected first-return time under an arbitrary initial
distribution. The result is infinite when a transient state has positive
initial probability. Otherwise recurrent-state expectations use shared class
stationary distributions and Kac's formula.

Complexity: excluding 'distributionWeights', @O(n^3 + s)@ worst-case time,
@O(n^2 + s)@ temporary space, and @O(1)@ result space for @n@ states and
initial stored support size @s@. The matrix may retain @O(n + E)@ graph-cache
space.
-}
expectation ::
    ( FiniteState state
    , Distribution distribution
    , DistributionState distribution ~ state
    ) =>
    TransitionMatrix state ->
    distribution ->
    Either LinearSystemError Expectation
expectation matrix initial =
    expectationUnderEither initial (expectationGivenInitialState matrix)
