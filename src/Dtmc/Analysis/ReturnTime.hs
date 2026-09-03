{- |
Module      : Dtmc.Analysis.ReturnTime
Description : Exact, bounded, eventual, and expected first-return times.

First-return quantities for DTMCs. For state @i@,
@T_i = inf { t >= 1 | X_t = i }@, so time zero is never a return. Scalar
exact-time and bounded queries work through any locally finite 'Transition'.
Eventual and expected queries use a finite 'TransitionMatrix'.

Exact-time and strictly bounded queries use finite recurrences. Eventual
queries use support classification and checked 'Double' linear solves.
Expected return times use class stationary distributions and Kac's formula.
Results are not clamped or renormalised.
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

{- | Exact first-return probability @P(T_i^+ = t | X_0 = i)@ through any
'Transition'. Time zero is exactly zero; a self-loop returns at time one.
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

{- | Strict bounded first-return probability @P(T_i^+ < c | X_0 = i)@ through
any 'Transition'. Bounds @0@ and @1@ are exactly zero.
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

{- | First-return probabilities
@f_i = P(T_i < infinity | X_0 = i)@ in state order. Recurrent states are
exactly @1@ from support classification.

For all transient states, one fundamental-matrix solve computes
@N = (I - Q)^-1@ and @f_i = 1 - 1/N(i,i)@. Transient results inherit solver
rounding and are not clamped to @[0,1]@. Returns 'Left' if the transient system
fails the numerical contract.

Worst-case time: @O(n^3)@; temporary space: @O(n^2)@; result space: @O(n)@.
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

{- | The probability of returning to one state after at least one transition.
A recurrent-state query returns exactly @1@ without forcing the
fundamental-matrix solve. The first transient query costs @O(n^3)@ worst case;
partial application shares that solve, making later lookups @O(1)@.

Transient queries inherit the numerical behavior and errors of
@eventualProbabilitiesByState@.
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

{- | The expected first-return time for one state. A transient state returns
'InfiniteExpectation' without a numerical solve. For a recurrent state @i@,
Kac's formula gives @E_i T_i^+ = 1 / pi_i@, where @pi@ is the stationary
distribution of its closed communicating class.

The first recurrent query computes the stationary distributions of all closed
classes in @O(n^3)@ worst-case time and @O(n^2)@ temporary space. Partial
application shares that table, making later lookups @O(1)@. A transient query
does not force the table. Numerical failures are those of
'stationaryDistributions'.
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

{- | Expected return times for all recurrent states from one set of
class-stationary solves. The array also contains infinity at transient
coordinates, although callers decide transience structurally before lookup.
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

{- | Probability of a finite-threshold event in the first-return time
@T_i = inf { t >= 1 | X_t = i }@.

'EqualTo' and the lower tails reuse the direct exact/bounded recurrences.
'GreaterThan' and 'AtLeast' use surviving mass directly, include the atom at
infinity, and avoid complement cancellation. Because time zero is excluded,
'EqualTo' @0@, 'LessThan' @1@, and 'AtMost' @0@ are exactly zero, while
'GreaterThan' @0@, 'AtLeast' @0@, and 'AtLeast' @1@ are exactly one.

The initial state is sampled from the supplied distribution; each path then
measures return to its own sampled state. The query works through any locally
finite 'Transition'. Results use ordinary 'Double' arithmetic without
clamping or renormalisation.
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

{- | Probability of a finite-threshold first-return event conditioned on
@X_0 = i@. The return target is that same initial state @i@.
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

{- | Probability under an arbitrary initial distribution of eventually
returning to the sampled initial state after at least one transition.
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

{- | Expected first-return time under an arbitrary initial distribution. The
result is infinite when a transient state has positive initial probability.
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
