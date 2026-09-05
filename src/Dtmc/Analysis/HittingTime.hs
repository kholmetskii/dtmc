{- |
Module      : Dtmc.Analysis.HittingTime
Description : Exact, bounded, eventual, competing, and expected hitting times.

Hitting-time quantities for DTMCs. Scalar exact-time and bounded queries work
through any locally finite 'Transition', including kernels on infinite state
spaces. Eventual, competing, and expected queries use a finite
'TransitionMatrix'. For a target set @A@,
@H_A = inf { t >= 0 | X_t in A }@.

Exact-time and strictly bounded queries use finite recurrences. Eventual,
competing, and expected queries use support reachability and checked 'Double'
linear solves. Results are not clamped or renormalised.

Unless stated otherwise, complexity bounds exclude 'FiniteState' method
costs. Bounds over abstract distributions, transitions, or target predicates
also identify the excluded typeclass-method and predicate costs.

For finite-matrix bounds, @n@ is the state count and @E@ the support-edge
count.
-}
module Dtmc.Analysis.HittingTime (
    LinearSystemError (..),
    Expectation (..),
    probability,
    probabilityGivenInitialState,
    eventualProbability,
    eventualProbabilityGivenInitialState,
    raceProbability,
    raceProbabilityGivenInitialState,
    expectation,
    expectationGivenInitialState,
) where

import Data.Array qualified as Array
import Data.Array.Unboxed qualified as Unboxed
import Data.Map.Strict (
    Map,
 )
import Data.Map.Strict qualified as Map
import Data.Maybe (
    fromMaybe,
 )
import Dtmc.Analysis.Classification.Internal (
    backwardReachable,
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
    rowSums,
    solveIminusQVector,
    subMatrix,
 )
import Dtmc.Distribution (
    Distribution (..),
 )
import Dtmc.Dynamics.Internal (
    pushSparseWeights,
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
 )
import Dtmc.State.Internal (
    stateCardinalityInt,
    stateFromInt,
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

toState :: (FiniteState state) => Int -> state
toState index =
    fromMaybe
        (error "Dtmc.Analysis.HittingTime: graph vertex out of bounds")
        (stateFromInt index)

-- Use a mask so duplicates collapse and per-state membership stays O(1).
indexMask :: Int -> [Int] -> Unboxed.UArray Int Bool
indexMask dim indices =
    Unboxed.accumArray (||) False (0, dim - 1) [(i, True) | i <- indices]

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

{- | Compute the exact scalar hitting-time probability
@P(H_A = t | X_0 = i)@ through any 'Transition'. The target set is represented
by a membership predicate, which also works when the state space is infinite.
Hitting includes time zero, and newly hit mass is removed after every step.

Complexity: excluding 'transitionLaw' and predicate evaluation,
@O(k (w + e log(u + 1) + u) + 1)@ time, @O(w + u)@ temporary space, and
@O(1)@ result space, where @k = t@ and @w@, @e@, and @u@ bound per-step
survivor states, traversed transition edges, and accumulated destinations.
-}
exactProbabilityAt ::
    ( Transition kernel
    , Ord (TransitionState kernel)
    ) =>
    Natural ->
    kernel ->
    (TransitionState kernel -> Bool) ->
    TransitionState kernel ->
    Double
exactProbabilityAt time kernel isTarget initialState
    | time == 0 = if isTarget initialState then 1 else 0
    | isTarget initialState = 0
    | otherwise = go time (Map.singleton initialState 1)
  where
    go 0 _ = 0
    go _ survivors | Map.null survivors = 0
    go remaining survivors =
        let (next, hitMass) = advanceUntilTarget kernel isTarget survivors
         in if remaining == 1
                then hitMass
                else go (remaining - 1) next

{- | Compute the strict bounded scalar hitting probability
@P(H_A < c | X_0 = i)@ through any 'Transition'. At @c = 0@ the result is
zero; at a positive bound an initial target gives one.

Complexity: excluding 'transitionLaw' and predicate evaluation,
@O(k (w + e log(u + 1) + u) + 1)@ time, @O(w + u)@ temporary space, and
@O(1)@ result space, where @k = c@ and @w@, @e@, and @u@ bound per-step
survivor states, traversed transition edges, and accumulated destinations.
-}
lowerTailProbability ::
    ( Transition kernel
    , Ord (TransitionState kernel)
    ) =>
    Natural ->
    kernel ->
    (TransitionState kernel -> Bool) ->
    TransitionState kernel ->
    Double
lowerTailProbability bound kernel isTarget initialState
    | bound == 0 = 0
    | isTarget initialState = 1
    | otherwise = go (bound - 1) (Map.singleton initialState 1) 0
  where
    go 0 _ total = total
    go _ survivors total | Map.null survivors = total
    go remaining survivors total =
        let (next, hitMass) = advanceUntilTarget kernel isTarget survivors
            cumulative = total + hitMass
         in cumulative `seq` go (remaining - 1) next cumulative

{- | Compute hitting probabilities
@h_i = P(H_A < infinity | X_0 = i)@ in state order. Target order and
duplicates are ignored; an empty target set gives an all-zero vector.

The result is the minimal non-negative solution of @h_i = 1@ on @A@ and
@h_i = sum_j P(i,j) h_j@ elsewhere. Target entries are exactly @1@, and
states from which @A@ is unreachable are exactly @0@. Remaining entries solve
@(I - P[D,D])x = P[D,A]1@ and inherit floating-point error.

Returns 'Left' if the interior solve fails the numerical contract.

Complexity: @O(n^3 + a)@ worst-case time, @O(n^2 + a)@ temporary space,
@O(n + E)@ retained graph-cache space, and @O(n)@ result space for @n@
states, @E@ support edges, and @a@ supplied targets.
-}
eventualProbabilitiesByState ::
    (FiniteState state) =>
    TransitionMatrix state ->
    [state] ->
    Either LinearSystemError (S.R (Cardinality state))
eventualProbabilitiesByState p targets =
    -- The ordinary hitting problem is the competing problem with no competing
    -- boundary (@H_B = infinity@), so it reuses the same single solve.
    raceProbabilitiesByState p targets []

{- | Compute the probability of ever hitting the target set from one state.
This has the same edge cases, numerical behaviour, and errors as
@eventualProbabilitiesByState@.

Partially applying the matrix and target set shares one lazy all-state solve:
the first forced query computes the table, and later lookups read it directly.

Complexity: the first forced query takes @O(n^3 + a)@ worst-case time and
@O(n^2 + a)@ temporary space for @n@ states and @a@ supplied targets. It may
retain an @O(n)@ all-state result and @O(n + E)@ graph cache; subsequent
shared lookups take @O(1)@ time and space. The scalar result occupies
@O(1)@ space.
-}
eventualProbabilityGivenInitialState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    [state] ->
    state ->
    Either LinearSystemError Double
eventualProbabilityGivenInitialState p targets =
    \i -> (`LA.atIndex` toIndex i) <$> probabilities
  where
    probabilities = S.extract <$> eventualProbabilitiesByState p targets

{- | Compute competing hitting probabilities
@h_i = P(H_A < H_B | X_0 = i)@ in state order, for a successful boundary @A@
(first argument) and a competing boundary @B@ (second argument). Hitting times
include time zero, @H_A = inf { t >= 0 | X_t in A }@ and likewise for @B@, and
the comparison is /strict/: @A@ must be reached strictly before @B@.

Target order and duplicate states are ignored in both lists.

Overlap and ties. The two boundaries need not be disjoint. Reaching a state in
both @A@ and @B@ ties the two hitting times (@H_A = H_B@), and a tie fails the
strict inequality, so the competing boundary claims every shared state. Writing
the effective successful set as @A' = A \\ B@:

* states in @A'@ have value exactly @1@;
* states in @B@ -- including states shared with @A@ -- have value exactly @0@;
* if @A@ and @B@ are identical, every result is @0@.

Empty boundaries. Because @H_(empty) = infinity@:

* an empty successful set gives an all-zero vector;
* an empty competing set agrees with 'eventualProbabilitiesByState' on the same
  successful set;
* two empty sets give an all-zero vector.

Reachability is structural. Using the strict-positive support graph, one reverse
traversal from @A'@ that is forbidden to pass through @B@ marks the states that
can reach @A'@ before @B@. This never consults floating-point probabilities, so
a state that cannot reach @A'@ at all, or can reach it only by first entering
@B@, is assigned exactly @0@ without a solve. For the remaining interior states
@D@ the values are the minimal non-negative solution of
@h_i = sum_j P(i,j) h_j@, i.e. @(I - P[D,D]) x = P[D,A'] 1@; these entries come
from the shared 'Double' linear solver, inherit its rounding, and are not
clamped to @[0, 1]@, renormalised, or given any tolerance.

Returns 'Left' if the interior solve fails the numerical contract. For a valid
transition matrix the system is nonsingular in exact arithmetic, but may still
be too ill-conditioned for a reliable 'Double' result. The all-state result
performs at most one linear solve.

Complexity: @O(n^3 + a + b)@ worst-case time and @O(n^2 + a + b)@
temporary space for @n@ states and boundary-list lengths @a@ and @b@. The
matrix retains @O(n + E)@ graph-cache space, and the result occupies @O(n)@
space.
-}
raceProbabilitiesByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    [state] ->
    [state] ->
    Either LinearSystemError (S.R (Cardinality state))
raceProbabilitiesByState p successful competing = do
    solved <-
        if null interiorIdx
            then Right []
            else
                LA.toList
                    <$> solveIminusQVector
                        (subMatrix interiorIdx interiorIdx matrix)
                        (rowSums (subMatrix interiorIdx effectiveIdx matrix))
    let interiorValues :: Unboxed.UArray Int Double
        interiorValues =
            Unboxed.accumArray
                (\_ x -> x)
                0
                (0, dim - 1)
                (zip interiorIdx solved)
        valueAt i
            | inEffective i = 1
            | canReach i = interiorValues Unboxed.! i
            | otherwise = 0
    pure (S.vector [valueAt i | i <- [0 .. dim - 1]])
  where
    dim = stateCardinalityInt @state
    -- Masks keep boundary and solution lookup constant-time during assembly.
    competingMask = indexMask dim (map toIndex competing)
    inCompeting i = competingMask Unboxed.! i
    successfulMask = indexMask dim (map toIndex successful)
    -- Effective successful set A' = A \ B. A state on both boundaries is a
    -- tie, and a tie loses, so the competing boundary claims it.
    inEffective i = successfulMask Unboxed.! i && not (inCompeting i)
    effectiveIdx = [i | i <- [0 .. dim - 1], inEffective i]
    -- One reverse traversal to A', forbidden to cross B, marks every state
    -- that can reach A' before B; this is structural, not numerical.
    reachMask =
        indexMask
            dim
            ( map
                toIndex
                ( backwardReachable
                    p
                    (not . inCompeting . toIndex)
                    (map toState effectiveIdx)
                )
            )
    canReach i = reachMask Unboxed.! i
    interiorIdx = [i | i <- [0 .. dim - 1], not (inEffective i), canReach i]
    matrix = S.extract (unTransitionMatrix p)

{- | Compute the probability of hitting the successful boundary strictly
before the competing boundary from one state, @P(H_A < H_B | X_0 = i)@. This
has the same overlap, empty-set, structural, numerical, and error behaviour as
@raceProbabilitiesByState@.

Partially applying the matrix and both boundaries shares one lazy all-state
solve.

Complexity: the first forced query takes @O(n^3 + a + b)@ worst-case time and
@O(n^2 + a + b)@ temporary space for @n@ states and boundary-list lengths
@a@ and @b@. It may retain an @O(n)@ all-state result and @O(n + E)@ graph
cache; subsequent shared lookups take @O(1)@ time and space. The scalar result
occupies @O(1)@ space.
-}
raceProbabilityGivenInitialState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    [state] ->
    [state] ->
    state ->
    Either LinearSystemError Double
raceProbabilityGivenInitialState p successful competing =
    \i -> (`LA.atIndex` toIndex i) <$> probabilities
  where
    probabilities =
        S.extract <$> raceProbabilitiesByState p successful competing

{- | Compute expected hitting times @E(H_A | X_0 = i)@ in state order.
Targets have exact expectation zero. A non-target state has
'InfiniteExpectation' exactly when the target is not hit with probability
one; this is decided from support reachability, not a floating-point
comparison. An empty target set therefore gives 'InfiniteExpectation' for
every state.

Finite entries are the solution of
@eta_i = 1 + sum_(j not in A) P(i,j) eta_j@. They inherit solver rounding and
are not clamped. Returns 'Left' if the finite-state system fails the numerical
contract.

Complexity: @O(n^3 + a)@ worst-case time, @O(n^2 + a)@ temporary space,
@O(n + E)@ retained graph-cache space, and @O(n)@ result space for @n@
states, @E@ support edges, and @a@ supplied targets.
-}
expectationsByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    [state] ->
    Either LinearSystemError [Expectation]
expectationsByState p targets = do
    solved <-
        if null certainIdx
            then Right []
            else
                LA.toList
                    <$> solveIminusQVector
                        (subMatrix certainIdx certainIdx matrix)
                        (LA.konst 1 (length certainIdx))
    let certainValues :: Unboxed.UArray Int Double
        certainValues =
            Unboxed.accumArray
                (\_ x -> x)
                0
                (0, dim - 1)
                (zip certainIdx solved)
        valueAt i
            | inTarget i = FiniteExpectation 0
            | doomedMask Unboxed.! i = InfiniteExpectation
            | otherwise = FiniteExpectation (certainValues Unboxed.! i)
    pure [valueAt i | i <- [0 .. dim - 1]]
  where
    dim = stateCardinalityInt @state
    targetMask = indexMask dim (map toIndex targets)
    inTarget i = targetMask Unboxed.! i
    -- One reverse traversal replaces a reachability query for every state.
    reachMask =
        indexMask dim (map toIndex (backwardReachable p (const True) targets))
    unreachable =
        [i | i <- [0 .. dim - 1], not (inTarget i), not (reachMask Unboxed.! i)]
    -- Reaching an unreachable state without crossing the target makes the
    -- hitting probability less than one.
    doomed =
        backwardReachable
            p
            (not . inTarget . toIndex)
            (map toState unreachable)
    doomedMask = indexMask dim (map toIndex doomed)
    certainIdx =
        [i | i <- [0 .. dim - 1], not (inTarget i), not (doomedMask Unboxed.! i)]
    matrix = S.extract (unTransitionMatrix p)

{- | Compute the expected time to hit the target set from one state. This has
the same edge cases, numerical behaviour, and errors as @expectationsByState@.

Partial application shares one lazy all-state table: the first forced query
computes the table, and later lookups read it directly.

Complexity: the first forced query takes @O(n^3 + a)@ worst-case time and
@O(n^2 + a)@ temporary space for @n@ states and @a@ supplied targets. It may
retain an @O(n)@ all-state result and @O(n + E)@ graph cache; subsequent
shared lookups take @O(1)@ time and space. The scalar result occupies
@O(1)@ space.
-}
expectationGivenInitialState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    [state] ->
    state ->
    Either LinearSystemError Expectation
expectationGivenInitialState p targets =
    \i -> (Array.! toIndex i) <$> table
  where
    -- Back the shared table with a boxed array so each state query is O(1)
    -- List indexing was linear; the array keeps shared queries constant-time.
    table =
        Array.listArray (0, dim - 1) <$> expectationsByState p targets
    dim = stateCardinalityInt @state

-- Direct survivor mass @P(H_A > t)@ through a locally finite transition.
-- Newly hit paths are removed at every step, so paths that never hit remain
-- in the map and the result naturally includes the infinity atom.
upperTailProbability ::
    ( Transition kernel
    , Ord (TransitionState kernel)
    ) =>
    Natural ->
    kernel ->
    (TransitionState kernel -> Bool) ->
    TransitionState kernel ->
    Double
upperTailProbability time kernel isTarget initialState
    | isTarget initialState = 0
    | otherwise = go time (Map.singleton initialState 1)
  where
    go 0 survivors = sum (Map.elems survivors)
    go _ survivors | Map.null survivors = 0
    go remaining survivors =
        let (next, _) = advanceUntilTarget kernel isTarget survivors
         in next `seq` go (remaining - 1) next

{- | Compute, under an arbitrary initial distribution, the probability of a
finite-threshold event in the hitting time
@H_A = inf { t >= 0 | X_t in A }@. The initial distribution supplies the
probability measure @P_mu@.

'EqualTo' and the two lower tails reuse the direct exact/bounded recurrences.
'GreaterThan' and 'AtLeast' use surviving mass directly rather than subtracting
a cumulative probability from one. Consequently upper tails include the atom
at infinity and avoid cancellation when the surviving probability is small.
'AtLeast' @0@ is exactly one.

This function works through any locally finite 'Transition'. Results use
ordinary 'Double' arithmetic without clamping or renormalisation.

For the complexity bounds, @s@ is the initial stored support size, @k@ the
event threshold, and @w@, @e@, and @u@ are per-step upper bounds on survivor
states, traversed transition edges, and accumulated destinations.

Complexity: excluding 'distributionWeights', 'transitionLaw', and predicate
evaluation, @O(s (k (w + e log(u + 1) + u) + 1))@ time,
@O(s + w + u)@ temporary space, and @O(1)@ result space.
-}
probability ::
    ( Distribution distribution
    , Transition kernel
    , DistributionState distribution ~ TransitionState kernel
    , Ord (TransitionState kernel)
    ) =>
    DiscreteEvent ->
    kernel ->
    (TransitionState kernel -> Bool) ->
    distribution ->
    Double
probability event kernel isTarget initial =
    probabilityUnder initial (probabilityGivenInitialState event kernel isTarget)

{- | Compute the probability of a finite-threshold hitting-time event
conditioned on @X_0 = i@. Hitting includes time zero and follows the event
boundary behaviour documented by 'probability'.

For the complexity bounds, @k@ is the event threshold and @w@, @e@, and @u@
are per-step upper bounds on survivor states, traversed transition edges, and
accumulated destinations.

Complexity: excluding 'transitionLaw' and predicate evaluation,
@O(k (w + e log(u + 1) + u) + 1)@ time, @O(w + u)@ temporary space, and
@O(1)@ result space.
-}
probabilityGivenInitialState ::
    ( Transition kernel
    , Ord (TransitionState kernel)
    ) =>
    DiscreteEvent ->
    kernel ->
    (TransitionState kernel -> Bool) ->
    TransitionState kernel ->
    Double
probabilityGivenInitialState event kernel isTarget initialState =
    case event of
        EqualTo time ->
            exactProbabilityAt time kernel isTarget initialState
        LessThan bound ->
            lowerTailProbability bound kernel isTarget initialState
        AtMost time ->
            lowerTailProbability (time + 1) kernel isTarget initialState
        GreaterThan time ->
            upperTailProbability time kernel isTarget initialState
        AtLeast 0 -> 1
        AtLeast time ->
            upperTailProbability (time - 1) kernel isTarget initialState

{- | Compute, under an arbitrary initial distribution, the probability of ever
hitting the target set. Target order and duplicates do not affect the result;
an empty target set gives zero. Numerical failures come from the checked
all-state linear solve.

Complexity: excluding 'distributionWeights', @O(n^3 + a + s)@ worst-case
time, @O(n^2 + a + s)@ temporary space, and @O(1)@ result space for @n@
states, @a@ supplied targets, and initial stored support size @s@. The matrix
may retain @O(n + E)@ graph-cache space.
-}
eventualProbability ::
    ( FiniteState state
    , Distribution distribution
    , DistributionState distribution ~ state
    ) =>
    TransitionMatrix state ->
    [state] ->
    distribution ->
    Either LinearSystemError Double
eventualProbability matrix targets initial =
    probabilityUnderEither initial (eventualProbabilityGivenInitialState matrix targets)

{- | Compute, under an arbitrary initial distribution, the probability of
hitting the successful boundary strictly before the competing boundary.
Overlap ties lose, and empty-boundary behaviour matches
'raceProbabilityGivenInitialState'.

Complexity: excluding 'distributionWeights', @O(n^3 + a + b + s)@
worst-case time, @O(n^2 + a + b + s)@ temporary space, and @O(1)@ result
space for @n@ states, boundary-list lengths @a@ and @b@, and initial stored
support size @s@. The matrix may retain @O(n + E)@ graph-cache space.
-}
raceProbability ::
    ( FiniteState state
    , Distribution distribution
    , DistributionState distribution ~ state
    ) =>
    TransitionMatrix state ->
    [state] ->
    [state] ->
    distribution ->
    Either LinearSystemError Double
raceProbability matrix successful competing initial =
    probabilityUnderEither initial (raceProbabilityGivenInitialState matrix successful competing)

{- | Compute the expected hitting time under an arbitrary initial
distribution. It is 'InfiniteExpectation' exactly when a state of positive
initial probability does not hit the target almost surely. Numerical failures
come from the checked all-state linear solve.

Complexity: excluding 'distributionWeights', @O(n^3 + a + s)@ worst-case
time, @O(n^2 + a + s)@ temporary space, and @O(1)@ result space for @n@
states, @a@ supplied targets, and initial stored support size @s@. The matrix
may retain @O(n + E)@ graph-cache space.
-}
expectation ::
    ( FiniteState state
    , Distribution distribution
    , DistributionState distribution ~ state
    ) =>
    TransitionMatrix state ->
    [state] ->
    distribution ->
    Either LinearSystemError Expectation
expectation matrix targets initial =
    expectationUnderEither initial (expectationGivenInitialState matrix targets)
