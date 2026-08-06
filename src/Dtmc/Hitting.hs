{- |
Module      : Dtmc.Hitting
Description : Hitting probabilities, expected hitting times, return times.

Hitting and first-return quantities for a finite DTMC. For a target set @A@,
@H_A = inf { t >= 0 | X_t in A }@; for state @i@,
@T_i = inf { t >= 1 | X_t = i }@.

Reachability decides unreachable zeros, recurrent returns, and infinite means
from the strict-positive support graph. Other values, including non-target
hitting probabilities that are mathematically one, use 'Double' LU solves of
@(I - Q)x = b@. No tolerance, conditioning test, residual check, or clamping
is applied, so results may leave their mathematical ranges or become
non-finite. A backend-reported singular system raises an error; for a valid
transition matrix the systems are nonsingular in exact arithmetic.

All-state functions return empty results for a @0 x 0@ matrix. State-specific
functions cannot be called because no 'Finite 0' state exists.
-}
module Dtmc.Hitting (
    MeanTime (..),
    hittingProbabilities,
    hittingProbability,
    hittingBeforeProbabilities,
    hittingBeforeProbability,
    expectedHittingTimes,
    expectedHittingTime,
    returnProbabilities,
    returnProbability,
    expectedReturnTimes,
    expectedReturnTime,
) where

import Data.Array qualified as Array
import Data.Array.Unboxed qualified as Unboxed
import Data.Finite (
    Finite,
    finite,
    finites,
    getFinite,
 )
import Data.Proxy (
    Proxy (..),
 )
import Dtmc.Classification (
    backwardReachable,
    recurrentState,
    transientStates,
 )
import Dtmc.Distribution.Internal (
    unDistribution,
 )
import Dtmc.Hitting.Internal (
    MeanTime (..),
 )
import Dtmc.Internal.LinearSystem (
    fundamental,
    rowSums,
    solveIminusQVector,
    subMatrix,
 )
import Dtmc.TransitionMatrix (
    rowAt,
 )
import Dtmc.TransitionMatrix.Internal (
    TransitionMatrix,
    unTransitionMatrix,
 )
import GHC.TypeNats (
    KnownNat,
    natVal,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

toIndex :: Finite n -> Int
toIndex = fromIntegral . getFinite

toFinite :: (KnownNat n) => Int -> Finite n
toFinite = finite . fromIntegral

-- Use a mask so duplicates collapse and per-state membership stays O(1).
indexMask :: Int -> [Int] -> Unboxed.UArray Int Bool
indexMask dim indices =
    Unboxed.accumArray (||) False (0, dim - 1) [(i, True) | i <- indices]

{- | Hitting probabilities
@h_i = P(H_A < infinity | X_0 = i)@ in state order. Target order and
duplicates are ignored; an empty target set gives an all-zero vector.

The result is the minimal non-negative solution of @h_i = 1@ on @A@ and
@h_i = sum_j P(i,j) h_j@ elsewhere. Target entries are exactly @1@, and
states from which @A@ is unreachable are exactly @0@. Remaining entries solve
@(I - P[D,D])x = P[D,A]1@ and inherit floating-point error.

Raises an error if the backend rejects the interior solve. Worst-case time:
@O(n^3)@; temporary space: @O(n^2)@; result space: @O(n)@.
-}
hittingProbabilities ::
    (KnownNat n) =>
    TransitionMatrix n ->
    [Finite n] ->
    S.R n
hittingProbabilities p targets =
    -- The ordinary hitting problem is the competing problem with no competing
    -- boundary (@H_B = infinity@), so it reuses the same single solve.
    hittingBeforeProbabilities p targets []

{- | The probability of ever hitting the target set from one state. This has
the same edge cases, numerical behavior, and errors as 'hittingProbabilities'.

Partially applying the matrix and target set shares one lazy all-state solve:
the first forced query costs @O(n^3)@ worst case and later lookups cost
@O(1)@.
-}
hittingProbability ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    [Finite n] ->
    Finite n ->
    Double
hittingProbability p targets =
    \i -> probabilities `LA.atIndex` toIndex i
  where
    probabilities = S.extract (hittingProbabilities p targets)

{- | Competing hitting probabilities
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
* an empty competing set agrees with 'hittingProbabilities' on the same
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

Raises an error if the backend rejects the interior solve; for a valid
transition matrix that system is nonsingular in exact arithmetic. The all-state
result performs at most one linear solve.

Worst-case time: @O(n^3)@; temporary space: @O(n^2)@; result space: @O(n)@.
-}
hittingBeforeProbabilities ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    [Finite n] ->
    [Finite n] ->
    S.R n
hittingBeforeProbabilities p successful competing =
    S.vector [valueAt i | i <- [0 .. dim - 1]]
  where
    dim = fromIntegral (natVal (Proxy @n))
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
                    (map toFinite effectiveIdx)
                )
            )
    canReach i = reachMask Unboxed.! i
    interiorIdx = [i | i <- [0 .. dim - 1], not (inEffective i), canReach i]
    matrix = S.extract (unTransitionMatrix p)
    solved
        | null interiorIdx = []
        | otherwise =
            case solveIminusQVector
                (subMatrix interiorIdx interiorIdx matrix)
                (rowSums (subMatrix interiorIdx effectiveIdx matrix)) of
                Just x -> LA.toList x
                Nothing ->
                    error
                        "Dtmc.Hitting.hittingBeforeProbabilities: interior system \
                        \singular; impossible for a valid transition matrix"
    interiorValues :: Unboxed.UArray Int Double
    interiorValues =
        Unboxed.accumArray (\_ x -> x) 0 (0, dim - 1) (zip interiorIdx solved)
    valueAt i
        | inEffective i = 1
        | canReach i = interiorValues Unboxed.! i
        | otherwise = 0

{- | The probability of hitting the successful boundary strictly before the
competing boundary from one state, @P(H_A < H_B | X_0 = i)@. This has the same
overlap, empty-set, structural, numerical, and error behaviour as
'hittingBeforeProbabilities'.

Partially applying the matrix and both boundaries shares one lazy all-state
solve: the first forced query costs @O(n^3)@ worst case and later lookups cost
@O(1)@.
-}
hittingBeforeProbability ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    [Finite n] ->
    [Finite n] ->
    Finite n ->
    Double
hittingBeforeProbability p successful competing =
    \i -> probabilities `LA.atIndex` toIndex i
  where
    probabilities = S.extract (hittingBeforeProbabilities p successful competing)

{- | Expected hitting times @E(H_A | X_0 = i)@ in state order. Targets have
exact mean zero. A non-target state is 'InfiniteMean' exactly when the target
is not hit with probability one; this is decided from support reachability,
not a floating-point comparison. An empty target set therefore gives
'InfiniteMean' for every state.

Finite entries are the solution of
@eta_i = 1 + sum_(j not in A) P(i,j) eta_j@. They inherit solver rounding and
are not clamped. Raises an error if the backend rejects the finite-state
system.

Worst-case time: @O(n^3)@; temporary space: @O(n^2)@; result space: @O(n)@.
-}
expectedHittingTimes ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    [Finite n] ->
    [MeanTime]
expectedHittingTimes p targets =
    [valueAt i | i <- [0 .. dim - 1]]
  where
    dim = fromIntegral (natVal (Proxy @n))
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
            (map toFinite unreachable)
    doomedMask = indexMask dim (map toIndex doomed)
    certainIdx =
        [i | i <- [0 .. dim - 1], not (inTarget i), not (doomedMask Unboxed.! i)]
    matrix = S.extract (unTransitionMatrix p)
    solved
        | null certainIdx = []
        | otherwise =
            case solveIminusQVector
                (subMatrix certainIdx certainIdx matrix)
                (LA.konst 1 (length certainIdx)) of
                Just x -> LA.toList x
                Nothing ->
                    error
                        "Dtmc.Hitting.expectedHittingTimes: certain system \
                        \singular; impossible for a valid transition matrix"
    -- A mask-backed table avoids linear list lookup while assembling results.
    certainValues :: Unboxed.UArray Int Double
    certainValues =
        Unboxed.accumArray (\_ x -> x) 0 (0, dim - 1) (zip certainIdx solved)
    valueAt i
        | inTarget i = FiniteMean 0
        | doomedMask Unboxed.! i = InfiniteMean
        | otherwise = FiniteMean (certainValues Unboxed.! i)

{- | The expected time to hit the target set from one state. This has the same
edge cases, numerical behavior, and errors as 'expectedHittingTimes'.

Partial application shares one lazy all-state table: the first forced query
costs @O(n^3)@ worst case and later lookups cost @O(1)@.
-}
expectedHittingTime ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    [Finite n] ->
    Finite n ->
    MeanTime
expectedHittingTime p targets =
    \i -> table Array.! toIndex i
  where
    -- Back the shared table with a boxed array so each state query is O(1)
    -- (list @!!@ was O(index)); the single solve is still shared across queries.
    table = Array.listArray (0, dim - 1) (expectedHittingTimes p targets)
    dim = fromIntegral (natVal (Proxy @n))

{- | First-return probabilities
@f_i = P(T_i < infinity | X_0 = i)@ in state order. Recurrent states are
exactly @1@ from support classification.

For all transient states, one fundamental-matrix solve computes
@N = (I - Q)^-1@ and @f_i = 1 - 1/N(i,i)@. Transient results inherit solver
rounding and are not clamped to @[0,1]@. Raises an error if the backend rejects
the transient system.

Worst-case time: @O(n^3)@; temporary space: @O(n^2)@; result space: @O(n)@.
-}
returnProbabilities ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    S.R n
returnProbabilities p =
    S.vector [valueAt i | i <- finites]
  where
    dim = fromIntegral (natVal (Proxy @n))
    transient = transientStates p
    transientIdx = map toIndex transient
    matrix = S.extract (unTransitionMatrix p)
    -- One fundamental matrix supplies every transient return probability.
    transientReturns
        | null transient = []
        | otherwise =
            case fundamental (subMatrix transientIdx transientIdx matrix) of
                Just nMatrix ->
                    [ 1 - 1 / (nMatrix `LA.atIndex` (k, k))
                    | k <- [0 .. length transient - 1]
                    ]
                Nothing ->
                    error
                        "Dtmc.Hitting.returnProbabilities: transient system \
                        \singular or numerically ill-conditioned"
    -- Index results once so output assembly does not scan the transient list.
    transientValues :: Unboxed.UArray Int Double
    transientValues =
        Unboxed.accumArray (\_ x -> x) 0 (0, dim - 1) (zip transientIdx transientReturns)
    valueAt i
        | recurrentState p i = 1
        | otherwise = transientValues Unboxed.! toIndex i

{- | The probability of returning to one state after at least one transition.
A recurrent-state query returns exactly @1@ without forcing the
fundamental-matrix solve. The first transient query costs @O(n^3)@ worst case;
partial application shares that solve, making later lookups @O(1)@.

Transient queries inherit the numerical behavior and errors of
'returnProbabilities'.
-}
returnProbability ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    Finite n ->
    Double
returnProbability p =
    \i ->
        if recurrentState p i
            then 1
            else probabilities `LA.atIndex` toIndex i
  where
    probabilities = S.extract (returnProbabilities p)

{- | Expected first-return times @E(T_i | X_0 = i)@ in state order. Transient
states are exactly 'InfiniteMean'. Each recurrent state uses the first-step
identity
@m_i = 1 + sum_j P(i,j) eta_j({i})@, where @eta_j({i})@ is the expected
hitting time of the singleton target @{i}@ from @j@.

The implementation performs a separate singleton hitting-time solve for each
recurrent state. For @r@ recurrent states, time is @O(n^2 + r n^3)@, at most
@O(n^4)@; each solve uses @O(n^2)@ temporary space and the result uses
@O(n)@ space.
-}
expectedReturnTimes ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    [MeanTime]
-- Reuse the hitting-time path; an all-state Kac calculation would require
-- stationary-distribution machinery not otherwise present in this module.
expectedReturnTimes p =
    map (expectedReturnTime p) finites

{- | The expected first-return time for one state. A transient state returns
'InfiniteMean' without a linear solve. A recurrent state performs one
singleton hitting-time solve and uses only stored transition probabilities
greater than zero; zero and tolerated negative entries contribute nothing.

A recurrent query takes @O(n^3)@ time and @O(n^2)@ temporary space. After
classification is cached, a transient query takes @O(1)@. Recurrent queries
inherit the numerical behavior and errors of 'expectedHittingTimes'.
-}
expectedReturnTime ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    Finite n ->
    MeanTime
expectedReturnTime p i
    | recurrentState p i = expectedReturnTimeFrom p i
    | otherwise = InfiniteMean

expectedReturnTimeFrom ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    Finite n ->
    MeanTime
expectedReturnTimeFrom p i =
    foldl' addTerm (FiniteMean 1) (zip finites row)
  where
    eta = expectedHittingTime p [i]
    row = LA.toList (S.extract (unDistribution (rowAt p i)))
    addTerm acc (j, pij)
        | pij <= 0 = acc
        | otherwise =
            case (acc, eta j) of
                (FiniteMean total, FiniteMean hittingTime) ->
                    FiniteMean (total + pij * hittingTime)
                _ -> InfiniteMean
