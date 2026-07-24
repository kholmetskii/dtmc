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

{- | An expected number of transitions, represented either by a 'Double' or an
exact infinite case. Library results use 'InfiniteMean' based on support-graph
reachability rather than floating-point overflow.

'FiniteMean' performs no validation: callers can construct negative,
non-finite, or @NaN@ values. Derived ordering places every 'FiniteMean'
constructor before 'InfiniteMean'; comparisons between finite constructors
inherit the behavior of 'Double', including @NaN@.
-}
data MeanTime
    = -- | A mathematically non-negative finite mean, subject to solver rounding.
      FiniteMean Double
    | -- | The target or return is not reached with probability one.
      InfiniteMean
    deriving (Eq, Ord, Show)

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
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    [Finite n] ->
    S.R n
hittingProbabilities p targets =
    S.vector [valueAt i | i <- [0 .. dim - 1]]
  where
    dim = fromIntegral (natVal (Proxy @n))
    -- Masks keep target and solution lookup constant-time during assembly.
    targetMask = indexMask dim (map toIndex targets)
    inTarget i = targetMask Unboxed.! i
    targetIdx = [i | i <- [0 .. dim - 1], inTarget i]
    -- One reverse traversal avoids a reachability search for every state.
    reachMask =
        indexMask dim (map toIndex (backwardReachable p (const True) targets))
    canReach i = reachMask Unboxed.! i
    interiorIdx = [i | i <- [0 .. dim - 1], not (inTarget i), canReach i]
    matrix = S.extract (unTransitionMatrix p)
    solved
        | null interiorIdx = []
        | otherwise =
            case solveIminusQVector
                (subMatrix interiorIdx interiorIdx matrix)
                (rowSums (subMatrix interiorIdx targetIdx matrix)) of
                Just x -> LA.toList x
                Nothing ->
                    error
                        "Dtmc.Hitting.hittingProbabilities: interior system \
                        \singular; impossible for a valid transition matrix"
    interiorValues :: Unboxed.UArray Int Double
    interiorValues =
        Unboxed.accumArray (\_ x -> x) 0 (0, dim - 1) (zip interiorIdx solved)
    valueAt i
        | inTarget i = 1
        | canReach i = interiorValues Unboxed.! i
        | otherwise = 0

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
