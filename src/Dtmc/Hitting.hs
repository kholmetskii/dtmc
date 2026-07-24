{- |
Module      : Dtmc.Hitting
Description : Hitting probabilities, expected hitting times, return times.

The quantitative hitting theory of a finite chain: for a target set
@A@ of states, the probability @h_iA@ of ever reaching @A@ from @i@, the
expected number of steps @eta_iA@ to do so, and the expected time to
/return/ to a state. Everything rests on first-step decompositions:
conditioning on the first transition turns each
quantity into a linear system over the states outside the target, and for a
finite chain the solution sets can be pinned down exactly.

The module keeps the exact and the approximate strictly apart. Which states
have @h_iA = 1@, @h_iA = 0@, or @eta_iA@ infinite is decided by /exact/
combinatorics on the support graph (reachability, never by comparing a
float to one). Only the genuinely fractional values -- probabilities
strictly between @0@ and @1@, finite means -- come from floating-point
solves of @(I - Q) x = b@ in "Dtmc.Internal.LinearSystem", and inherit its
rounding. Infinity is likewise a constructor ('InfiniteMean'), not an IEEE
value.
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

{- | An expected number of steps: a finite mean, or provably infinite. The
infinite case is a constructor rather than an IEEE infinity so that "the
mean does not exist" is an exact statement (decided combinatorially), never
the outcome of float arithmetic. The derived 'Ord' agrees with the order of
the extended reals: finite means compare by value and every finite mean is
below 'InfiniteMean'.
-}
data MeanTime
    = FiniteMean Double
    | InfiniteMean
    deriving (Eq, Ord, Show)

-- Convert a bounded 'Finite' @n@ index back to a raw @Int@.
toIndex :: Finite n -> Int
toIndex = fromIntegral . getFinite

-- Convert a raw @Int@ index into the bounded 'Finite' @n@ index.
toFinite :: (KnownNat n) => Int -> Finite n
toFinite = finite . fromIntegral

-- A Boolean membership mask over states @0 .. dim-1@: 'True' at exactly the
-- listed indices (duplicates collapse). Backs O(1) membership tests, keeping
-- the surrounding per-state assembly linear instead of quadratic list @elem@.
indexMask :: Int -> [Int] -> Unboxed.UArray Int Bool
indexMask dim indices =
    Unboxed.accumArray (||) False (0, dim - 1) [(i, True) | i <- indices]

{- | The vector of hitting probabilities @h_iA = P(H_A < infinity | X_0 = i)@,
one entry per state, for the target set @A@ (order and duplicates in the
target list are irrelevant; an empty target is never hit, so @h = 0@).

/Theorem (first-step decomposition)./ @(h_iA)@ is the minimal
non-negative solution of @h_iA = 1@ on @A@ and
@h_iA = sum_j P_ij h_jA@ off @A@.

The system as stated can have many solutions (any state that cannot reach
@A@ satisfies its equation with /any/ constant), so the implementation
first splits the states exactly: @h = 1@ on @A@; @h = 0@ on the states
from which @A@ is unreachable in the support graph (the minimal choice);
and only the remaining interior @D@ -- states off @A@ that can reach @A@ --
is solved as @(I - P_DD) x = P_DA 1@.

/Proof that the interior solve is exact./ On @D@ the restricted system has
a unique solution, which is then forced to be the minimal one: every
@i@ in @D@ reaches @A@ by some support path of length at most @|D|@, each
step of which has positive probability, so the probability of remaining
inside @D@ for @|D|@ consecutive steps is uniformly below one. Hence
@||(P_DD)^{|D|}|| < 1@ in the row-sum norm, the spectral radius of @P_DD@
is below one, and @I - P_DD@ is invertible.

Entries on @A@ are exactly @1@ and entries on the unreachable set exactly
@0@; only the interior entries carry floating-point error from the solve.
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
    -- Target membership as an O(1) mask (duplicates collapse).
    targetMask = indexMask dim (map toIndex targets)
    inTarget i = targetMask Unboxed.! i
    -- Unique target indices, ascending, for the solve's column block.
    targetIdx = [i | i <- [0 .. dim - 1], inTarget i]
    -- One reverse traversal replaces the per-state reachability search:
    -- @reachMask@ marks every state that can reach the target.
    reachMask =
        indexMask dim (map toIndex (backwardReachable p (const True) targets))
    canReach i = reachMask Unboxed.! i
    -- Interior D: off the target and able to reach it.
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
    -- Solved interior values indexed by state for O(1) lookup.
    interiorValues :: Unboxed.UArray Int Double
    interiorValues =
        Unboxed.accumArray (\_ x -> x) 0 (0, dim - 1) (zip interiorIdx solved)
    valueAt i
        | inTarget i = 1
        | canReach i = interiorValues Unboxed.! i
        | otherwise = 0

{- | The probability of ever hitting the target set from one supplied state.
This is an indexed view of 'hittingProbabilities'; partial application to a
matrix and target set shares the one global solve across subsequent state
queries.
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

{- | The expected hitting times @eta_iA = E(H_A | X_0 = i)@ of the target set
@A@, one entry per state in state order.

/Theorem./ @eta_iA = 0@ on @A@; @eta_iA@ is infinite exactly
when @h_iA < 1@; and on the states with @h_iA = 1@ the family is the
minimal non-negative solution of
@eta_iA = 1 + sum_{j not in A} P_ij eta_jA@.

Which states have @h_iA < 1@ is decided exactly, not by comparing the
solved @h@ to one: let @Z@ be the states that cannot reach @A@ at all.
Then @h_iA < 1@ iff @i@ can reach @Z@ by a support path avoiding @A@.

/Proof./ If such a path exists it has positive probability, and after it
the chain can never hit @A@, so @h_iA < 1@. Conversely, on the event that
@A@ is never hit the trajectory stays in the complement of @A@ forever,
and (finitely many states) visits some state @j@ infinitely often; were
@A@ reachable from @j@, each visit would give a probability bounded away
from zero of hitting @A@ within @n@ steps, forcing a hit almost surely.
So @j@ lies in @Z@, reached by a path avoiding @A@. Hence if no such path
exists, the never-hit event has probability zero and @h_iA = 1@.

The certain states @B@ (off @A@, not doomed) step only into @A@ or @B@ --
an edge from certain @i@ to doomed @j@ would make @i@ doomed -- so the
system @(I - P_BB) x = 1@ is closed, and @I - P_BB@ is invertible by the
same spectral-radius argument as in 'hittingProbabilities' (from every
state of @B@ the chain exits to @A@ with uniform positive probability).
Minimality again follows from uniqueness on @B@.
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
    -- One reverse traversal: states that can reach the target. Off the target,
    -- its complement is Z -- the states that cannot reach the target at all.
    reachMask =
        indexMask dim (map toIndex (backwardReachable p (const True) targets))
    unreachable =
        [i | i <- [0 .. dim - 1], not (inTarget i), not (reachMask Unboxed.! i)]
    -- Doomed: off-target states that can reach Z along a target-avoiding path
    -- (reverse reachability of Z within the off-target subgraph, O(V + E)). On
    -- these h < 1, so eta is infinite. Held as an O(1) membership mask.
    doomed =
        backwardReachable
            p
            (not . inTarget . toIndex)
            (map toFinite unreachable)
    doomedMask = indexMask dim (map toIndex doomed)
    -- Certain B: off the target and not doomed.
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
    -- Solved certain means indexed by state for O(1) lookup. Every off-target,
    -- non-doomed state is certain, so the partition is total and no lookup can
    -- miss.
    certainValues :: Unboxed.UArray Int Double
    certainValues =
        Unboxed.accumArray (\_ x -> x) 0 (0, dim - 1) (zip certainIdx solved)
    valueAt i
        | inTarget i = FiniteMean 0
        | doomedMask Unboxed.! i = InfiniteMean
        | otherwise = FiniteMean (certainValues Unboxed.! i)

{- | The expected time to hit the target set from one supplied state. This is
an indexed view of 'expectedHittingTimes'; partial application shares the
table and its linear solve across subsequent state queries.
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

{- | The vector of first-return probabilities
@f_i = P(T_i < infinity | X_0 = i)@, one entry per state, where @T_i@ is
the first return time to @i@ after at least one step.

Recurrent states are assigned exactly @1@ from the support-graph
classification. For all transient states at once, let @Q@ be the transient
block and @N = (I - Q)^-1@ its fundamental matrix. Since @N_ii@ is the
expected number of visits to @i@ starting from @i@, including time zero,
the renewal identity @N_ii = 1 + f_i N_ii@ gives
@f_i = 1 - 1 / N_ii@. Thus the whole vector needs one matrix solve rather
than one singleton hitting solve per state.
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
    -- Transient return probabilities f_i = 1 - 1/N_ii, all from one solve.
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
    -- Indexed by state for O(1) lookup, mirroring 'hittingProbabilities';
    -- recurrent entries default to 0 and are never read.
    transientValues :: Unboxed.UArray Int Double
    transientValues =
        Unboxed.accumArray (\_ x -> x) 0 (0, dim - 1) (zip transientIdx transientReturns)
    valueAt i
        | recurrentState p i = 1
        | otherwise = transientValues Unboxed.! toIndex i

{- | The probability of returning to one supplied state after at least one
step. This is an indexed view of 'returnProbabilities'; a recurrent state
returns @1@ immediately from the support-graph classification, so a query
confined to recurrent states never forces the fundamental-matrix solve.
Partial application shares that single solve across transient-state queries.
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

{- | The expected first-return times @m_i = E(T_i | X_0 = i)@, one entry per
state. The implementation uses only the already established first-step
decomposition: after the first move @i -> j@, returning to @i@ is the same
as hitting the singleton target @{i}@ from @j@. Thus

@m_i = 1 + sum_j P_ij eta_j{i}@,

where @eta_i{i} = 0@ handles an immediate self-loop. A positive-probability
successor with infinite singleton hitting mean makes @m_i@ infinite; zero
probability times infinity contributes zero.

For a finite chain a state is transient exactly when its mean return time is
infinite, so transient states are read off as 'InfiniteMean' directly from
the support-graph classification and skip the solve; only recurrent states
run the @O(n^3)@ hitting-time solve. The worst case is still @O(n^4)@ -- an
irreducible chain, where every state is recurrent -- but chains with
transient states are markedly cheaper.
-}
expectedReturnTimes ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    [MeanTime]
expectedReturnTimes p =
    map (expectedReturnTime p) finites

{- | The expected first-return time for one supplied state. A transient state
has infinite mean return time by the finite-chain classification, returned
immediately without any solve; a recurrent state performs only the one
required singleton hitting solve, @O(n^3)@ worst case, rather than computing
return means for every possible starting state.
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

{-
More efficient all-states alternative, intentionally kept out of the active
implementation until stationary distributions and Kac's return-time theorem
have been introduced. It computes every expected return time in O(n^3): solve
one stationary system per closed recurrent class, assign 1 / pi_i to its
members, and assign InfiniteMean to transient states. To activate this code,
also import communicatingClasses from Dtmc.Classification.

stationaryDistribution :: LA.Matrix Double -> Maybe (LA.Vector Double)
stationaryDistribution q
    | LA.rows q == 0 = Just (LA.fromList [])
    | otherwise =
        LA.flatten <$> LA.linearSolve coefficients (LA.asColumn rhs)
  where
    size = LA.rows q
    balanceRows = take (size - 1) (LA.toRows (LA.tr q - LA.ident size))
    coefficients = LA.fromRows (balanceRows ++ [LA.konst 1 size])
    rhs = LA.fromList (replicate (size - 1) 0 ++ [1])

expectedReturnTimesViaStationary ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    [MeanTime]
expectedReturnTimesViaStationary p =
    [valueAt i | i <- finites]
  where
    matrix = S.extract (unTransitionMatrix p)
    recurrentClasses =
        [ members
        | members@(representative : _) <- communicatingClasses p
        , recurrentState p representative
        ]
    recurrentValue = concatMap solveClass recurrentClasses
    solveClass members =
        case stationaryDistribution (subMatrix indices indices matrix) of
            Just piVector -> zip members (map meanFromMass (LA.toList piVector))
            Nothing ->
                error
                    "Dtmc.Hitting.expectedReturnTimesViaStationary: stationary system \
                    \singular or numerically ill-conditioned"
      where
        indices = map toIndex members
    meanFromMass mass
        | mass > 0 = FiniteMean (1 / mass)
        | otherwise =
            error
                "Dtmc.Hitting.expectedReturnTimesViaStationary: \
                \non-positive stationary mass"
    valueAt i = fromMaybe InfiniteMean (lookup i recurrentValue)
-}
