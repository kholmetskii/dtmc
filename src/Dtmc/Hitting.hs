-- |
-- Module      : Dtmc.Hitting
-- Description : Hitting probabilities, expected hitting times, return times.
--
-- The quantitative hitting theory of a finite chain: for a target set
-- @A@ of states, the probability @h_iA@ of ever reaching @A@ from @i@, the
-- expected number of steps @eta_iA@ to do so, and the expected time to
-- /return/ to a state. Everything rests on first-step decompositions:
-- conditioning on the first transition turns each
-- quantity into a linear system over the states outside the target, and for a
-- finite chain the solution sets can be pinned down exactly.
--
-- The module keeps the exact and the approximate strictly apart. Which states
-- have @h_iA = 1@, @h_iA = 0@, or @eta_iA@ infinite is decided by /exact/
-- combinatorics on the support graph (reachability, never by comparing a
-- float to one). Only the genuinely fractional values -- probabilities
-- strictly between @0@ and @1@, finite means -- come from floating-point
-- solves of @(I - Q) x = b@ in "Dtmc.Internal.Block", and inherit its
-- rounding. Infinity is likewise a constructor ('InfiniteMean'), not an IEEE
-- value.
module Dtmc.Hitting (
    MeanHittingTime (..),
    hittingProbabilities,
    expectedHittingTimes,
    returnProbability,
    expectedReturnTime,
) where

import Data.Finite (
    Finite,
    finites,
    getFinite,
 )
import Data.List (
    nub,
 )
import Data.Maybe (
    fromMaybe,
 )
import Dtmc.Classification (
    accessibleIn,
    supportEdgeIn,
    supportGraphOf,
 )
import Dtmc.Internal.Types ( 
    unDistribution 
 )
import Dtmc.Internal.Block (
    rowSums,
    solveIminusQVector,
    subMatrix,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    rowAt,
    unTransitionMatrix,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

-- | An expected number of steps: a finite mean, or provably infinite. The
-- infinite case is a constructor rather than an IEEE infinity so that "the
-- mean does not exist" is an exact statement (decided combinatorially), never
-- the outcome of float arithmetic. The derived 'Ord' agrees with the order of
-- the extended reals: finite means compare by value and every finite mean is
-- below 'InfiniteMean'.
data MeanHittingTime
    = FiniteMean Double
    | InfiniteMean
    deriving (Eq, Ord, Show)

-- Convert a bounded 'Finite' @n@ index back to a raw @Int@.
toIndex :: Finite n -> Int
toIndex = fromIntegral . getFinite

-- | The vector of hitting probabilities @h_iA = P(H_A < infinity | X_0 = i)@,
-- one entry per state, for the target set @A@ (order and duplicates in the
-- target list are irrelevant; an empty target is never hit, so @h = 0@).
--
-- /Theorem (first-step decomposition)./ @(h_iA)@ is the minimal
-- non-negative solution of @h_iA = 1@ on @A@ and
-- @h_iA = sum_j P_ij h_jA@ off @A@.
--
-- The system as stated can have many solutions (any state that cannot reach
-- @A@ satisfies its equation with /any/ constant), so the implementation
-- first splits the states exactly: @h = 1@ on @A@; @h = 0@ on the states
-- from which @A@ is unreachable in the support graph (the minimal choice);
-- and only the remaining interior @D@ -- states off @A@ that can reach @A@ --
-- is solved as @(I - P_DD) x = P_DA 1@.
--
-- /Proof that the interior solve is exact./ On @D@ the restricted system has
-- a unique solution, which is then forced to be the minimal one: every
-- @i@ in @D@ reaches @A@ by some support path of length at most @|D|@, each
-- step of which has positive probability, so the probability of remaining
-- inside @D@ for @|D|@ consecutive steps is uniformly below one. Hence
-- @||(P_DD)^{|D|}|| < 1@ in the row-sum norm, the spectral radius of @P_DD@
-- is below one, and @I - P_DD@ is invertible.
--
-- Entries on @A@ are exactly @1@ and entries on the unreachable set exactly
-- @0@; only the interior entries carry floating-point error from the solve.
hittingProbabilities ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    [Finite n] ->
    S.R n
hittingProbabilities p targets =
    S.vector [valueAt i | i <- finites]
  where
    sg = supportGraphOf p
    targetSet = nub targets
    inTarget i = i `elem` targetSet
    canReach i = any (accessibleIn sg i) targetSet
    interior = [i | i <- finites, not (inTarget i), canReach i]
    matrix = S.extract (unTransitionMatrix p)
    interiorIdx = map toIndex interior
    targetIdx = map toIndex targetSet
    solved
        | null interior = []
        | otherwise =
            case solveIminusQVector
                (subMatrix interiorIdx interiorIdx matrix)
                (rowSums (subMatrix interiorIdx targetIdx matrix)) of
                Just x -> LA.toList x
                Nothing ->
                    error
                        "Dtmc.Hitting.hittingProbabilities: interior system \
                        \singular; impossible for a valid transition matrix"
    interiorValue = zip interior solved
    valueAt i
        | inTarget i = 1
        | otherwise = fromMaybe 0 (lookup i interiorValue)

-- | The expected hitting times @eta_iA = E(H_A | X_0 = i)@ of the target set
-- @A@, as a total function of the start state (backed by a table computed
-- once, so partial application shares the work).
--
-- /Theorem./ @eta_iA = 0@ on @A@; @eta_iA@ is infinite exactly
-- when @h_iA < 1@; and on the states with @h_iA = 1@ the family is the
-- minimal non-negative solution of
-- @eta_iA = 1 + sum_{j not in A} P_ij eta_jA@.
--
-- Which states have @h_iA < 1@ is decided exactly, not by comparing the
-- solved @h@ to one: let @Z@ be the states that cannot reach @A@ at all.
-- Then @h_iA < 1@ iff @i@ can reach @Z@ by a support path avoiding @A@.
--
-- /Proof./ If such a path exists it has positive probability, and after it
-- the chain can never hit @A@, so @h_iA < 1@. Conversely, on the event that
-- @A@ is never hit the trajectory stays in the complement of @A@ forever,
-- and (finitely many states) visits some state @j@ infinitely often; were
-- @A@ reachable from @j@, each visit would give a probability bounded away
-- from zero of hitting @A@ within @n@ steps, forcing a hit almost surely.
-- So @j@ lies in @Z@, reached by a path avoiding @A@. Hence if no such path
-- exists, the never-hit event has probability zero and @h_iA = 1@.
--
-- The certain states @B@ (off @A@, not doomed) step only into @A@ or @B@ --
-- an edge from certain @i@ to doomed @j@ would make @i@ doomed -- so the
-- system @(I - P_BB) x = 1@ is closed, and @I - P_BB@ is invertible by the
-- same spectral-radius argument as in 'hittingProbabilities' (from every
-- state of @B@ the chain exits to @A@ with uniform positive probability).
-- Minimality again follows from uniqueness on @B@.
expectedHittingTimes ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    [Finite n] ->
    Finite n ->
    MeanHittingTime
expectedHittingTimes p targets =
    \i -> table !! toIndex i
  where
    sg = supportGraphOf p
    targetSet = nub targets
    inTarget i = i `elem` targetSet
    offTarget = [i | i <- finites, not (inTarget i)]
    -- Z: cannot reach the target at all.
    unreachable = [i | i <- offTarget, not (any (accessibleIn sg i) targetSet)]
    -- Doomed: can reach Z by a support path avoiding the target, i.e. the
    -- backward closure of Z along edges from non-target states. On these
    -- states h < 1, hence eta is infinite.
    doomed = grow unreachable
    grow current
        | null newcomers = current
        | otherwise = grow (current ++ newcomers)
      where
        newcomers =
            [ i
            | i <- offTarget
            , i `notElem` current
            , any (supportEdgeIn sg i) current
            ]
    -- B: off the target and certain to hit it.
    certain = [i | i <- offTarget, i `notElem` doomed]
    matrix = S.extract (unTransitionMatrix p)
    certainIdx = map toIndex certain
    solved
        | null certain = []
        | otherwise =
            case solveIminusQVector
                (subMatrix certainIdx certainIdx matrix)
                (LA.konst 1 (length certainIdx)) of
                Just x -> LA.toList x
                Nothing ->
                    error
                        "Dtmc.Hitting.expectedHittingTimes: certain system \
                        \singular; impossible for a valid transition matrix"
    certainValue = zip certain solved
    table = [valueAt i | i <- finites]
    valueAt i
        | inTarget i = FiniteMean 0
        | i `elem` doomed = InfiniteMean
        | otherwise =
            maybe
                (error "Dtmc.Hitting.expectedHittingTimes: state escaped the partition")
                FiniteMean
                (lookup i certainValue)

-- | The return probability @f_i = P(T_i < infinity | X_0 = i)@, where @T_i@
-- is the time of the /first return/ to @i@ (at least one step -- unlike the
-- hitting probability of @{i}@ from @i@, which is trivially one).
--
-- /Theorem (first-step decomposition)./ @f_i = sum_j P_ij h_j{i}@.
--
-- /Proof./ Condition on the first step: the chain moves to @j@ with
-- probability @P_ij@, after which returning to @i@ is exactly hitting @{i}@
-- from @j@. The term @j = i@ is covered by the boundary value
-- @h_i{i} = 1@: a self-loop is an immediate return.
--
-- Recurrence of a state is classically defined by @f_i = 1@; for finite
-- chains this agrees
-- with the closed-class criterion of 'Dtmc.Classification.recurrentState',
-- and the spec checks the two implementations against each other. Note the
-- result is a floating-point value from the hitting solve: recurrence of a
-- state is decided exactly by 'Dtmc.Classification.recurrentState', not by
-- comparing this number to one.
returnProbability ::
    (KnownNat n) =>
    TransitionMatrix n ->
    Finite n ->
    Double
returnProbability p i =
    LA.dot
        (S.extract (unDistribution (rowAt p i)))
        (S.extract (hittingProbabilities p [i]))

-- | The expected return time @m_i = E(T_i | X_0 = i)@, where @T_i@ is the
-- time of the /first return/ to @i@ (at least one step, unlike the hitting
-- time of @{i}@, which is zero when starting there).
--
-- /Theorem (first-return decomposition)./
-- @m_i = 1 + sum_j P_ij eta_j{i}@, with the convention @0 * infinity = 0@:
-- successors with zero probability do not contribute, and any successor with
-- positive probability and infinite expected hitting time makes @m_i@
-- infinite.
--
-- /Proof./ Condition on the first step: from @i@ the chain moves to @j@ with
-- probability @P_ij@, after which the time to reach @i@ is the hitting time
-- of @{i}@ from @j@ (with @eta_i{i} = 0@ covering a self-loop). Adding the
-- one step already taken gives the formula.
--
-- For a finite chain @m_i@ is finite exactly when @i@ is recurrent
-- ('Dtmc.Classification.recurrentState'): a recurrent class is closed and
-- finite, so every state in it hits @i@ almost surely in uniformly bounded
-- expected time, while a transient @i@ has a successor with
-- @h_j{i} < 1@ (otherwise the return probability would be one). The spec
-- verifies this equivalence -- combinatorics against linear algebra -- on
-- random chains.
expectedReturnTime ::
    (KnownNat n) =>
    TransitionMatrix n ->
    Finite n ->
    MeanHittingTime
expectedReturnTime p i =
    foldl' addTerm (FiniteMean 1) (zip finites row)
  where
    eta = expectedHittingTimes p [i]
    row = LA.toList (S.extract (unDistribution (rowAt p i)))
    addTerm acc (j, pij)
        | pij <= 0 = acc
        | otherwise = case (acc, eta j) of
            (FiniteMean s, FiniteMean e) -> FiniteMean (s + pij * e)
            _ -> InfiniteMean
