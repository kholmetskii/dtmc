{- |
Module      : Dtmc.Analysis.Limiting
Description : Limits of the n-step transition matrix.

Long-run behaviour of @P^n@ for a finite chain. For a target @j@ in a
recurrent class @C@ that is aperiodic,

@lim_(n -> infinity) (P^n)(i,j) = h(i,C) pi^C(j)@,

where @h(i,C)@ is the probability of ever entering @C@ from @i@ and @pi^C@ is
the stationary distribution carried by @C@. A transient target has limit
zero. The limit therefore exists exactly when every recurrent class is
aperiodic; a periodic class makes the entries oscillate forever.

For an irreducible chain of period @d@ the whole matrix still has @d@
subsequential limits, one for each residue of @n@ modulo @d@, and
'cyclicLimits' returns them.

Both results assemble already-solved quantities -- eventual hitting
probabilities and per-class stationary distributions -- rather than powering
the matrix, so no truncation or convergence threshold is involved. Deciding
periodicity is combinatorial and exact; the entries inherit the numerical
behaviour of the underlying solves.
-}
module Dtmc.Analysis.Limiting (
    LinearSystemError (..),
    converges,
    limitingMatrix,
    cyclicLimits,
) where

import Data.Finite (
    getFinite,
 )
import Data.Map.Strict qualified as Map
import Data.Maybe (
    fromMaybe,
 )
import Dtmc.Analysis.Classification (
    Irreducible,
    classClosed,
    classPeriod,
    classesOf,
    classify,
    cyclicClasses,
    unIrreducible,
 )
import Dtmc.Analysis.HittingTime qualified as Hitting
import Dtmc.Analysis.LinearSystem (
    LinearSystemError (..),
 )
import Dtmc.Analysis.Stationary (
    classStationaryDistributions,
    stationaryDistribution,
 )
import Dtmc.Distribution.Vector (
    unDistributionVector,
 )
import Dtmc.State (
    FiniteState,
    finiteStates,
    stateIndex,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

toIndex :: (FiniteState state) => state -> Int
toIndex = fromIntegral . getFinite . stateIndex

{- | Whether @P^n@ converges entrywise, that is whether every recurrent class
is aperiodic. Transient classes are irrelevant: their columns tend to zero
whatever their period.

The test is combinatorial, taken from the support graph, so it involves no
arithmetic and no tolerance. The empty chain converges vacuously.

Time: @O(n^2 + n log n + E)@ on an unforced matrix, @O(c)@ afterwards.
-}
converges :: (FiniteState state) => TransitionMatrix state -> Bool
converges p =
    all
        ((== Just 1) . classPeriod)
        [c | c <- classesOf (classify p), classClosed c]

{- | The entrywise limit of @P^n@, with rows and columns in the canonical
order of the 'FiniteState' instance. 'Nothing' when the limit does not exist,
which by 'converges' is exactly when some recurrent class is periodic.

Entry @(i,j)@ is @h(i,C) pi^C(j)@ for @j@ in the recurrent class @C@, and an
exact zero for transient @j@ -- the class distributions vanish there, so the
zero costs no arithmetic. Each recurrent class contributes one hitting solve
and one stationary solve.

Rows sum to one mathematically: from any state the chain enters some
recurrent class almost surely.

Time: @O(n^3)@. Result space: @O(n^2)@.
-}
limitingMatrix ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError (Maybe [[Double]])
limitingMatrix p
    | not (converges p) = Right Nothing
    | otherwise = do
        classes <- classStationaryDistributions p
        contributions <- traverse contribution classes
        pure (Just (foldr addMatrices zeros contributions))
  where
    states = finiteStates :: [state]
    dim = length states
    zeros = replicate dim (replicate dim 0)
    addMatrices = zipWith (zipWith (+))
    contribution (members, classDistribution) = do
        entering <-
            traverse
                (Hitting.eventualProbabilityGivenInitialState p members)
                states
        let
            mass = LA.toList (S.extract (unDistributionVector classDistribution))
        pure [[h * m | m <- mass] | h <- entering]

{- | The @d@ subsequential limits of an irreducible chain of period @d@:
element @r@ is @lim_(n -> infinity) P^(n d + r)@. An aperiodic chain has
@d = 1@ and a single element, which is then the ordinary limit.

Writing @A_0, ..., A_(d-1)@ for the cyclic classes, entry @(i,j)@ of the
@r@-th limit is @d pi(j)@ when @j@ lies @r@ steps after @i@ in the cycle,
that is when the phase of @j@ minus the phase of @i@ is @r@ modulo @d@, and
exactly zero otherwise. The factor @d@ is there because the mass that
stationarity spreads over the whole space is concentrated on one cyclic class
at a time.

Time: @O(n^3)@. Result space: @O(d n^2)@.
-}
cyclicLimits ::
    (FiniteState state) =>
    Irreducible state ->
    Either LinearSystemError [[[Double]]]
cyclicLimits witness = do
    stationary <- stationaryDistribution witness
    let mass = LA.toList (S.extract (unDistributionVector stationary))
        entry r i j
            | (phaseOf j - phaseOf i) `mod` period == r =
                fromIntegral period * (mass !! toIndex j)
            | otherwise = 0
    pure
        [ [[entry r i j | j <- states] | i <- states]
        | r <- [0 .. period - 1]
        ]
  where
    states = finiteStates
    -- A closed class always contains a cycle, so an irreducible chain always
    -- has a defined period; the fallback is the aperiodic reading and cannot
    -- be reached through a valid witness.
    cyclic = fromMaybe [states] (cyclicClasses (unIrreducible witness))
    period = length cyclic
    phases = Map.fromList [(s, r) | (r, group) <- zip [0 ..] cyclic, s <- group]
    phaseOf s = Map.findWithDefault 0 s phases
