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

For any finite chain, let @d@ be the least common multiple of its recurrent
class periods. The whole matrix has @d@ subsequential limits, one for each
residue of @n@ modulo @d@, and 'cyclicLimits' returns them. An empty chain uses
@d = 1@.

'limitingMatrix' assembles eventual hitting probabilities and per-class
stationary distributions without powering the matrix. 'cyclicLimits' applies
that convergent calculation to @P^d@ and advances each limit by one transition.
No truncation or convergence threshold is involved. Class periods are decided
combinatorially; entries inherit the numerical behaviour of the underlying
solves and matrix products.
-}
module Dtmc.Analysis.Limiting (
    LinearSystemError (..),
    converges,
    limitingMatrix,
    cyclicLimits,
) where

import Dtmc.Analysis.Classification (
    classClosed,
    classPeriod,
    classesOf,
    classify,
 )
import Dtmc.Analysis.HittingTime qualified as Hitting
import Dtmc.Analysis.LinearSystem (
    LinearSystemError (..),
 )
import Dtmc.Analysis.Stationary (
    stationaryDistributions,
 )
import Dtmc.Distribution.Vector (
    unDistributionVector,
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
    finiteStates,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    matrixPower,
    unTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )

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
    | otherwise = Just <$> convergentLimit p

{- | The limiting matrix of a chain already known to have only aperiodic
recurrent classes. Keeping this separate lets 'cyclicLimits' apply the same
class-and-hitting decomposition to @P^d@ without a redundant convergence
branch.
-}
convergentLimit ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError [[Double]]
convergentLimit p = do
    classes <- stationaryDistributions p
    contributions <- traverse contribution classes
    pure (foldr addMatrices zeros contributions)
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

{- | The @d@ subsequential limits of any finite chain, where @d@ is the least
common multiple of its recurrent class periods: element @r@ is
@lim_(n -> infinity) P^(n d + r)@. When every recurrent class is aperiodic,
@d = 1@ and the single result is the ordinary limiting matrix. The empty chain
also returns one empty matrix.

The powered chain @Q = P^d@ has only aperiodic recurrent classes. Its ordinary
limit is the residue-zero result; right-multiplying successively by @P@ gives
the remaining residues. This accounts automatically for multiple recurrent
classes, transient-state hitting probabilities, and entry phases.

Time: @O(n^3 (d + log(d + 1)))@. Result space: @O(d n^2)@.
-}
cyclicLimits ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError [[[Double]]]
cyclicLimits p = do
    atZero <- convergentLimit (matrixPower commonPeriod p)
    let initial = S.matrix (concat atZero) :: S.Sq (Cardinality state)
    pure (successiveLimits commonPeriod initial)
  where
    commonPeriod =
        foldr
            lcm
            1
            [ classPeriodValue
            | recurrentClass <- classesOf (classify p)
            , classClosed recurrentClass
            , Just classPeriodValue <- [classPeriod recurrentClass]
            ]

    successiveLimits ::
        Natural ->
        S.Sq (Cardinality state) ->
        [[[Double]]]
    successiveLimits 0 _ = []
    successiveLimits remaining current =
        LA.toLists (S.extract current)
            : successiveLimits
                (remaining - 1)
                (current S.<> unTransitionMatrix p)
