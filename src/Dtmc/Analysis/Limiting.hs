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
that decomposition to @P^d@ and rotates its recurrent phase classes. No
truncation or convergence threshold is involved. Class periods are decided
combinatorially; entries inherit the numerical behaviour of the underlying
solves and matrix products.

Unless stated otherwise, complexity bounds exclude 'FiniteState' method
costs. For those bounds, @n@ is the state count, @E@ the support-edge count,
@c@ the number of recurrent classes, and @d@ the least common multiple of
their periods.
-}
module Dtmc.Analysis.Limiting (
    LinearSystemError (..),
    converges,
    limitingMatrix,
    cyclicLimits,
) where

import Data.Array.Unboxed qualified as Unboxed
import Data.List qualified as List
import Dtmc.Analysis.Classification (
    classClosed,
    classMembers,
    classPeriod,
    classesOf,
    classify,
    transientStates,
 )
import Dtmc.Analysis.LinearSystem (
    LinearSystemError (..),
 )
import Dtmc.Analysis.LinearSystem.Internal (
    rowSums,
    solveIminusQ,
    subMatrix,
 )
import Dtmc.Analysis.Stationary (
    stationaryDistributions,
 )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector,
    unDistributionVector,
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.State.Internal (
    stateCardinalityInt,
    stateIndexInt,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    power,
 )
import Dtmc.Transition.Matrix.Internal (
    unTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )

{- | Test whether @P^n@ converges entrywise, equivalently whether every
recurrent class is aperiodic. Transient classes are irrelevant because their
columns tend to zero whatever their period.

The test is combinatorial, taken from the support graph, so it involves no
arithmetic and no tolerance. The empty chain converges vacuously.

Complexity: on an unforced matrix, the first query takes
@O(n^2 + (n + E) log(n + 1))@ time and @O(n^2 + n + E)@ temporary space and
retains @O(n + E)@ graph and classification caches. With those facts cached,
it takes @O(c)@ time, @O(c)@ temporary space, and @O(1)@ result space.
-}
converges :: (FiniteState state) => TransitionMatrix state -> Bool
converges p =
    all
        ((== Just 1) . classPeriod)
        [c | c <- classesOf (classify p), classClosed c]

{- | Compute the entrywise limit of @P^n@, with rows and columns in the
canonical order of the 'FiniteState' instance. Return 'Nothing' exactly when
some recurrent class is periodic and the limit does not exist.

Entry @(i,j)@ is @h(i,C) pi^C(j)@ for @j@ in the recurrent class @C@, and an
exact zero for transient @j@ -- the class distributions vanish there, so the
zero costs no arithmetic. All recurrent-class hitting probabilities are the
columns of one linear system; each recurrent class still contributes one
stationary solve.

Rows sum to one mathematically because a finite chain enters some recurrent
class almost surely. Numerical failures come from the checked hitting and
stationary solves. A periodic chain returns 'Nothing' without those solves.

Complexity: @O(n^3)@ worst-case time, @O(n^2)@ temporary space,
@O(n + E)@ retained graph-cache space, and @O(n^2)@ result space.
-}
limitingMatrix ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError (Maybe [[Double]])
limitingMatrix p
    | not (converges p) = Right Nothing
    | otherwise = Just <$> convergentLimit p

{- | Compute the limiting matrix of a chain already known to have only
aperiodic recurrent classes. Keeping this separate lets 'cyclicLimits' apply
the same class-and-hitting decomposition to @P^d@ without a redundant
convergence branch.

The caller is responsible for the aperiodicity precondition. Numerical
failures come from the checked hitting and stationary solves.

Complexity: @O(n^3)@ worst-case time, @O(n^2)@ temporary space,
@O(n + E)@ retained graph-cache space, and @O(n^2)@ result space.
-}
convergentLimit ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError [[Double]]
convergentLimit p = do
    (classes, entering) <- limitDecomposition p
    if null classes
        then pure (replicate dim (replicate dim 0))
        else
            pure
                ( LA.toLists
                    ( entering
                        LA.<> LA.fromRows
                            [ S.extract (unDistributionVector distribution)
                            | (_, distribution) <- classes
                            ]
                    )
                )
  where
    dim = stateCardinalityInt @state

{- | Compute the recurrent-class stationary distributions and the matrix
@H@ whose entry @(i,C)@ is the probability of eventually entering class @C@
from @i@.

Writing @T@ for the transient states, all transient rows are obtained from the
single multiple-right-hand-side system

@(I - P[T,T]) H[T,*] = B@,

where @B(i,C) = sum_(j in C) P(i,j)@. Recurrent rows are exact zero/one
boundary values supplied from classification.

An empty chain returns no classes and an empty entering matrix. Numerical
failures come from the transient hitting solve or a class stationary solve.

Complexity: @O(n^3)@ worst-case time, @O(n^2)@ temporary space,
@O(n + E)@ retained graph-cache space, and @O(c n)@ result space.
-}
limitDecomposition ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either
        LinearSystemError
        ([([state], DistributionVector state)], LA.Matrix Double)
limitDecomposition p = do
    classes <- stationaryDistributions p
    transientSolution <-
        if null transient
            then Right (LA.konst 0 (0, classCount))
            else
                solveIminusQ
                    (subMatrix transientIndices transientIndices matrix)
                    ( LA.fromColumns
                        [ rowSums (subMatrix transientIndices (map toIndex members) matrix)
                        | (members, _) <- classes
                        ]
                    )
    let entering =
            LA.fromLists
                [ [ entryAt stateIndex classIndex
                  | classIndex <- [0 .. classCount - 1]
                  ]
                | stateIndex <- [0 .. dim - 1]
                ]
        entryAt stateIndex classIndex
            | transientPosition Unboxed.! stateIndex >= 0 =
                transientSolution
                    `LA.atIndex` (transientPosition Unboxed.! stateIndex, classIndex)
            | otherwise =
                if recurrentClass Unboxed.! stateIndex == classIndex then 1 else 0
    pure (classes, entering)
  where
    dim = stateCardinalityInt @state
    matrix = S.extract (unTransitionMatrix p)
    transient = transientStates p
    transientIndices = map toIndex transient
    closedClasses =
        [ classMembers recurrentClass'
        | recurrentClass' <- classesOf (classify p)
        , classClosed recurrentClass'
        ]
    classCount = length closedClasses

    transientPosition :: Unboxed.UArray Int Int
    transientPosition =
        Unboxed.accumArray
            (\_ position -> position)
            (-1)
            (0, dim - 1)
            (zip transientIndices [0 ..])

    recurrentClass :: Unboxed.UArray Int Int
    recurrentClass =
        Unboxed.accumArray
            (\_ classIndex -> classIndex)
            (-1)
            (0, dim - 1)
            [ (toIndex member, classIndex)
            | (classIndex, members) <- zip [0 ..] closedClasses
            , member <- members
            ]

    toIndex = stateIndexInt

{- | Compute the @d@ subsequential limits of any finite chain, where @d@ is
the least common multiple of its recurrent class periods. Element @r@ is
@lim_(n -> infinity) P^(n d + r)@. When every recurrent class is aperiodic,
@d = 1@ and the single result is the ordinary limiting matrix. The empty chain
also returns one empty matrix.

The powered chain @Q = P^d@ has only aperiodic recurrent classes. Its closed
classes are the cyclic phases of the original recurrent classes. One batched
solve finds the probability of entering every phase at times divisible by
@d@. A transition under @P@ permutes those phase classes, so later residues
are assembled by rotating class indices rather than multiplying dense
matrices. This accounts automatically for multiple recurrent classes,
transient-state hitting probabilities, and entry phases.

Numerical failures come from the checked hitting and stationary solves. An
inconsistent phase-class successor relation produces 'SingularSystem'.

Complexity: @O(n^3 log(d + 1) + d n^2)@ time, @O(n^2)@ temporary space,
@O(n + E)@ retained graph-cache space, and @O(d n^2)@ result space.
-}
cyclicLimits ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError [[[Double]]]
cyclicLimits p
    | dim == 0 = Right [[]]
    | otherwise = do
        (phaseClasses, entering) <- limitDecomposition powered
        let classCount = length phaseClasses
            classBounds = (0, classCount - 1)
            phaseClassByState = classByState phaseClasses
        successors <- traverse (successorOf phaseClassByState . fst) phaseClasses
        if List.sort successors /= [0 .. classCount - 1]
            then Left SingularSystem
            else
                pure
                    ( successiveLimits
                        commonPeriod
                        classBounds
                        classCount
                        (Unboxed.listArray classBounds [0 .. classCount - 1])
                        ( Unboxed.array
                            classBounds
                            [(successor, source) | (source, successor) <- zip [0 ..] successors]
                        )
                        entering
                        phaseClassByState
                        (stationaryMass phaseClasses)
                    )
  where
    dim = stateCardinalityInt @state
    original = S.extract (unTransitionMatrix p)
    powered = power commonPeriod p
    commonPeriod =
        foldr
            lcm
            1
            [ classPeriodValue
            | recurrentClass <- classesOf (classify p)
            , classClosed recurrentClass
            , Just classPeriodValue <- [classPeriod recurrentClass]
            ]

    classByState phaseClasses =
        Unboxed.accumArray
            (\_ classIndex -> classIndex)
            (-1)
            (0, dim - 1)
            [ (stateIndexInt member, classIndex)
            | (classIndex, (members, _)) <- zip [0 ..] phaseClasses
            , member <- members
            ]

    stationaryMass phaseClasses =
        Unboxed.accumArray
            (\_ probability -> probability)
            0
            (0, dim - 1)
            [ (memberIndex, vector `LA.atIndex` memberIndex)
            | (members, distribution) <- phaseClasses
            , let vector = S.extract (unDistributionVector distribution)
            , member <- members
            , let memberIndex = stateIndexInt member
            ]

    successorOf phaseClassByState members =
        case destinations of
            successor : rest
                | successor >= 0 && all (== successor) rest -> Right successor
            _ -> Left SingularSystem
      where
        destinations =
            [ phaseClassByState Unboxed.! destination
            | member <- members
            , let source = stateIndexInt member
            , destination <- [0 .. dim - 1]
            , original `LA.atIndex` (source, destination) > 0
            ]

    successiveLimits ::
        Natural ->
        (Int, Int) ->
        Int ->
        Unboxed.UArray Int Int ->
        Unboxed.UArray Int Int ->
        LA.Matrix Double ->
        Unboxed.UArray Int Int ->
        Unboxed.UArray Int Double ->
        [[[Double]]]
    successiveLimits 0 _ _ _ _ _ _ _ = []
    successiveLimits remaining classBounds classCount origins predecessors entering classOf mass =
        [ [ limitEntry initial target
          | target <- [0 .. dim - 1]
          ]
        | initial <- [0 .. dim - 1]
        ]
            : successiveLimits
                (remaining - 1)
                classBounds
                classCount
                ( Unboxed.listArray
                    classBounds
                    [ predecessors Unboxed.! (origins Unboxed.! classIndex)
                    | classIndex <- [0 .. classCount - 1]
                    ]
                )
                predecessors
                entering
                classOf
                mass
      where
        limitEntry initial target
            | targetClass < 0 = 0
            | otherwise =
                entering `LA.atIndex` (initial, origins Unboxed.! targetClass)
                    * (mass Unboxed.! target)
          where
            targetClass = classOf Unboxed.! target
