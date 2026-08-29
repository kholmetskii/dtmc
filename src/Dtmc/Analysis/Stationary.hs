{- |
Module      : Dtmc.Analysis.Stationary
Description : Stationary distributions of finite chains.

Every finite irreducible DTMC has exactly one stationary distribution,
including periodic chains: 'stationaryDistribution' computes it from an
'Irreducible' witness. For a row-stochastic matrix @P@ the result @pi@
satisfies @transpose(P) pi = pi@ and @sum pi = 1@.

A reducible chain has one stationary distribution per recurrent class and no
canonical choice among them, so 'classStationaryDistributions' returns them
all. Every stationary distribution of the chain is a convex combination of
these, and a chain with two or more recurrent classes therefore has
infinitely many.

Both entry points share one solve: a dependent balance equation is replaced
by the normalization equation and the system goes through a checked LU
factorisation. Results are not clamped or renormalised. Non-finite, singular,
ill-conditioned, and high-residual systems return an explicit
'LinearSystemError'.
-}
module Dtmc.Analysis.Stationary (
    LinearSystemError (..),
    stationaryDistribution,
    classStationaryDistributions,
) where

import Data.Array.Unboxed qualified as Unboxed
import Data.Finite (
    getFinite,
 )
import Data.Proxy (
    Proxy (..),
 )
import Dtmc.Analysis.Classification (
    Irreducible,
    classClosed,
    classMembers,
    classesOf,
    classify,
    unIrreducible,
 )
import Dtmc.Analysis.LinearSystem (
    LinearSystemError (..),
 )
import Dtmc.Analysis.LinearSystem.Internal (
    solveLinearSystem,
    subMatrix,
 )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector (DistributionVector),
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
    stateIndex,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    unTransitionMatrix,
 )
import GHC.TypeNats (
    natVal,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

toIndex :: (FiniteState state) => state -> Int
toIndex = fromIntegral . getFinite . stateIndex

{- | The stationary vector of a non-empty stochastic block, as a dynamically
sized vector in the block's own ordering.

One balance equation is dependent on the others, so the last row of the
system is replaced by @sum pi = 1@. A @1 x 1@ block reduces to @1 * x = 1@
and needs no balance row at all.

Time: @O(m^3)@ for an @m x m@ block.
-}
stationaryOfBlock ::
    LA.Matrix Double ->
    Either LinearSystemError (LA.Vector Double)
stationaryOfBlock block =
    LA.flatten <$> solveLinearSystem coefficient rightHandSide
  where
    dimension = LA.rows block
    balanceRows = LA.toLists (LA.tr block - LA.ident dimension)
    coefficient =
        LA.fromLists
            ( take (dimension - 1) balanceRows
                ++ [replicate dimension 1]
            )
    rightHandSide =
        LA.fromLists
            ( replicate (dimension - 1) [0]
                ++ [[1]]
            )

{- | The unique stationary distribution of a certified finite irreducible
chain. The witness cannot exist for an empty chain. Coordinates follow the
canonical state order of the 'FiniteState' instance.

Time: @O(n^3)@. Space: @O(n^2)@.
-}
stationaryDistribution ::
    (FiniteState state) =>
    Irreducible state ->
    Either LinearSystemError (DistributionVector state)
stationaryDistribution witness = do
    solution <-
        stationaryOfBlock (S.extract (unTransitionMatrix (unIrreducible witness)))
    pure (DistributionVector (S.vector (LA.toList solution)))

{- | One stationary distribution per recurrent class, paired with the class it
lives on. Classes come in the order of 'Dtmc.Analysis.Classification.classify',
that is by least member.

Each distribution is returned over the whole state space, carrying exact zeros
outside its class. This is correct rather than merely convenient: a recurrent
class is closed, so a distribution supported on it satisfies @pi P = pi@ for
the full matrix, and no stationary distribution of a finite chain puts mass on
a transient state.

Every stationary distribution of the chain is a convex combination of these,
so the result has exactly one element for an irreducible chain -- equal to
'stationaryDistribution' of its witness -- and two or more elements exactly
when the chain has infinitely many stationary distributions.

Each class is solved separately, so a numerical failure names the whole call
rather than one class.

Time: @O(n^2 + sum_C |C|^3)@, which is at most @O(n^3)@. Result space:
@O(c n)@ for @c@ recurrent classes.
-}
classStationaryDistributions ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError [([state], DistributionVector state)]
classStationaryDistributions p =
    traverse distributionOn closedClasses
  where
    dim = fromIntegral (natVal (Proxy @(Cardinality state)))
    matrix = S.extract (unTransitionMatrix p)
    closedClasses =
        [classMembers c | c <- classesOf (classify p), classClosed c]
    distributionOn members = do
        solution <- stationaryOfBlock (subMatrix indices indices matrix)
        let placed :: Unboxed.UArray Int Double
            placed =
                Unboxed.accumArray
                    (\_ x -> x)
                    0
                    (0, dim - 1)
                    (zip indices (LA.toList solution))
        pure
            ( members
            , DistributionVector (S.vector [placed Unboxed.! i | i <- [0 .. dim - 1]])
            )
      where
        indices = map toIndex members
