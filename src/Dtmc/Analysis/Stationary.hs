{- |
Module      : Dtmc.Analysis.Stationary
Description : Stationary distributions of finite chains.

Every finite irreducible DTMC has exactly one stationary distribution,
including periodic chains. A reducible chain has one extremal stationary
distribution per recurrent class and no canonical choice among them.
'stationaryDistributions' returns these class-supported distributions. Every
stationary distribution of the chain is a convex combination of them, and a
chain with two or more recurrent classes therefore has infinitely many.

Each recurrent class is solved by Grassmann-Taksar-Heyman state reduction.
GTH never forms @transpose(P) - I@ and performs no subtraction at all: every
step adds, multiplies, or divides non-negative quantities, avoiding the
cancellation that damages a balance solve on a nearly uncoupled chain. It has
the same @O(m^3)@ asymptotic cost as an LU factorisation.

Every supplied block is irreducible because a closed communicating class
restricted to itself is irreducible. This is exactly the condition under which
no elimination step can divide by zero. An 'IllConditionedSystem' is therefore
unreachable here, unlike in the hitting- and return-time solves that share the
error type. Results are not clamped or renormalised; a non-finite input or
solution, and a residual @|pi P - pi|@ above @1e-9@, are still reported
explicitly. Complexity bounds exclude 'FiniteState' method costs. For the
top-level bound, @n@ is the state count and @E@ the support-edge count.
-}
module Dtmc.Analysis.Stationary (
    LinearSystemError (..),
    stationaryDistributions,
) where

import Control.Monad.ST (
    ST,
    runST,
 )
import Data.Array.MArray (
    newListArray,
    readArray,
    writeArray,
 )
import Data.Array.ST (
    STUArray,
 )
import Data.Array.Unboxed qualified as Unboxed
import Dtmc.Analysis.Classification (
    classClosed,
    classMembers,
    classesOf,
    classify,
 )
import Dtmc.Analysis.LinearSystem (
    LinearSystemError (..),
 )
import Dtmc.Analysis.LinearSystem.Internal (
    subMatrix,
 )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector (DistributionVector),
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
 )
import Dtmc.Transition.Matrix.Internal (
    unTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

toIndex :: (FiniteState state) => state -> Int
toIndex = stateIndexInt

{- | Compute the stationary vector of a non-empty irreducible stochastic
block, in the block's own ordering, by Grassmann-Taksar-Heyman state
reduction.

The reduction removes states one at a time. Censoring the chain on
@{1, ..., k-1}@ -- watching it only when it is outside @k@ -- leaves a Markov
chain with the same stationary distribution up to normalisation, and its
transition probabilities are

@P'(i,j) = P(i,j) + P(i,k) P(k,j) / S,   S = sum_(j < k) P(k,j)@,

read as: reach @j@ directly, or by entering @k@ and leaving it at @j@. This is
the same elimination order as Gaussian elimination, arranged so that no
subtraction appears. In particular @S@ is accumulated from the off-diagonal
entries rather than as @1 - P(k,k)@, which is what saves a state that holds
probability close to one: the difference would lose most of its significant
digits while the sum loses none.

The back substitution reads @x(k) = sum_(i < k) x(i) P(i,k)@ off the stored
factors, @x(k)@ being the expected number of visits to @k@ per visit to the
first state in the chain censored on @{1, ..., k}@ -- non-negative, so it does
not cancel either. Normalising gives @pi@.

An irreducible block makes @S > 0@ at every step because the censored chain is
again irreducible and its last state must be able to leave. A reducible block
that reaches @S = 0@ is reported as 'SingularSystem'. Non-finite inputs,
weights, or solutions and an excessive stationarity residual produce the
corresponding 'LinearSystemError'.

Complexity: @O(m^3)@ time, @O(m^2)@ temporary space, and @O(m)@ result
space for an @m x m@ block.
-}
stationaryOfBlock ::
    LA.Matrix Double ->
    Either LinearSystemError (LA.Vector Double)
stationaryOfBlock block
    | dimension == 0 = Left SingularSystem
    | not (all isFinite (LA.toList (LA.flatten block))) = Left NonFiniteSystem
    | otherwise =
        case reduceGth dimension (LA.toList (LA.flatten block)) of
            Nothing -> Left SingularSystem
            Just weights -> normalise weights
  where
    dimension = LA.rows block
    normalise weights
        | not (all isFinite weights) = Left NonFiniteSolution
        | total <= 0 = Left SingularSystem
        | residual > limit =
            Left
                ( ResidualTooLarge
                    { relativeResidual = residual
                    , residualLimit = limit
                    }
                )
        | otherwise = Right stationary
      where
        limit = 1e-9
        total = sum weights
        stationary = LA.fromList (map (/ total) weights)
        residual =
            foldr (max . abs) 0 (LA.toList (LA.tr block LA.#> stationary - stationary))

isFinite :: Double -> Bool
isFinite value = not (isNaN value || isInfinite value)

-- Allocating through a signature that quantifies the state thread keeps the
-- array type unambiguous without local annotations inside 'runST'.
newFlatArray :: (Int, Int) -> [Double] -> ST s (STUArray s Int Double)
newFlatArray = newListArray

{- | Compute the unnormalised GTH weights of an @n x n@ block supplied in
row-major order. Return 'Nothing' when an elimination step finds no positive
way out of the state being removed.

The caller must supply exactly @n^2@ entries. This helper performs no
finiteness, stochasticity, or irreducibility validation.

Complexity: @O(n^3)@ time, @O(n^2)@ temporary space, and @O(n)@ result
space.
-}
reduceGth :: Int -> [Double] -> Maybe [Double]
reduceGth n entries = runST $ do
    a <- newFlatArray (0, n * n - 1) entries
    let index i j = i * n + j

        exitMass k =
            sum <$> mapM (readArray a . index k) [0 .. k - 1]

        absorbRow k s i = do
            entering <- readArray a (index i k)
            let scaled = entering / s
            writeArray a (index i k) scaled
            mapM_
                ( \j -> do
                    leaving <- readArray a (index k j)
                    current <- readArray a (index i j)
                    writeArray a (index i j) (current + scaled * leaving)
                )
                [0 .. k - 1]

        eliminate k
            | k < 1 = pure True
            | otherwise = do
                s <- exitMass k
                if s <= 0
                    then pure False
                    else do
                        mapM_ (absorbRow k s) [0 .. k - 1]
                        eliminate (k - 1)

        -- visits holds x(0) .. x(k-1) in order.
        substitute k visits
            | k >= n = pure visits
            | otherwise = do
                terms <-
                    mapM
                        (\(i, x) -> (x *) <$> readArray a (index i k))
                        (zip [0 ..] visits)
                substitute (k + 1) (visits ++ [sum terms])

    feasible <- eliminate (n - 1)
    if feasible
        then Just <$> substitute 1 [1]
        else pure Nothing

{- | Compute the extremal stationary distributions, one per recurrent class
and paired with the class on which each lives. Classes come in the order of
'Dtmc.Analysis.Classification.classify', that is by least member. An empty
chain returns an empty list.

Each distribution is returned over the whole state space, carrying exact zeros
outside its class. This is correct rather than merely convenient: a recurrent
class is closed, so a distribution supported on it satisfies @pi P = pi@ for
the full matrix, and no stationary distribution of a finite chain puts mass on
a transient state.

Every stationary distribution of the chain is a convex combination of these.
The result has exactly one element precisely when the stationary distribution
is unique, including reducible chains with transient states and one recurrent
class. It has two or more elements exactly when the chain has infinitely many
stationary distributions.

Each class is solved separately. The first numerical failure aborts the
traversal, although the error does not identify its class.

Complexity: @O(n^2 + sum_C |C|^3)@ time, at most @O(n^3)@, and @O(n^2)@
temporary space. Result space is @O(c n)@ for @c@ recurrent classes. The
matrix may retain @O(n + E)@ graph-cache space.
-}
stationaryDistributions ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError [([state], DistributionVector state)]
stationaryDistributions p =
    traverse distributionOn closedClasses
  where
    dim = stateCardinalityInt @state
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
