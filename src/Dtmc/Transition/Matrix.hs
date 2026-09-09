{- |
Module      : Dtmc.Transition.Matrix
Description : Row-stochastic matrices over finite state types.

One-step transition probabilities for a DTMC over a 'FiniteState' type.
'fromRows' builds one from a grid of weights and 'fromKernel' from an
already-validated finite-state kernel; 'compose', 'identity', and 'power'
provide multi-step transitions. 'toRows' reads the stored probabilities back
as plain lists.
-}
module Dtmc.Transition.Matrix (
    -- * Representation
    TransitionMatrix,
    TransitionMatrixError (..),

    -- * Construction and inspection
    fromKernel,
    fromRows,
    toRows,
    rowAt,

    -- * Composition
    compose,
    identity,
    power,
) where

import Data.Bifunctor (
    first,
 )
import Data.Semigroup (
    mtimesDefault,
 )
import Dtmc.Distribution.Map.Internal (
    denseWeights,
 )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector,
 )
import Dtmc.Simplex (
    SimplexError,
 )
import Dtmc.Simplex.Internal (
    canonicaliseSimplexEntries,
 )
import Dtmc.State (
    FiniteState,
    finiteStates,
 )
import Dtmc.State.Internal (
    stateCardinalityInt,
 )
import Dtmc.Transition (
    Transition (transitionLaw),
 )
import Dtmc.Transition.Kernel (
    TransitionKernel,
 )
import Dtmc.Transition.Matrix.Internal (
    TransitionMatrix,
    matrixRowAt,
    unTransitionMatrix,
    unsafeTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.Natural (
    Natural,
 )

{- | Why a supplied grid of weights is not a transition matrix. Row and column
indices are zero-based and follow the canonical state order of the
'FiniteState' instance.
-}
data TransitionMatrixError
    = -- | A row failed simplex validation: its index and the underlying
      -- failure, whose coordinate index is the zero-based column.
      InRow Int SimplexError
    | -- | The state cardinality and the supplied number of rows.
      WrongRowCount Int Int
    | -- | A row of the wrong width: its index, the state cardinality, and the
      -- supplied width.
      WrongRowWidth Int Int Int
    deriving (Eq, Show)

{- | Construct a row-stochastic matrix from a grid of weights in canonical
state order, stopping at the first problem. Within each accepted row,
tolerated coordinate error is clamped to @[0, 1]@ and the repaired row is
normalised. The support graph and classification cache remain lazy, and the
empty @0 x 0@ matrix is accepted.

This inverts 'toRows' up to that repair and needs no @hmatrix@ value: the
shape is checked here and reported as 'WrongRowCount' or 'WrongRowWidth'
rather than raised by the array backend.

Complexity: @O(n^2)@ time and @O(n^2)@ temporary and result space.
-}
fromRows ::
    forall state.
    (FiniteState state) =>
    [[Double]] ->
    Either TransitionMatrixError (TransitionMatrix state)
fromRows rows
    | suppliedRows /= dimension = Left (WrongRowCount dimension suppliedRows)
    | otherwise =
        unsafeTransitionMatrix . (dimension LA.>< dimension) . concat
            <$> traverse canonicaliseRow (zip [0 ..] rows)
  where
    dimension = stateCardinalityInt @state
    suppliedRows = length rows

    canonicaliseRow (index, row)
        | width /= dimension = Left (WrongRowWidth index dimension width)
        | otherwise = first (InRow index) (canonicaliseSimplexEntries row)
      where
        width = length row

{- | Materialise a finite-state kernel as a dense transition matrix. Kernel
rows are already validated 'Dtmc.Distribution.Map.DistributionMap' values, so
this conversion is total and performs no additional clamping or
renormalisation. Missing coordinates become exact zeros. The support graph and
classification cache remain lazy, and the empty @0 x 0@ matrix is accepted.

Complexity: excluding evaluation of 'finiteStates' and the kernel laws,
@O(n^2)@ time and @O(n^2)@ temporary and result space.
-}
fromKernel ::
    forall state.
    (FiniteState state) =>
    TransitionKernel state ->
    TransitionMatrix state
fromKernel kernel =
    unsafeTransitionMatrix $
        (dimension LA.>< dimension)
            [ weight
            | source <- finiteStates
            , let distribution = transitionLaw kernel source
            , weight <- denseWeights finiteStates distribution
            ]
  where
    dimension = stateCardinalityInt @state

{- | Return all stored entries as rows in canonical state order. Exact zeros
are retained. This is a representation-neutral copy of the dense matrix and
does not force its support graph.

Complexity: @O(n^2)@ time and @O(n^2)@ temporary and result space.
-}
toRows :: TransitionMatrix state -> [[Double]]
toRows = LA.toLists . unTransitionMatrix

{- | Compose two transitions: @compose p q@ means take a @p@ step,
then a @q@ step, and stores the matrix product @P Q@.

The product is not revalidated. Row-stochastic matrices are closed under
multiplication mathematically, but floating-point rounding can accumulate.

Complexity: @O(n^3)@ worst-case time and @O(n^2)@ temporary and result space.
The support graph and classification cache are built lazily.
-}
compose ::
    TransitionMatrix state ->
    TransitionMatrix state ->
    TransitionMatrix state
compose = (<>)

{- | Return the @n x n@ identity: the zero-step transition that leaves every
state unchanged. For @n = 0@ this is the empty matrix.

Complexity: @O(1)@ construction time and @O(1)@ construction space. Forcing
the dense entries or support graph takes @O(n^2)@ time and @O(n^2)@ temporary
space; the support graph itself occupies @O(n)@ space.
-}
identity :: (FiniteState state) => TransitionMatrix state
identity = mempty

{- | Compute the @k@-step transition matrix @p^k@. Exponent zero returns
'identity'; positive exponents use repeated squaring through
'Data.Semigroup.mtimesDefault'.

Chapman-Kolmogorov gives @p^(m+n) = p^m p^n@ mathematically; computed matrices
may differ by floating-point rounding and are not revalidated.

Complexity: @O(n^2 + n^3 log(k + 1))@ time and @O(n^2)@ temporary and result
space.
-}
power ::
    (FiniteState state) =>
    Natural ->
    TransitionMatrix state ->
    TransitionMatrix state
power = mtimesDefault

{- | Return the stored row for a state: its next-state distribution.
'FiniteState' indexing makes the lookup total. The row is wrapped without
revalidation, so any floating-point drift from matrix arithmetic is
preserved.

Complexity: excluding 'Dtmc.State.stateIndex', @O(n)@ time and @O(n)@ result
space for state cardinality @n@.
-}
rowAt ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    DistributionVector state
rowAt = matrixRowAt
