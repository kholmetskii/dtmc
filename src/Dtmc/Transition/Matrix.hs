{- |
Module      : Dtmc.Transition.Matrix
Description : Row-stochastic matrices over finite state types.

One-step transition probabilities for a DTMC over a 'FiniteState' type.
'mkTransitionMatrix' validates each row with the simplex tolerance;
'mulTransitionMatrix', 'identityMatrix', and 'matrixPower' provide multi-step
transitions.
-}
module Dtmc.Transition.Matrix (
    TransitionMatrix,
    TransitionMatrixError (..),
    mkTransitionMatrix,
    unTransitionMatrix,
    mulTransitionMatrix,
    identityMatrix,
    matrixPower,
    rowAt,
) where

import Data.Bifunctor (
    first,
 )
import Data.Foldable (
    traverse_,
 )
import Data.Semigroup (
    mtimesDefault,
 )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector,
 )
import Dtmc.Simplex (
    SimplexError,
 )
import Dtmc.Simplex.Internal (
    validateSimplex,
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
 )
import Dtmc.Transition.Matrix.Internal (
    TransitionMatrix,
    matrixRowAt,
    unTransitionMatrix,
    unsafeTransitionMatrix,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )

{- | A row failed simplex validation. Carries its zero-based index and the
underlying error, whose coordinate index is the zero-based column.
-}
data TransitionMatrixError
    = -- | The row index and its simplex failure.
      InRow Int SimplexError
    deriving (Eq, Show)

{- | Validate every row as a next-state distribution, stopping at the first
failure. Accepted entries are preserved without clamping or renormalisation;
the support graph remains lazy. The empty @0 x 0@ matrix is accepted.

Time: @O(n^2)@. Additional validation space: @O(n)@.
-}
mkTransitionMatrix ::
    (FiniteState state) =>
    S.Sq (Cardinality state) ->
    Either TransitionMatrixError (TransitionMatrix state)
mkTransitionMatrix matrix =
    unsafeTransitionMatrix matrix <$ traverse_ validateRow (zip [0 ..] (S.toRows matrix))
  where
    validateRow (index, row) =
        first (InRow index) (validateSimplex row)

{- | Compose two transitions: @mulTransitionMatrix p q@ means take a @p@ step,
then a @q@ step, and stores the matrix product @P Q@.

The product is not revalidated. Exact row-stochastic matrices are closed under
multiplication, but tolerated input error and floating-point rounding can
accumulate, so the result may fail 'mkTransitionMatrix' if checked again.

Time: @O(n^3)@. Result space: @O(n^2)@; its support graph is built lazily.
-}
mulTransitionMatrix ::
    (FiniteState state) =>
    TransitionMatrix state ->
    TransitionMatrix state ->
    TransitionMatrix state
mulTransitionMatrix = (<>)

{- | The @n x n@ identity: the zero-step transition that leaves every state
unchanged. For @n = 0@ this is the empty matrix.

Time and result space: @O(n^2)@.
-}
identityMatrix :: (FiniteState state) => TransitionMatrix state
identityMatrix = mempty

{- | The @k@-step transition matrix @p^k@. Exponent zero returns
'identityMatrix'; positive exponents use repeated squaring through
'Data.Semigroup.mtimesDefault'.

Chapman-Kolmogorov gives @p^(m+n) = p^m p^n@ mathematically; computed matrices
may differ by floating-point rounding and are not revalidated.

Time: @O(n^2 + n^3 log(k + 1))@.
-}
matrixPower ::
    (FiniteState state) =>
    Natural ->
    TransitionMatrix state ->
    TransitionMatrix state
matrixPower = mtimesDefault

{- | The stored row for a state: its next-state distribution. 'FiniteState'
indexing makes the lookup total. The row is wrapped without revalidation, so
any floating-point drift from matrix arithmetic is preserved.
-}
rowAt ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    DistributionVector state
rowAt = matrixRowAt
