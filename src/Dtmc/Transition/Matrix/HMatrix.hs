{- |
Module      : Dtmc.Transition.Matrix.HMatrix
Description : Explicit hmatrix interoperability for transition matrices.

Low-level construction and projection for callers already using @hmatrix@.
Most users should construct a 'Dtmc.Transition.Kernel.TransitionKernel' and
convert it with 'Dtmc.Transition.Matrix.fromKernel'; inspect matrices with
'Dtmc.Transition.Matrix.toRows'.
-}
module Dtmc.Transition.Matrix.HMatrix (
    TransitionMatrixError (..),
    mkTransitionMatrix,
    unTransitionMatrix,
) where

import Data.Bifunctor (
    first,
 )
import Dtmc.Simplex (
    SimplexError,
 )
import Dtmc.Simplex.Internal (
    canonicaliseSimplex,
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
 )
import Dtmc.Transition.Matrix.Internal (
    unTransitionMatrix,
    unsafeTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

{- | A row failed simplex validation. Carries its zero-based index and the
underlying error, whose coordinate index is the zero-based column.
-}
data TransitionMatrixError
    = -- | The row index and its simplex failure.
      InRow Int SimplexError
    deriving (Eq, Show)

{- | Construct a row-stochastic matrix from a statically sized @hmatrix@
matrix, stopping at the first invalid row. Within each accepted row, tolerated
coordinate error is clamped to @[0, 1]@ and the repaired row is normalised.
The support graph remains lazy. The empty @0 x 0@ matrix is accepted.

Time and result space: @O(n^2)@.
-}
mkTransitionMatrix ::
    (FiniteState state) =>
    S.Sq (Cardinality state) ->
    Either TransitionMatrixError (TransitionMatrix state)
mkTransitionMatrix matrix =
    unsafeTransitionMatrix . S.matrix . concatMap (LA.toList . S.extract)
        <$> traverse canonicaliseRow (zip [0 ..] (S.toRows matrix))
  where
    canonicaliseRow (index, row) =
        first (InRow index) (canonicaliseSimplex row)
