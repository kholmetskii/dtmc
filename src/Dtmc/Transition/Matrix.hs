{- |
Module      : Dtmc.Transition.Matrix
Description : Row-stochastic matrices over finite state types.

One-step transition probabilities for a DTMC over a 'FiniteState' type.
'fromKernel' turns an already-validated finite-state kernel into a dense
matrix; 'mulTransitionMatrix', 'identityMatrix', and 'matrixPower' provide
multi-step transitions. Explicit @hmatrix@ interoperability lives in
"Dtmc.Transition.Matrix.HMatrix".
-}
module Dtmc.Transition.Matrix (
    -- * Representation
    TransitionMatrix,

    -- * Construction and inspection
    fromKernel,
    toRows,
    rowAt,

    -- * Composition
    mulTransitionMatrix,
    identityMatrix,
    matrixPower,
) where

import Data.Semigroup (
    mtimesDefault,
 )
import Dtmc.Distribution.Map.Internal (
    denseWeights,
 )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector,
 )
import Dtmc.State (
    FiniteState,
    finiteStates,
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
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )

{- | Materialise a finite-state kernel as a dense transition matrix. Kernel
rows are already validated 'Dtmc.Distribution.Map.DistributionMap' values, so
this conversion is total and performs no additional clamping or
renormalisation. Missing coordinates become exact zeros. The support graph
remains lazy, and the empty @0 x 0@ matrix is accepted.

Time and result space: @O(n^2)@.
-}
fromKernel ::
    (FiniteState state) =>
    TransitionKernel state ->
    TransitionMatrix state
fromKernel kernel =
    unsafeTransitionMatrix $
        S.matrix
            [ weight
            | source <- finiteStates
            , let distribution = transitionLaw kernel source
            , weight <- denseWeights finiteStates distribution
            ]

{- | Return all stored entries as rows in canonical state order. Exact zeros
are retained. This is a representation-neutral copy of the dense matrix and
does not force its support graph.

Time and result space: @O(n^2)@.
-}
toRows :: (FiniteState state) => TransitionMatrix state -> [[Double]]
toRows = LA.toLists . S.extract . unTransitionMatrix

{- | Compose two transitions: @mulTransitionMatrix p q@ means take a @p@ step,
then a @q@ step, and stores the matrix product @P Q@.

The product is not revalidated. Row-stochastic matrices are closed under
multiplication mathematically, but floating-point rounding can accumulate.

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
