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

Complexity: excluding evaluation of 'finiteStates' and the kernel laws,
@O(n^2)@ time and @O(n^2)@ temporary and result space.
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

Complexity: @O(n^2)@ time and @O(n^2)@ temporary and result space.
-}
toRows :: (FiniteState state) => TransitionMatrix state -> [[Double]]
toRows = LA.toLists . S.extract . unTransitionMatrix

{- | Compose two transitions: @mulTransitionMatrix p q@ means take a @p@ step,
then a @q@ step, and stores the matrix product @P Q@.

The product is not revalidated. Row-stochastic matrices are closed under
multiplication mathematically, but floating-point rounding can accumulate.

Complexity: @O(n^3)@ worst-case time and @O(n^2)@ temporary and result space.
The support graph is built lazily.
-}
mulTransitionMatrix ::
    (FiniteState state) =>
    TransitionMatrix state ->
    TransitionMatrix state ->
    TransitionMatrix state
mulTransitionMatrix = (<>)

{- | Return the @n x n@ identity: the zero-step transition that leaves every
state unchanged. For @n = 0@ this is the empty matrix.

Complexity: @O(1)@ construction time and @O(1)@ construction space. Forcing
the dense entries or support graph takes @O(n^2)@ time and @O(n^2)@ temporary
space; the support graph itself occupies @O(n)@ space.
-}
identityMatrix :: (FiniteState state) => TransitionMatrix state
identityMatrix = mempty

{- | Compute the @k@-step transition matrix @p^k@. Exponent zero returns
'identityMatrix'; positive exponents use repeated squaring through
'Data.Semigroup.mtimesDefault'.

Chapman-Kolmogorov gives @p^(m+n) = p^m p^n@ mathematically; computed matrices
may differ by floating-point rounding and are not revalidated.

Complexity: @O(n^2 + n^3 log(k + 1))@ time and @O(n^2)@ temporary and result
space.
-}
matrixPower ::
    (FiniteState state) =>
    Natural ->
    TransitionMatrix state ->
    TransitionMatrix state
matrixPower = mtimesDefault

{- | Return the stored row for a state: its next-state distribution.
'FiniteState' indexing makes the lookup total. The row is wrapped without
revalidation, so any floating-point drift from matrix arithmetic is
preserved.

Complexity: excluding 'Dtmc.State.stateIndex', @O(n^2)@ worst-case time and
@O(n^2)@ result space for state cardinality @n@.
-}
rowAt ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    DistributionVector state
rowAt = matrixRowAt
