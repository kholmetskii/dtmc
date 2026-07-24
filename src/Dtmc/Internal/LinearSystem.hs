{- |
Module      : Dtmc.Internal.LinearSystem
Description : @(I - Q)@ linear solves and sub-block extraction.

Dynamic linear algebra for DTMC hitting and return calculations: extract
blocks indexed by runtime state sets and solve systems of the form
@(I - Q) x = b@. Public modules convert bounded state indices to 'Int' and
keep dynamically sized matrices inside the implementation.

All arithmetic uses 'Double'. This module adds no numeric tolerance,
conditioning test, or residual check.
-}
module Dtmc.Internal.LinearSystem (
    subMatrix,
    rowSums,
    solveIminusQ,
    solveIminusQVector,
    fundamental,
) where

import Numeric.LinearAlgebra qualified as LA

{- | The block of @m@ selected by the row and column indices. Their order and
multiplicity are preserved; an empty list produces a zero-sized dimension.

Row indices must be in @{0 .. rows(m)-1}@ and column indices in
@{0 .. cols(m)-1}@; otherwise the backend raises an error.

Time and result space: @O(RC)@ for @R@ selected rows and @C@ selected columns.
-}
subMatrix :: [Int] -> [Int] -> LA.Matrix Double -> LA.Matrix Double
subMatrix rowIdx colIdx m =
    m LA.?? (LA.Pos (LA.idxs rowIdx), LA.Pos (LA.idxs colIdx))

{- | The vector of row sums of @m@, i.e. @m@ applied to a vector of ones. For
a block @P[D, A]@ of a transition matrix this is the one-step probability
of transitioning from each state in @D@ directly into @A@. Results use
ordinary floating-point summation and are not clamped to @[0, 1]@.

Time: @O(RC)@. Result space: @O(R)@.
-}
rowSums :: LA.Matrix Double -> LA.Vector Double
rowSums m = m LA.#> LA.konst 1 (LA.cols m)

{- | Solve @(I - Q) X = B@ by LU decomposition. @Q@ must be @n x n@, @B@
must be @n x r@, and all entries must be finite; incompatible dimensions
raise a backend error.

Returns 'Nothing' when the backend reports singularity, including @n = 0@.
Current DTMC callers choose @Q@ with spectral radius below one, which
guarantees invertibility in exact arithmetic. This wrapper applies no
tolerance or residual check, so an ill-conditioned system may yield an
inaccurate result.

Time: @O(n^3 + n^2 r)@. Space: @O(n^2 + nr)@.
-}
solveIminusQ :: LA.Matrix Double -> LA.Matrix Double -> Maybe (LA.Matrix Double)
solveIminusQ q =
    LA.linearSolve (LA.ident (LA.rows q) - q)

{- | 'solveIminusQ' for a vector of length @n@. It has the same shape,
singularity, and floating-point behaviour.

Time: @O(n^3)@. Space: @O(n^2)@.
-}
solveIminusQVector :: LA.Matrix Double -> LA.Vector Double -> Maybe (LA.Vector Double)
solveIminusQVector q b =
    LA.flatten <$> solveIminusQ q (LA.asColumn b)

{- | Compute @(I - Q)^-1@ by solving @(I - Q) G = I@. When @Q@ is a
transient-to-transient transition block with spectral radius below one, this
is the fundamental matrix @sum_{k=0}^infinity Q^k@.

Requires a non-empty square matrix with finite entries and inherits the
singularity and floating-point behaviour of 'solveIminusQ'.

Time: @O(n^3)@. Space: @O(n^2)@.
-}
fundamental :: LA.Matrix Double -> Maybe (LA.Matrix Double)
fundamental q =
    solveIminusQ q (LA.ident (LA.rows q))
