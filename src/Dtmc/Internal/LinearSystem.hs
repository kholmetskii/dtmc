{- |
Module      : Dtmc.Internal.LinearSystem
Description : @(I - Q)@ linear solves and sub-block extraction.

Dynamic linear algebra for DTMC hitting and return calculations: extract
blocks indexed by runtime state sets and solve systems of the form
@(I - Q) x = b@. Public modules convert bounded state indices to 'Int' and
keep dynamically sized matrices inside the implementation.

All arithmetic uses 'Double'. Every solve validates finiteness, rejects a
reciprocal condition estimate below @1e-12@, and verifies a scaled residual
against @1e-9@.
-}
module Dtmc.Internal.LinearSystem (
    LinearSystemError (..),
    subMatrix,
    rowSums,
    solveIminusQ,
    solveIminusQVector,
    fundamental,
) where

import Numeric.LinearAlgebra qualified as LA

{- | A numerical linear system could not be accepted safely.

The reciprocal condition estimate and relative residual are dimensionless.
Smaller reciprocal condition estimates indicate greater sensitivity; smaller
relative residuals indicate a better computed solution.
-}
data LinearSystemError
    = -- | A solver reported that the coefficient matrix is singular.
      SingularSystem
    | -- | The coefficient matrix is too sensitive for the numerical contract.
      IllConditionedSystem
        { reciprocalConditionEstimate :: Double
        }
    | -- | A coefficient or right-hand-side entry was @NaN@ or infinite.
      NonFiniteSystem
    | -- | The solver produced a @NaN@ or infinite result.
      NonFiniteSolution
    | -- | The computed solution did not satisfy the equations closely enough.
      ResidualTooLarge
        { relativeResidual :: Double
        , residualLimit :: Double
        }
    deriving (Eq, Show)

conditionLimit :: Double
conditionLimit = 1e-12

residualLimitValue :: Double
residualLimitValue = 1e-9

allFinite :: LA.Matrix Double -> Bool
allFinite = all isFinite . LA.toList . LA.flatten

isFinite :: Double -> Bool
isFinite value = not (isNaN value || isInfinite value)

infinityNorm :: LA.Matrix Double -> Double
infinityNorm = maximum . (0 :) . map (sum . map abs) . LA.toLists

relativeSystemResidual ::
    LA.Matrix Double ->
    LA.Matrix Double ->
    LA.Matrix Double ->
    Double
relativeSystemResidual coefficient solution rightHandSide =
    infinityNorm (coefficient LA.<> solution - rightHandSide)
        / max
            1
            ( infinityNorm coefficient * infinityNorm solution
                + infinityNorm rightHandSide
            )

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

Returns a 'LinearSystemError' when the input or result is non-finite, the backend
reports singularity, the reciprocal condition estimate is below @1e-12@, or
the scaled infinity-norm residual exceeds @1e-9@. Current DTMC callers choose
@Q@ with spectral radius below one, which guarantees invertibility in exact
arithmetic, but does not guarantee a reliable 'Double' result.

Time: @O(n^3 + n^2 r)@. Space: @O(n^2 + nr)@.
-}
solveIminusQ ::
    LA.Matrix Double ->
    LA.Matrix Double ->
    Either LinearSystemError (LA.Matrix Double)
solveIminusQ q rightHandSide
    | not (allFinite q && allFinite rightHandSide) = Left NonFiniteSystem
    | otherwise =
        case LA.linearSolve coefficient rightHandSide of
            Nothing -> Left SingularSystem
            Just solution
                | not (allFinite solution) -> Left NonFiniteSolution
                | not (isFinite reciprocalCondition) -> Left NonFiniteSystem
                | reciprocalCondition < conditionLimit ->
                    Left (IllConditionedSystem reciprocalCondition)
                | residual > residualLimitValue ->
                    Left
                        ( ResidualTooLarge
                            { relativeResidual = residual
                            , residualLimit = residualLimitValue
                            }
                        )
                | otherwise -> Right solution
              where
                residual =
                    relativeSystemResidual coefficient solution rightHandSide
  where
    coefficient = LA.ident (LA.rows q) - q
    reciprocalCondition = LA.rcond coefficient

{- | 'solveIminusQ' for a vector of length @n@. It has the same shape,
singularity, and floating-point behaviour.

Time: @O(n^3)@. Space: @O(n^2)@.
-}
solveIminusQVector ::
    LA.Matrix Double ->
    LA.Vector Double ->
    Either LinearSystemError (LA.Vector Double)
solveIminusQVector q b =
    LA.flatten <$> solveIminusQ q (LA.asColumn b)

{- | Compute @(I - Q)^-1@ by solving @(I - Q) G = I@. When @Q@ is a
transient-to-transient transition block with spectral radius below one, this
is the fundamental matrix @sum_{k=0}^infinity Q^k@.

Requires a non-empty square matrix with finite entries and inherits the
validation and error behaviour of 'solveIminusQ'.

Time: @O(n^3)@. Space: @O(n^2)@.
-}
fundamental :: LA.Matrix Double -> Either LinearSystemError (LA.Matrix Double)
fundamental q =
    solveIminusQ q (LA.ident (LA.rows q))
