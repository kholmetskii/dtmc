{- |
Module      : Dtmc.Analysis.LinearSystem.Internal
Description : @(I - Q)@ linear solves and sub-block extraction.

Dynamic linear algebra for DTMC hitting and return calculations: extract
blocks indexed by runtime state sets and solve systems of the form
@(I - Q) x = b@. Public modules convert bounded state indices to 'Int' and
keep dynamically sized matrices inside the implementation.

All arithmetic uses 'Double'. Every solve validates finiteness, rejects a
reciprocal condition estimate below @1e-12@, and verifies a scaled residual
against @1e-9@.
-}
module Dtmc.Analysis.LinearSystem.Internal (
    subMatrix,
    rowSums,
    solveLinearSystem,
    solveIminusQ,
    solveIminusQVector,
    fundamental,
) where

import Data.Vector.Storable qualified as Storable
import Dtmc.Analysis.LinearSystem (
    LinearSystemError (..),
 )
import Numeric.LinearAlgebra qualified as LA

conditionLimit :: Double
conditionLimit = 1e-12

residualLimitValue :: Double
residualLimitValue = 1e-9

allFinite :: LA.Matrix Double -> Bool
allFinite = Storable.all isFinite . LA.flatten

isFinite :: Double -> Bool
isFinite value = not (isNaN value || isInfinite value)

infinityNorm :: LA.Matrix Double -> Double
infinityNorm = maximumAbsoluteSum . LA.toRows

oneNorm :: LA.Matrix Double -> Double
oneNorm = maximumAbsoluteSum . LA.toColumns

maximumAbsoluteSum :: [LA.Vector Double] -> Double
maximumAbsoluteSum =
    foldl' max 0
        . map (Storable.foldl' (\acc value -> acc + abs value) 0)

maximumAbsolute :: LA.Vector Double -> Double
maximumAbsolute = Storable.foldl' (\acc value -> max acc (abs value)) 0

reciprocalConditionForIminusQ ::
    LA.Matrix Double ->
    LA.Vector Double ->
    Double
reciprocalConditionForIminusQ coefficient inverseTimesOne
    | Storable.all isFinite inverseTimesOne
        && isFinite provenLowerBound
        && provenLowerBound >= conditionLimit =
        provenLowerBound
    | otherwise = LA.rcond coefficient
  where
    inverseInfinityNorm = maximumAbsolute inverseTimesOne
    -- The inverse one-norm is at most n times its infinity norm, and every
    -- spectral norm is at most the square root of one-norm times infinity norm.
    provenLowerBound =
        1
            / ( inverseInfinityNorm
                    * sqrt
                        ( fromIntegral (LA.rows coefficient)
                            * oneNorm coefficient
                            * infinityNorm coefficient
                        )
              )

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

validateSolution ::
    LA.Matrix Double ->
    LA.Matrix Double ->
    LA.Matrix Double ->
    Double ->
    Either LinearSystemError (LA.Matrix Double)
validateSolution coefficient rightHandSide solution reciprocalCondition
    | not (allFinite solution) = Left NonFiniteSolution
    | not (isFinite reciprocalCondition) = Left NonFiniteSystem
    | reciprocalCondition < conditionLimit =
        Left (IllConditionedSystem reciprocalCondition)
    | residual > residualLimitValue =
        Left
            ( ResidualTooLarge
                { relativeResidual = residual
                , residualLimit = residualLimitValue
                }
            )
    | otherwise = Right solution
  where
    residual = relativeSystemResidual coefficient solution rightHandSide

{- | Solve @A X = B@ by LU decomposition. @A@ must be a non-empty square
matrix, @B@ must have the same number of rows, and all entries must be finite;
incompatible dimensions raise a backend error.

Returns a 'LinearSystemError' when the input or result is non-finite, the
backend reports singularity, the reciprocal condition estimate is below
@1e-12@, or the scaled infinity-norm residual exceeds @1e-9@.

Complexity: @O(n^3 + n^2 r)@ time, @O(n^2 + n r)@ temporary space, and
@O(n r)@ result space for an @n x n@ coefficient matrix and @r@ right-hand
sides.
-}
solveLinearSystem ::
    LA.Matrix Double ->
    LA.Matrix Double ->
    Either LinearSystemError (LA.Matrix Double)
solveLinearSystem coefficient rightHandSide
    | not (allFinite coefficient && allFinite rightHandSide) =
        Left NonFiniteSystem
    | otherwise =
        case LA.linearSolve coefficient rightHandSide of
            Nothing -> Left SingularSystem
            Just solution ->
                validateSolution
                    coefficient
                    rightHandSide
                    solution
                    reciprocalCondition
  where
    reciprocalCondition = LA.rcond coefficient

{- | Extract the block of @m@ selected by the row and column indices. Their
order and multiplicity are preserved; an empty list produces a zero-sized
dimension.

Row indices must be in @{0 .. rows(m)-1}@ and column indices in
@{0 .. cols(m)-1}@; otherwise the backend raises an error.

Complexity: @O(R + C + R C)@ time, @O(R + C)@ temporary space, and
@O(R C)@ result space for @R@ selected rows and @C@ selected columns.
-}
subMatrix :: [Int] -> [Int] -> LA.Matrix Double -> LA.Matrix Double
subMatrix rowIdx colIdx m =
    m LA.?? (LA.Pos (LA.idxs rowIdx), LA.Pos (LA.idxs colIdx))

{- | Compute the vector of row sums of @m@, equivalently @m@ applied to a
vector of ones. For a block @P[D, A]@ of a transition matrix, this is the
one-step probability of transitioning from each state in @D@ directly into
@A@. Results use ordinary floating-point summation and are not clamped to
@[0, 1]@.

Complexity: @O(R C + R + C)@ time, @O(C)@ temporary space, and @O(R)@
result space for an @R x C@ matrix.
-}
rowSums :: LA.Matrix Double -> LA.Vector Double
rowSums m = m LA.#> LA.konst 1 (LA.cols m)

{- | Solve @(I - Q) X = B@ by one LU decomposition. @Q@ must be a
transient-to-transient, substochastic @n x n@ matrix, @B@ must be @n x r@,
and all entries must be finite; incompatible dimensions raise a backend
error.

Returns a 'LinearSystemError' when the input or result is non-finite, the
backend reports singularity, the reciprocal condition estimate is below
@1e-12@, or the scaled infinity-norm residual exceeds @1e-9@.

The solve appends a column of ones to @B@. Since @(I - Q)^-1@ is
non-negative, the extra result gives its exact infinity norm in exact
arithmetic. Matrix-norm inequalities then provide a conservative lower bound
on the reciprocal 2-norm condition number. A bound above the acceptance limit
proves that the solution is sufficiently conditioned without a separate
singular-value decomposition; uncertain cases retain the backend condition
estimate and therefore the same rejection threshold.

Complexity: @O(n^3 + n^2 r)@ time, @O(n^2 + n r)@ temporary space, and
@O(n r)@ result space for @r@ right-hand sides.
-}
solveIminusQ ::
    LA.Matrix Double ->
    LA.Matrix Double ->
    Either LinearSystemError (LA.Matrix Double)
solveIminusQ q rightHandSide
    | not (allFinite q && allFinite rightHandSide) = Left NonFiniteSystem
    | otherwise =
        case LA.linearSolve coefficient augmentedRightHandSide of
            Nothing -> Left SingularSystem
            Just augmentedSolution ->
                let solution =
                        LA.takeColumns
                            rightHandSideColumns
                            augmentedSolution
                    inverseTimesOne =
                        LA.flatten
                            ( LA.dropColumns
                                rightHandSideColumns
                                augmentedSolution
                            )
                 in validateSolution
                        coefficient
                        rightHandSide
                        solution
                        ( reciprocalConditionForIminusQ
                            coefficient
                            inverseTimesOne
                        )
  where
    coefficient = LA.ident (LA.rows q) - q
    rightHandSideColumns = LA.cols rightHandSide
    augmentedRightHandSide =
        rightHandSide LA.||| LA.asColumn (LA.konst 1 (LA.rows q))

{- | Solve @(I - Q) x = b@ for a vector of length @n@. This has the same
shape requirements, numerical behaviour, and errors as 'solveIminusQ'.

Complexity: @O(n^3)@ time, @O(n^2)@ temporary space, and @O(n)@ result
space.
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

This requires a non-empty square matrix with finite entries and inherits the
validation and error behaviour of 'solveIminusQ'.

Complexity: @O(n^3)@ time and @O(n^2)@ temporary and result space.
-}
fundamental :: LA.Matrix Double -> Either LinearSystemError (LA.Matrix Double)
fundamental q =
    solveIminusQ q (LA.ident (LA.rows q))
