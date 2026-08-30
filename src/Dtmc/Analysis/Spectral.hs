{- |
Module      : Dtmc.Analysis.Spectral
Description : Eigenvalues of a transition matrix and the convergence rate.

The spectrum of a finite transition matrix, and the two scalars derived from
it that govern how fast @P^n@ settles.

Every stochastic matrix has @1@ as an eigenvalue and no eigenvalue of modulus
greater than @1@. For an irreducible chain the eigenvalues of modulus @1@ are
exactly the @d@-th roots of unity, where @d@ is the period, so such a chain is
aperiodic precisely when @1@ is the only one. The second largest modulus then
bounds the geometric rate at which @P^n@ approaches the limit computed by
"Dtmc.Analysis.Limiting".

Everything here is a floating-point eigenvalue computation and is approximate,
in contrast to the exact combinatorial period and classification in
"Dtmc.Analysis.Classification". Prefer those where a structural answer is
wanted: 'Dtmc.Analysis.Classification.aperiodic' decides periodicity without
arithmetic, while a spectral test would need a tolerance.

There is deliberately no diagonalisability predicate. Diagonalisability is a
measure-zero property -- an arbitrarily small perturbation makes a defective
matrix diagonalisable -- so no floating-point test can decide it. The
diagonalisation formula @(P^n)(i,j) = sum_k Q(i,k) Q^-1(k,j) lambda_k^n@ holds
under that hypothesis but is not offered as a function: 'matrixPower' computes
the same matrix by repeated squaring, faster and without an eigenvector solve.
-}
module Dtmc.Analysis.Spectral (
    spectrum,
    secondLargestModulus,
    spectralGap,
) where

import Data.Complex (
    Complex,
    magnitude,
 )
import Data.List (
    sortOn,
 )
import Data.Ord (
    Down (Down),
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    unTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

{- | The eigenvalues of the transition matrix, ordered by decreasing modulus.
An @n@-state chain has @n@ of them, listed with algebraic multiplicity, and
they may be genuinely complex: the three-state chain of section 4.2 has
@(-5 +- i sqrt 5) / 10@.

Two facts hold for every stochastic matrix. Since each row sums to one,
@P 1 = 1@, so @1@ is an eigenvalue and the first element of the result. And
every eigenvalue satisfies @|lambda| <= 1@: if @P x = lambda x@, pick a
coordinate @i@ of largest modulus, then

@|lambda| |x_i| = |sum_j P(i,j) x_j| <= sum_j P(i,j) |x_j| <= |x_i|@.

The computation is a LAPACK eigenvalue decomposition and is approximate. The
empty chain has an empty spectrum.

Time: @O(n^3)@.
-}
spectrum ::
    (FiniteState state) =>
    TransitionMatrix state ->
    [Complex Double]
spectrum p =
    sortOn (Down . magnitude) (LA.toList (LA.eigenvalues matrix))
  where
    matrix = S.extract (unTransitionMatrix p)

{- | The largest modulus among the eigenvalues once one copy of the Perron
root has been dropped -- the /second largest eigenvalue modulus/. 'Nothing'
for a chain with fewer than two states, which has no second eigenvalue.

For an irreducible aperiodic chain this is strictly below one and
@|(P^n)(i,j) - pi(j)|@ decays like its @n@-th power, so it is the rate at
which the limit of "Dtmc.Analysis.Limiting" is approached. A value of one
signals that no such uniform rate exists, which happens in two ways: an
irreducible chain of period @d > 1@ carries @d@ eigenvalues on the unit
circle, and a chain with several recurrent classes carries one eigenvalue
@1@ per class.

Being a computed modulus, this is approximate and near-unit values should not
be compared against one without a tolerance. To decide periodicity exactly,
use 'Dtmc.Analysis.Classification.aperiodic'.

Time: @O(n^3)@.
-}
secondLargestModulus ::
    (FiniteState state) =>
    TransitionMatrix state ->
    Maybe Double
secondLargestModulus p =
    case spectrum p of
        (_ : second : _) -> Just (magnitude second)
        _ -> Nothing

{- | @1@ minus 'secondLargestModulus': the spectral gap. Larger means faster
convergence of @P^n@; a gap of zero means no geometric rate, and a small
negative value is floating-point noise around a modulus of one rather than a
meaningful quantity.

'Nothing' exactly when 'secondLargestModulus' is.
-}
spectralGap ::
    (FiniteState state) =>
    TransitionMatrix state ->
    Maybe Double
spectralGap p =
    (1 -) <$> secondLargestModulus p
