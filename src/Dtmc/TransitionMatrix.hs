{- |
Module      : Dtmc.TransitionMatrix
Description : Row-stochastic transition matrices and their monoid.

One-step transition probabilities for a DTMC on @n@ states.
'mkTransitionMatrix' validates each row with the simplex tolerance;
'mulTransitionMatrix', 'identityMatrix', and 'matrixPower' provide
multi-step transitions.
-}
module Dtmc.TransitionMatrix (
    TransitionMatrix,
    TransitionMatrixError (..),
    mkTransitionMatrix,
    unTransitionMatrix,
    mulTransitionMatrix,
    identityMatrix,
    matrixPower,
    rowAt,
    transitionProbability,
    transitionProbabilityN,
) where

import Data.Bifunctor (
    first,
 )
import Data.Finite (
    Finite,
    getFinite,
 )
import Data.Foldable (
    traverse_,
 )
import Data.Semigroup (
    mtimesDefault,
 )
import Dtmc.Distribution (
    probabilityAt,
 )
import Dtmc.Distribution.Internal (
    Distribution (Distribution),
 )
import Dtmc.Simplex (
    SimplexError,
 )
import Dtmc.Simplex.Internal (
    validateSimplex,
 )
import Dtmc.TransitionMatrix.Internal (
    TransitionMatrix,
    unTransitionMatrix,
    unsafeTransitionMatrix,
 )
import GHC.TypeNats (
    KnownNat,
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
mkTransitionMatrix :: (KnownNat n) => S.Sq n -> Either TransitionMatrixError (TransitionMatrix n)
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
mulTransitionMatrix :: (KnownNat n) => TransitionMatrix n -> TransitionMatrix n -> TransitionMatrix n
mulTransitionMatrix = (<>)

{- | The @n x n@ identity: the zero-step transition that leaves every state
unchanged. For @n = 0@ this is the empty matrix.

Time and result space: @O(n^2)@.
-}
identityMatrix :: (KnownNat n) => TransitionMatrix n
identityMatrix = mempty

{- | The @k@-step transition matrix @p^k@. Exponent zero returns
'identityMatrix'; positive exponents use repeated squaring through
'Data.Semigroup.mtimesDefault'.

Chapman-Kolmogorov gives @p^(m+n) = p^m p^n@ mathematically; computed matrices
may differ by floating-point rounding and are not revalidated.

Time: @O(n^2 + n^3 log(k + 1))@.
-}
matrixPower :: (KnownNat n) => Natural -> TransitionMatrix n -> TransitionMatrix n
matrixPower = mtimesDefault

{- | The stored row for state @i@: its next-state distribution. The 'Finite'
index makes the lookup total. The row is wrapped without revalidation, so any
floating-point drift from matrix arithmetic is preserved.
-}
rowAt :: (KnownNat n) => TransitionMatrix n -> Finite n -> Distribution n
rowAt p index =
    Distribution (S.toRows (unTransitionMatrix p) !! fromIntegral (getFinite index))

{- | The one-step transition probability @P(i, j) = P(X_1 = j | X_0 = i)@ from
source @i@ to destination @j@. By construction it reads the destination
coordinate of the source row:

@
transitionProbability p i j == 'probabilityAt' ('rowAt' p i) j
@

Both 'Finite' indices are bounded, so the query is total. The stored entry is
returned exactly, with no clamping, renormalisation, or revalidation.

Time: @O(n)@, dominated by extracting the source row.
-}
transitionProbability ::
    (KnownNat n) =>
    TransitionMatrix n ->
    -- | source @i@
    Finite n ->
    -- | destination @j@
    Finite n ->
    Double
transitionProbability p i = probabilityAt (rowAt p i)

{- | The @k@-step transition probability
@P^k(i, j) = P(X_k = j | X_0 = i)@, read from the @k@-th matrix power:

@
transitionProbabilityN k p i j == transitionProbability ('matrixPower' k p) i j
@

At @k = 0@ the power is the identity, so the result is the Kronecker delta:
one when @i == j@ and zero otherwise. Chapman-Kolmogorov holds in exact
arithmetic; the computed 'Double' is neither revalidated nor clamped, so
rounding from the repeated-squaring power is preserved.

Time: @O(n^2 + n^3 log(k + 1))@, dominated by 'matrixPower'.
-}
transitionProbabilityN ::
    (KnownNat n) =>
    Natural ->
    TransitionMatrix n ->
    -- | source @i@
    Finite n ->
    -- | destination @j@
    Finite n ->
    Double
transitionProbabilityN k p  = transitionProbability (matrixPower k p) 
