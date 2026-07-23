-- |
-- Module      : Dtmc.Dynamics
-- Description : Deterministic forward dynamics and multi-step transitions.
--
-- The analytic (non-random) evolution of a chain: how a distribution is pushed
-- forward by the transition matrix, and how one-step matrices compose into
-- @k@-step matrices. Row-stochastic @p@ acts on a column distribution @mu@ by
-- @mu' = transpose p #> mu@, i.e. @mu'(j) = sum_i mu(i) * p(i,j)@. Powers and
-- products reuse the 'Monoid' structure of t'TransitionMatrix'.
module Dtmc.Dynamics (
    evolve,
    evolveN,
    identityMatrix,
    matrixPower,
) where

import Data.Semigroup (
    mtimesDefault,
 )
import Dtmc.Internal.Types (
    Distribution (..),
    TransitionMatrix (..),
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (Natural)

-- | One-step push-forward of a distribution: @evolve mu p = transpose p #> mu@.
-- Maps the law of the current state to the law of the next state.
evolve :: (KnownNat n) => Distribution n -> TransitionMatrix n -> Distribution n
evolve (Distribution v) p =
    Distribution
        { unDistribution = S.tr (unTransitionMatrix p) S.#> v
        }

-- | The @n*n@ identity as a transition matrix: the zero-step transition
-- (@mempty@), which leaves any distribution unchanged.
identityMatrix :: (KnownNat n) => TransitionMatrix n
identityMatrix = mempty

-- | The @k@-step transition matrix @p^k@, formed by @k@-fold monoidal product
-- (@matrixPower 0 = identityMatrix@).
matrixPower :: (KnownNat n) => Natural -> TransitionMatrix n -> TransitionMatrix n
matrixPower = mtimesDefault

-- | @k@-step push-forward: @evolveN k mu p = evolve mu (p^k)@, the law of the
-- state after @k@ transitions.
evolveN ::
    (KnownNat n) =>
    Natural ->
    Distribution n ->
    TransitionMatrix n ->
    Distribution n
evolveN k mu p =
    evolve mu (matrixPower k p)
