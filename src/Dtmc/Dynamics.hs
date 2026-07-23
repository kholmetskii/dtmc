-- |
-- Module      : Dtmc.Dynamics
-- Description : Deterministic forward evolution of distributions.
--
-- The analytic (non-random) push-forward of a distribution by the chain: how
-- the law of the current state becomes the law of the next. Row-stochastic @p@
-- acts on a column distribution @mu@ by @mu' = transpose p #> mu@, i.e.
-- @mu'(j) = sum_i mu(i) * p(i,j)@. The @k@-step version reuses 'matrixPower'
-- from "Dtmc.TransitionMatrix".
module Dtmc.Dynamics (
    evolve,
    evolveN,
) where

import Dtmc.Internal.Types (
    Distribution (..),
    TransitionMatrix (..),
 )
import Dtmc.TransitionMatrix (
    matrixPower,
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
