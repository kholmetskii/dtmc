{-# LANGUAGE TypeFamilies #-}

{- |
Module      : Dtmc.Analysis.Initial.Internal
Description : Mixing state-conditioned analysis results under an initial law.

Private helpers for turning quantities conditioned on @X_0 = i@ into the
corresponding quantity under an arbitrary finite-support initial distribution.
-}
module Dtmc.Analysis.Initial.Internal (
    probabilityUnder,
    probabilityUnderEither,
    expectationUnderEither,
) where

import Control.Monad (
    foldM,
 )
import Dtmc.Analysis.Expectation (
    Expectation (..),
 )
import Dtmc.Distribution (
    Distribution (..),
 )

probabilityUnder ::
    (Distribution distribution) =>
    distribution ->
    (DistributionState distribution -> Double) ->
    Double
probabilityUnder initial query =
    sum
        [ weight * query state
        | (state, weight) <- distributionWeights initial
        ]

probabilityUnderEither ::
    (Distribution distribution) =>
    distribution ->
    (DistributionState distribution -> Either error Double) ->
    Either error Double
probabilityUnderEither initial query =
    sum
        <$> traverse
            (\(state, weight) -> (weight *) <$> query state)
            (distributionWeights initial)

expectationUnderEither ::
    (Distribution distribution) =>
    distribution ->
    (DistributionState distribution -> Either error Expectation) ->
    Either error Expectation
expectationUnderEither initial query =
    foldM add (FiniteExpectation 0) (distributionWeights initial)
  where
    add InfiniteExpectation _ = Right InfiniteExpectation
    add (FiniteExpectation total) (state, weight) = do
        value <- query state
        pure $
            case value of
                InfiniteExpectation
                    | weight > 0 -> InfiniteExpectation
                    | otherwise -> FiniteExpectation total
                FiniteExpectation x ->
                    FiniteExpectation (total + weight * x)
