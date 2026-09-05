{-# LANGUAGE TypeFamilies #-}

{- |
Module      : Dtmc.Analysis.Initial.Internal
Description : Mixing state-conditioned analysis results under an initial law.

Private helpers for turning quantities conditioned on @X_0 = i@ into the
corresponding quantity under an arbitrary finite-support initial distribution.
Weights and query results are combined with ordinary 'Double' arithmetic;
these helpers do not validate, clamp, or renormalise them.
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

{- | Compute the initial-law mixture of state-conditioned probabilities.
Every stored state is queried, and an empty weight list sums to zero.

Complexity: excluding 'distributionWeights' and query evaluations, @O(s)@
time, @O(s)@ temporary space, and @O(1)@ result space for @s@ stored weights.
-}
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

{- | Compute the initial-law mixture of fallible state-conditioned
probabilities. Queries are evaluated in stored order, and the first 'Left'
is returned without evaluating later queries. On success, values are combined
with ordinary 'Double' arithmetic.

Complexity: excluding 'distributionWeights' and query evaluations, @O(s)@
worst-case time, @O(s)@ temporary space, and @O(1)@ result space for @s@
stored weights.
-}
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

{- | Compute the initial-law mixture of fallible non-negative expectations.
A positive-weight 'InfiniteExpectation' makes the result infinite; a
zero-weight infinity is ignored. After the result becomes infinite, remaining
weights are traversed without evaluating their queries. Before that point,
the first query error is returned.

Finite values use ordinary 'Double' arithmetic without validation.

Complexity: excluding 'distributionWeights' and query evaluations, @O(s)@
worst-case time, @O(s)@ temporary space, and @O(1)@ result space for @s@
stored weights.
-}
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
