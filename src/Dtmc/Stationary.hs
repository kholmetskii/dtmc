{- |
Module      : Dtmc.Stationary
Description : Stationary distributions of finite irreducible chains.

Stationary analysis for a finite transition matrix certified by an
'Irreducible' witness. Every finite irreducible DTMC has exactly one stationary
distribution, including periodic chains. For a row-stochastic matrix @P@, the
result @pi@ satisfies @transpose(P) pi = pi@ and @sum pi = 1@.

The calculation uses ordinary 'Double' arithmetic. It replaces one dependent
balance equation with the normalization equation and performs a checked LU
solve. The result is not clamped or renormalised. Non-finite, singular,
ill-conditioned, and high-residual systems return an explicit
'LinearSystemError'.
-}
module Dtmc.Stationary (
    LinearSystemError (..),
    stationaryDistribution,
) where

import Dtmc.Classification.Internal ( Irreducible, unIrreducible )
import Dtmc.Distribution.Vector.Internal (
    DistributionVector (DistributionVector),
 )
import Dtmc.Internal.LinearSystem (
    LinearSystemError (..),
    solveLinearSystem,
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.Transition.Matrix (
    unTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

{- | Compute the unique stationary distribution of a certified finite
irreducible chain.

The witness cannot exist for an empty chain, so the normalization row always
replaces one balance equation. The returned coordinates follow the canonical
state order of the 'FiniteState' instance.

Time: @O(n^3)@. Space: @O(n^2)@.
-}
stationaryDistribution ::
    (FiniteState state) =>
    Irreducible state ->
    Either LinearSystemError (DistributionVector state)
stationaryDistribution witness = do
    solution <- solveLinearSystem coefficient rightHandSide
    pure (DistributionVector (S.vector (LA.toList (LA.flatten solution))))
  where
    transition =
        S.extract (unTransitionMatrix (unIrreducible witness))
    dimension = LA.rows transition
    balanceRows =
        LA.toLists (LA.tr transition - LA.ident dimension)
    coefficient =
        LA.fromLists
            ( take (dimension - 1) balanceRows
                ++ [replicate dimension 1]
            )
    rightHandSide =
        LA.fromLists
            ( replicate (dimension - 1) [0]
                ++ [[1]]
            )
