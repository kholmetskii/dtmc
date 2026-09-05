{- |
Module      : Dtmc.Dynamics.Internal
Description : Unchecked sparse forward-dynamics primitive.

Sparse weight propagation shared by public evolution and path-analysis
algorithms. Inputs may be sub-probability maps; no simplex invariant is
required or restored here.
-}
module Dtmc.Dynamics.Internal (
    pushSparseWeights,
) where

import Data.Map.Strict (
    Map,
 )
import Data.Map.Strict qualified as Map
import Dtmc.Distribution.Map.Internal (
    unDistributionMap,
 )
import Dtmc.Transition (
    Transition (..),
 )

{- | Push a finite, possibly sub-probability weight map through one locally
finite kernel step. Exact zero results are removed. No validation, clamping,
or renormalisation is performed.

For the complexity bounds, @s@ is the number of source states, @e@ the number
of traversed support edges, @u@ the number of distinct destinations
encountered, and @r@ the number retained after exact-zero removal.

Complexity: excluding 'transitionLaw' evaluation,
@O(s + e log(u + 1) + u)@ time, @O(u)@ temporary space, and @O(r)@ result
space.
-}
pushSparseWeights ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    Map (TransitionState kernel) Double ->
    kernel ->
    Map (TransitionState kernel) Double
pushSparseWeights weights kernel =
    Map.filter (/= 0) (Map.foldlWithKey' pushState Map.empty weights)
  where
    pushState accumulated state stateWeight =
        Map.foldlWithKey'
            ( \next nextState transitionWeight ->
                Map.insertWith
                    (+)
                    nextState
                    (stateWeight * transitionWeight)
                    next
            )
            accumulated
            (unDistributionMap (transitionLaw kernel state))
