{- |
Module      : Dtmc.Dynamics.Internal
Description : Unchecked sparse forward-dynamics primitive.
-}
module Dtmc.Dynamics.Internal (
    pushSparseWeights,
) where

import Data.Map.Strict (
    Map,
 )
import Data.Map.Strict qualified as Map
import Dtmc.Distribution.Internal (
    unSparseDistribution,
 )
import Dtmc.Kernel (
    Transition (..),
 )

{- | Push a finite, possibly sub-probability weight map through one locally
finite kernel step. Exact zero results are removed. No validation, clamping,
or renormalisation is performed.
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
            (unSparseDistribution (transitionLaw kernel state))
