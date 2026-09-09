{- |
Module      : Dtmc.Analysis.Classification.Reachability.Internal
Description : Solver-oriented typed support-graph reachability.

Internal adapter between finite-state values and the integer support graph.
-}
module Dtmc.Analysis.Classification.Reachability.Internal (
    backwardReachable,
) where

import Data.Maybe (fromMaybe)
import Dtmc.State (FiniteState)
import Dtmc.State.Internal (stateFromInt, stateIndexInt)
import Dtmc.Transition.Matrix.Internal (TransitionMatrix, tmSupport)
import Dtmc.Transition.Matrix.Internal.Graph qualified as G

toState :: (FiniteState state) => Int -> state
toState index =
    fromMaybe
        (error "Dtmc.Analysis.Classification: graph vertex out of bounds")
        (stateFromInt index)

{- | Return states from which an allowed seed is reachable along a support
path containing only states accepted by @allowed@. Disallowed seeds are
ignored; the result is duplicate-free and ordered by state index.

For the complexity bounds, @n@ is the state count, @E@ the support-edge count,
@s@ the number of supplied seeds, and @r@ the number of returned states.

Complexity: excluding @n@ evaluations of @allowed@, 'FiniteState' method
costs, and shared support-graph construction, @O(n + E + s)@ time,
@O(n + E + s)@ temporary space, and @O(r)@ result space. The first reverse
traversal also retains @O(n + E)@ predecessor-cache space.
-}
backwardReachable ::
    (FiniteState state) =>
    TransitionMatrix state ->
    (state -> Bool) ->
    [state] ->
    [state]
backwardReachable matrix allowed seeds =
    map toState
        ( G.backwardReachable
            (tmSupport matrix)
            (allowed . toState)
            (map stateIndexInt seeds)
        )
