{-# LANGUAGE ExplicitNamespaces #-}

{- |
Module      : Dtmc.Analysis.Classification
Description : Communication, irreducibility, periodicity, and recurrence.

Qualitative DTMC properties derived from the support graph of @P@: there is an
edge @i -> j@ exactly when the stored @P(i,j) > 0@. The comparison has no
tolerance, so a tiny positive value is a transition while zero or a negative
value is not. Results depend on which entries are positive, not their
magnitudes. Recurrence statements assume a finite, valid transition matrix.
Queries accept named state constructors through 'FiniteState'; state lists are
returned in the canonical order of that instance.

For the complexity bounds, @n@ is the number of states and @E@ is the number
of strictly positive entries. The stated per-operation bounds exclude
'FiniteState' method costs and construction of the shared support graph. Its
first use adds @O(n^2)@ time and temporary space and retains @O(n + E)@ cache
space. Strong components, closedness, periods, and phases are also computed
lazily and shared by later queries on the same matrix.

A @0 x 0@ matrix has no communicating classes and is neither irreducible nor
aperiodic here.
-}
module Dtmc.Analysis.Classification (
    -- * Reachability
    supportEdge,
    accessible,
    reachesAny,
    communicates,

    -- * Communicating classes
    communicatingClasses,
    irreducible,

    -- * Periodicity
    period,
    aperiodic,
    cyclicClasses,

    -- * Recurrence and transience
    recurrentState,
    transientState,
    recurrentStates,
    transientStates,

    -- * Classification summary
    type CommClass (..),
    type Classification,
    classesOf,
    isIrreducible,
    isAperiodic,
    isErgodic,
    chainPeriod,
    recurrentStatesOf,
    transientStatesOf,
    absorbingStates,
    classify,
) where

import Data.Array qualified as Array
import Data.Maybe (
    fromMaybe,
 )
import Dtmc.Analysis.Classification.Internal (
    type Classification (..),
    type CommClass (..),
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.State.Internal (
    stateFromInt,
    stateIndexInt,
 )
import Dtmc.Transition.Matrix.Internal (TransitionMatrix, tmSupport)
import Dtmc.Transition.Matrix.Internal.Graph qualified as G
import Numeric.Natural (Natural)

toState :: (FiniteState state) => Int -> state
toState index =
    fromMaybe
        (error "Dtmc.Analysis.Classification: graph vertex out of bounds")
        (stateFromInt index)

toIndex :: (FiniteState state) => state -> Int
toIndex = stateIndexInt

{- | Test whether @P(i,j) > 0@, so that the support graph contains the direct
edge @i -> j@. No tolerance is applied.

Complexity: excluding shared support-graph construction, @O(d_i)@ time and
@O(1)@ temporary and result space for out-degree @d_i@ of state @i@.
-}
supportEdge ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    state ->
    Bool
supportEdge p i j = G.hasEdge (tmSupport p) (toIndex i) (toIndex j)

{- | Test whether @j@ is reachable from @i@ in zero or more transitions.
Every state is therefore reachable from itself, even without a self-loop.

Complexity: excluding shared support-graph construction, @O(n + E)@ time,
@O(n)@ temporary space, and @O(1)@ result space.
-}
accessible ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    state ->
    Bool
accessible p i j = G.reachable (tmSupport p) (toIndex i) (toIndex j)

{- | Test whether any supplied target is reachable from @i@ in zero or more
transitions. An empty target list gives 'False'; including @i@ gives 'True'.

The graph is traversed once rather than once per target.

Complexity: excluding shared support-graph construction, @O(n + E + t)@
worst-case time, @O(n + t)@ temporary space, and @O(1)@ result space for @t@
supplied targets.
-}
reachesAny ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    [state] ->
    Bool
reachesAny p i targets =
    G.reachesAny (tmSupport p) (toIndex i) (map toIndex targets)

{- | Test whether @i@ and @j@ communicate, meaning that each is reachable
from the other. This is an equivalence relation on the state space.

Complexity: excluding shared support-graph construction, the first query
takes @O(n + E + n log(n + 1))@ time and @O(n + E)@ temporary space and
retains @O(n)@ component-cache space; subsequent queries take @O(1)@ time.
Temporary and result space per cached query are @O(1)@.
-}
communicates ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    state ->
    Bool
communicates p i j =
    G.sameComponent (tmSupport p) (toIndex i) (toIndex j)

{- | Return the communicating classes, equivalently the strongly connected
components of the support graph. States within each class are ascending, and
classes are ordered by their least member.

This is a focused projection of the complete 'classify' report.

Complexity: excluding shared support-graph construction, the first full
evaluation takes @O(n + E + n log(n + 1))@ time and @O(n + E)@ temporary
space and retains @O(n)@ component cache; subsequent evaluations take
@O(n)@ time and temporary space. Result space is @O(n)@.
-}
communicatingClasses :: (FiniteState state) => TransitionMatrix state -> [[state]]
communicatingClasses = map classMembers . classesOf . classify

{- | Test whether every state communicates with every other state. The empty
chain is not irreducible.

Complexity: excluding shared support-graph construction, the first query
takes @O(n + E + n log(n + 1))@ time and @O(n + E)@ temporary space and
retains @O(n)@ component-cache space; subsequent queries take @O(1)@ time.
Temporary and result space per cached query are @O(1)@.
-}
irreducible :: TransitionMatrix state -> Bool
irreducible = graphIrreducible . tmSupport

graphIrreducible :: G.Graph -> Bool
graphIrreducible graph =
    case G.components graph of
        [component] -> not (null component)
        _ -> False

{- | Return the period of @i@:
@gcd { k >= 1 | (P^k)(i,i) > 0 }@. Returns 'Nothing' when @i@ has no
positive-length return path, necessarily a singleton class without a
self-transition.

Complexity: excluding shared support-graph construction, the first query
takes @O((n + E) log(n + 1))@ time and @O(n + E)@ temporary space and retains
@O(n)@ period-cache space; subsequent queries take @O(1)@ time. Temporary
and result space per cached query are @O(1)@.
-}
period ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Maybe Natural
period p i = G.periodOf (tmSupport p) (toIndex i)

{- | Test whether every communicating class has period @1@. The empty chain
and a chain containing a class with undefined period are not aperiodic under
this definition.

Complexity: excluding shared support-graph construction, the first query
takes @O((n + E) log(n + 1))@ time and @O(n + E)@ temporary space and retains
@O(n)@ component and period cache; later queries take @O(c)@ time for @c@
communicating classes. Temporary and result space per cached query are
@O(1)@.
-}
aperiodic :: TransitionMatrix state -> Bool
aperiodic = graphAperiodic . tmSupport

graphAperiodic :: G.Graph -> Bool
graphAperiodic graph =
    not (null components)
        && all ((== Just 1) . G.componentPeriod graph) components
  where
    components = G.components graph

{- | Partition an irreducible chain of period @d@ into cyclic classes
@C_0, ..., C_(d-1)@. Every transition from @C_r@ enters
@C_((r+1) mod d)@; @C_0@ contains the least state. Returns 'Nothing' for a
reducible chain or an undefined period.

Complexity: excluding shared support-graph construction, the first full
evaluation takes @O((n + E) log(n + 1))@ time and @O(n + E)@ temporary space
and retains @O(n)@ component, period, and phase cache; later evaluations take
@O(n)@ time and temporary space. Result space is @O(n)@.
-}
cyclicClasses :: (FiniteState state) => TransitionMatrix state -> Maybe [[state]]
cyclicClasses p
    | not (irreducible p) = Nothing
    | otherwise =
        case G.periodOf g 0 of
            Nothing -> Nothing
            Just d ->
                let dInt = fromIntegral d
                    -- One pass buckets every vertex by its phase (@O(V + d)@),
                    -- rather than scanning all vertices once per phase.
                    buckets =
                        Array.accumArray
                            (flip (:))
                            []
                            (0, dInt - 1)
                            [(G.phaseOf g v, toState v) | v <- [0 .. G.graphDim g - 1]]
                 in Just [reverse (buckets Array.! r) | r <- [0 .. dInt - 1]]
  where
    g = tmSupport p

{- | Test whether the chain returns to @i@ with probability one when started
there. For a finite DTMC this holds exactly when @i@ belongs to a closed
communicating class.

Complexity: excluding shared support-graph construction, the first query
takes @O((n + E) log(n + 1))@ time and @O(n + E)@ temporary space and retains
@O(n)@ component-closedness cache; subsequent queries take @O(1)@ time.
Temporary and result space per cached query are @O(1)@.
-}
recurrentState ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Bool
recurrentState p i = G.inClosedComponent (tmSupport p) (toIndex i)

{- | Test whether @i@ is transient, meaning that its return probability is
less than one. This is the negation of 'recurrentState' for a finite DTMC.

Complexity: excluding shared support-graph construction, the first query
takes @O((n + E) log(n + 1))@ time and @O(n + E)@ temporary space and retains
@O(n)@ component-closedness cache; subsequent queries take @O(1)@ time.
Temporary and result space per cached query are @O(1)@.
-}
transientState ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Bool
transientState p i = not (recurrentState p i)

{- | Return the members of closed communicating classes, ordered by class and
state index. Every non-empty finite DTMC has at least one; the empty chain
returns the empty list.

This is a focused projection of the complete 'classify' report.

Complexity: excluding shared support-graph construction, the first full
evaluation takes @O((n + E) log(n + 1))@ time and @O(n + E)@ temporary space
and retains @O(n)@ component and closedness cache; later evaluations take
@O(n)@ time and temporary space. Result space is @O(n)@.
-}
recurrentStates :: (FiniteState state) => TransitionMatrix state -> [state]
recurrentStates = recurrentStatesOf . classify

{- | Return the members of non-closed communicating classes, ordered by class
and state index. The result is empty exactly when every class is closed.

This is a focused projection of the complete 'classify' report.

Complexity: excluding shared support-graph construction, the first full
evaluation takes @O((n + E) log(n + 1))@ time and @O(n + E)@ temporary space
and retains @O(n)@ component and closedness cache; later evaluations take
@O(n)@ time and temporary space. Result space is @O(n)@.
-}
transientStates :: (FiniteState state) => TransitionMatrix state -> [state]
transientStates = transientStatesOf . classify

{- | Build the complete class, period, recurrence, absorbing-state,
irreducibility, and aperiodicity report from one shared support graph. The
standalone whole-chain queries are focused projections of this report; scalar
reachability, period, and recurrence queries remain direct graph lookups.

Complexity: full evaluation on an unforced matrix takes
@O(n^2 + (n + E) log(n + 1))@ time, @O(n^2 + n + E)@ temporary space,
@O(n + E)@ retained graph-cache space, and @O(n)@ result space. With all
graph facts cached, it takes @O(n)@ time and @O(n)@ temporary and result
space.
-}
classify :: (FiniteState state) => TransitionMatrix state -> Classification state
classify p =
    Classification
        { classesOf = cs
        , isIrreducible = irreducible'
        , isAperiodic = aperiodic'
        , isErgodic = irreducible' && aperiodic'
        , chainPeriod = chainPeriod'
        , recurrentStatesOf = concatMap classMembers (filter classClosed cs)
        , transientStatesOf = concatMap classMembers (filter (not . classClosed) cs)
        , absorbingStates = [i | cc <- cs, classClosed cc, [i] <- [classMembers cc]]
        }
  where
    g = tmSupport p
    cs =
        [ CommClass
            { classMembers = map toState c
            , classPeriod = G.periodOf g v
            , classClosed = G.inClosedComponent g v
            }
        | c@(v : _) <- G.components g
        ]
    irreducible' = graphIrreducible g
    aperiodic' = graphAperiodic g
    chainPeriod' = case cs of
        [c] -> classPeriod c
        _ -> Nothing
