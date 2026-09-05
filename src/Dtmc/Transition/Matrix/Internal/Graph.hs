{- |
Module      : Dtmc.Transition.Matrix.Internal.Graph
Description : Support-graph reachability, components, periods, and phases.

A small DTMC-specific layer over "Data.Graph". It knows nothing about
probabilities: vertices are the integers @{0 .. n-1}@ and edges are the
positive entries of the transition matrix's support.

The graph is stored as adjacency lists in both directions. Keeping the
transpose makes forward and reverse traversals proportional to the graph
actually visited instead of requiring matrix row or column scans.

Unless stated otherwise, every vertex argument must be in @{0 .. V-1}@.
Passing an out-of-range vertex may raise an array-bounds error.

Query complexities assume that the required forward or reverse adjacency
array has already been forced. The first query after 'fromAdjacency' also pays
the documented cost of building that array.
-}
module Dtmc.Transition.Matrix.Internal.Graph (
    Graph,
    graphDim,
    fromAdjacency,
    hasEdge,
    reachable,
    reachesAny,
    backwardReachable,
    components,
    componentOf,
    sameComponent,
    isClosed,
    inClosedComponent,
    componentPeriod,
    periodOf,
    phaseOf,
) where

import Data.Array qualified as Array
import Data.Array.Unboxed qualified as Unboxed
import Data.Graph qualified as DG
import Data.IntMap.Strict qualified as IntMap
import Data.IntSet qualified as IntSet
import Data.List qualified as List
import Data.Sequence qualified as Sequence
import Data.Tree (Tree, flatten)
import Numeric.Natural (Natural)

{- | An immutable directed graph.

'graphSuccessors' contains outgoing neighbours. 'graphPredecessors' is the
transposed graph and therefore contains incoming neighbours. Both describe
the same logical edge set.

Strongly connected components, the component lookup table, the
closed-component table, and the per-vertex period and phase tables are lazy
derived fields. The component structure comes from 'DG.scc'; the period and
phase tables come from one BFS per component. Each field is computed on first
use.
-}
data Graph = Graph
    { graphDim :: Int
    {- ^ Number of vertices @V@.

    Complexity: @O(1)@ time and @O(1)@ space.
    -}
    , graphSuccessors :: DG.Graph
    -- ^ Original graph: the row for @u@ contains every @v@ with @u -> v@.
    , graphPredecessors :: DG.Graph
    -- ^ Transpose: the row for @v@ contains every @u@ with @u -> v@.
    , graphSccs :: [[Int]]
    -- ^ Normalised strongly connected components.
    , graphComponentOf :: Array.Array Int [Int]
    -- ^ Constant-time vertex-to-component lookup after SCC construction.
    , graphComponentId :: Unboxed.UArray Int Int
    {- ^ Constant-time vertex-to-component-index lookup: two vertices share a
    strongly connected component (communicate) iff they map to the same
    index. Backs 'sameComponent'. The table is also used while deriving
    component closedness; retaining it here adds @O(V)@ unboxed storage and
    makes same-component queries @O(1)@.
    -}
    , graphClosedComponentTable :: Unboxed.UArray Int Bool
    {- ^ Per-vertex closedness of its component: @True@ iff the vertex's
    strongly connected component is a sink of the condensation (no edge
    leaves it). Settled in one pass over all edges.
    -}
    , graphPeriodOf :: Array.Array Int (Maybe Natural)
    {- ^ Per-vertex period of its strongly connected component (@Nothing@ when
    the component has no cycles). Filled by one BFS per component.
    -}
    , graphPhaseOf :: Unboxed.UArray Int Int
    {- ^ Per-vertex phase within its component: the BFS level from the
    component's least vertex, modulo the period. Every edge /within a
    component/ advances the phase by one (modulo that period); edges leaving
    a component relate phases across different components and obey no such
    rule. Shares its BFS with 'graphPeriodOf'.
    -}
    }

{- | Build a graph from a complete Boolean adjacency association list. Only
entries whose value is 'True' become edges.

The dimension must be non-negative, and the list is expected to contain
exactly one entry for each pair in @{0 .. V-1}^2@. Completeness and uniqueness
are not validated: missing pairs act as 'False', and repeated 'True' entries
create duplicate edges. A negative dimension raises an error; an out-of-range
endpoint may fail when a lazy field is forced.

Complexity: @O(1)@ initial time and @O(1)@ initial space. For @A@ supplied
entries, forcing the forward adjacency array takes @O(V + A)@ time and
@O(V + E)@ result space; first forcing the transpose takes a further
@O(V + E)@ time and @O(V + E)@ result space. For a complete input, @A = V^2@.
-}
fromAdjacency :: Int -> [((Int, Int), Bool)] -> Graph
fromAdjacency dim entries
    | dim < 0 = error "Dtmc.Transition.Matrix.Internal.Graph.fromAdjacency: negative dimension"
    | otherwise =
        Graph
            { graphDim = dim
            , graphSuccessors = successors
            , graphPredecessors = DG.transposeG successors
            , graphSccs = sccs
            , graphComponentOf = componentTable
            , graphComponentId = componentIds
            , graphClosedComponentTable = closedComponentTable
            , graphPeriodOf = periodTable
            , graphPhaseOf = phaseTable
            }
  where
    successors =
        DG.buildG
            (vertexBounds dim)
            [pair | (pair, present) <- entries, present]

    sccs = normaliseComponents (DG.scc successors)

    componentTable =
        Array.array
            (vertexBounds dim)
            [ (vertex, component)
            | component <- sccs
            , vertex <- component
            ]

    -- A component is closed iff no edge leaves it. Record component ids so one
    -- pass over cross-component edges can mark every open component.
    componentIds :: Unboxed.UArray Int Int
    componentIds =
        Unboxed.array
            (vertexBounds dim)
            [ (vertex, componentIndex)
            | (componentIndex, component) <- zip [0 ..] sccs
            , vertex <- component
            ]

    openComponentIds :: IntSet.IntSet
    openComponentIds =
        IntSet.fromList
            [ componentIds Unboxed.! from
            | from <- [0 .. dim - 1]
            , to <- successors Array.! from
            , componentIds Unboxed.! from /= componentIds Unboxed.! to
            ]

    closedComponentTable :: Unboxed.UArray Int Bool
    closedComponentTable =
        Unboxed.listArray
            (vertexBounds dim)
            [ not (IntSet.member (componentIds Unboxed.! vertex) openComponentIds)
            | vertex <- [0 .. dim - 1]
            ]

    -- One BFS per component yields both its period and each vertex's phase.
    -- The list is a shared thunk, so the period and phase tables never repeat
    -- the traversal.
    componentPhases :: [(Maybe Natural, [(Int, Int)])]
    componentPhases = [componentPhasing successors component | component <- sccs]

    periodTable :: Array.Array Int (Maybe Natural)
    periodTable =
        Array.array
            (vertexBounds dim)
            [ (vertex, period)
            | (period, phases) <- componentPhases
            , (vertex, _) <- phases
            ]

    phaseTable :: Unboxed.UArray Int Int
    phaseTable =
        Unboxed.array
            (vertexBounds dim)
            [ (vertex, phase)
            | (_, phases) <- componentPhases
            , (vertex, phase) <- phases
            ]

vertexBounds :: Int -> (Int, Int)
vertexBounds dim = (0, dim - 1)

vertices :: Graph -> [Int]
vertices graph = [0 .. graphDim graph - 1]

{- | Test whether a direct edge leads from @from@ to @to@.

Algorithms should normally enumerate an adjacency row instead of repeatedly
calling 'hasEdge'.

Complexity: @O(outDegree(from))@ time and @O(1)@ space because a 'Data.Graph'
row is a list.
-}
hasEdge :: Graph -> Int -> Int -> Bool
hasEdge graph from to = to `elem` (graphSuccessors graph Array.! from)

{- | Test whether @to@ is reachable from @from@ in zero or more steps.

This delegates to 'DG.path' and performs a graph search rather than
retaining a quadratic transitive closure.

Complexity: @O(V + E)@ worst-case time and @O(V)@ traversal space per query.
-}
reachable :: Graph -> Int -> Int -> Bool
reachable graph = DG.path (graphSuccessors graph)

{- | Test whether @from@ can reach any supplied target in zero or more steps.

Targets are materialised as a Boolean membership array. The lazy reachable
stream is then consumed until it encounters a target, so the traversal can
terminate early.

Complexity: @O(V + T + E_r)@ worst-case time, where @T@ is the number of
supplied targets and @E_r@ is the portion of the graph examined before
termination; @O(V)@ space.
-}
reachesAny :: Graph -> Int -> [Int] -> Bool
reachesAny _ _ [] = False
reachesAny graph from targets =
    any (targetMask Unboxed.!) reachableVertices
  where
    targetMask :: Unboxed.UArray Int Bool
    targetMask =
        Unboxed.accumArray
            (||)
            False
            (vertexBounds (graphDim graph))
            [(target, True) | target <- targets]

    reachableVertices = DG.reachable (graphSuccessors graph) from

{- | Return the vertices that can reach a seed inside the subgraph induced by
@allowed@.

The allowed predicate is evaluated once per vertex. The transpose is
filtered to the induced subgraph, after which 'DG.dfs' performs one
multi-source traversal. A Boolean result mask restores ascending output
order without an @O(R log R)@ comparison sort.

Complexity: @O(V + E)@ time, excluding the cost of the @V@ predicate calls,
and @O(V + E)@ space for the filtered adjacency lists and traversal state.
-}
backwardReachable :: Graph -> (Int -> Bool) -> [Int] -> [Int]
backwardReachable graph allowed seeds =
    [vertex | vertex <- allVertices, reachedMask Unboxed.! vertex]
  where
    allVertices = vertices graph
    dim = graphDim graph

    allowedMask :: Unboxed.UArray Int Bool
    allowedMask =
        Unboxed.listArray
            (vertexBounds dim)
            (map allowed allVertices)

    isAllowed vertex = allowedMask Unboxed.! vertex

    allowedSeeds = List.filter isAllowed seeds

    restrictedPredecessors :: DG.Graph
    restrictedPredecessors =
        Array.listArray
            (vertexBounds dim)
            [ if isAllowed vertex
                then
                    List.filter
                        isAllowed
                        (graphPredecessors graph Array.! vertex)
                else []
            | vertex <- allVertices
            ]

    reached =
        concatMap flatten (DG.dfs restrictedPredecessors allowedSeeds)

    reachedMask :: Unboxed.UArray Int Bool
    reachedMask =
        Unboxed.accumArray
            (||)
            False
            (vertexBounds dim)
            [(vertex, True) | vertex <- reached]

{- | Return the strongly connected components. Vertices within each component
are in ascending order, and components are ordered by their least vertex.

Complexity: first full evaluation takes @O(V + E + V log V)@ time and @O(V)@
temporary and result space. Later projections take @O(1)@ time and @O(1)@
space before traversal of the cached result.
-}
components :: Graph -> [[Int]]
components = graphSccs

normaliseComponents :: [Tree Int] -> [[Int]]
normaliseComponents =
    List.sortOn componentKey . map (List.sort . flatten)
  where
    componentKey [] = -1
    componentKey (first : _) = first

{- | Return the strongly connected component containing a vertex.

Complexity: the first query takes @O(V + E + V log V)@ time and @O(V)@ cache
space; subsequent queries take @O(1)@ time and @O(1)@ space.
-}
componentOf :: Graph -> Int -> [Int]
componentOf graph vertex
    | vertex < 0 || vertex >= graphDim graph =
        error "Dtmc.Transition.Matrix.Internal.Graph.componentOf: vertex out of bounds"
    | otherwise = graphComponentOf graph Array.! vertex

{- | Test whether two vertices lie in the same strongly connected component;
that is, whether they communicate. This compares their cached component
indices instead of performing two reachability searches.

Complexity: the first query takes @O(V + E + V log V)@ time and @O(V)@ cache
space; subsequent queries take @O(1)@ time and @O(1)@ space.
-}
sameComponent :: Graph -> Int -> Int -> Bool
sameComponent graph a b
    | outOfRange a || outOfRange b =
        error "Dtmc.Transition.Matrix.Internal.Graph.sameComponent: vertex out of bounds"
    | otherwise =
        graphComponentId graph Unboxed.! a == graphComponentId graph Unboxed.! b
  where
    outOfRange v = v < 0 || v >= graphDim graph

{- | Test whether a vertex set is closed: no direct edge leaves it.

Duplicates are ignored, and the empty set is closed.

Complexity: @O(V + S + E_C)@ time, where @S@ is the supplied list length and
@E_C@ is the total out-degree of its vertices; @O(V)@ space for membership.
-}
isClosed :: Graph -> [Int] -> Bool
isClosed graph suppliedVertices =
    all staysInside uniqueVertices
  where
    dim = graphDim graph
    allVertices = vertices graph

    member :: Unboxed.UArray Int Bool
    member =
        Unboxed.accumArray
            (||)
            False
            (vertexBounds dim)
            [(vertex, True) | vertex <- suppliedVertices]

    uniqueVertices =
        [vertex | vertex <- allVertices, member Unboxed.! vertex]

    staysInside from =
        all
            (member Unboxed.!)
            (graphSuccessors graph Array.! from)

{- | Test whether a vertex lies in a closed strongly connected component: a
sink of the condensation with no outgoing edge. In finite-chain terms, this
is exactly recurrence, but that interpretation belongs to the Markov-chain
layer rather than this graph module.

This is the specialised, precomputed form of
@'isClosed' g ('componentOf' g v)@: the open/closed status of every component
is settled once by a pass over all edges and cached, so each later query is a
constant-time array read.

Complexity: the first query takes @O((V + E) log V)@ time and @O(V)@ cache
space; subsequent queries take @O(1)@ time and @O(1)@ space.
-}
inClosedComponent :: Graph -> Int -> Bool
inClosedComponent graph vertex
    | vertex < 0 || vertex >= graphDim graph =
        error "Dtmc.Transition.Matrix.Internal.Graph.inClosedComponent: vertex out of bounds"
    | otherwise = graphClosedComponentTable graph Unboxed.! vertex

{- | Return the period of a strongly connected component: the gcd of the
lengths of all its closed walks. 'Nothing' denotes an empty component or a
singleton with no self-loop.

The input is expected to be a genuine strongly connected component; the value
returned is the period of the component containing its first vertex, read from
the precomputed 'graphPeriodOf' table.

Complexity: an empty input takes @O(1)@ time and @O(1)@ space. The first
non-empty query takes @O((V + E) log V)@ time and @O(V)@ cache space;
subsequent queries take @O(1)@ time and @O(1)@ space.
-}
componentPeriod :: Graph -> [Int] -> Maybe Natural
componentPeriod _ [] = Nothing
componentPeriod graph (root : _) = periodOf graph root

{- | Return the period of the strongly connected component containing the
vertex. Returns 'Nothing' when that component has no cycles. The value is read
from the precomputed table.

Complexity: the first query takes @O((V + E) log V)@ time and @O(V)@ cache
space; subsequent queries take @O(1)@ time and @O(1)@ space.
-}
periodOf :: Graph -> Int -> Maybe Natural
periodOf graph vertex
    | vertex < 0 || vertex >= graphDim graph =
        error "Dtmc.Transition.Matrix.Internal.Graph.periodOf: vertex out of bounds"
    | otherwise = graphPeriodOf graph Array.! vertex

{- | Return the phase of a vertex within its strongly connected component: its
BFS level from the component's least vertex, modulo the component's period
@d@. Every edge @u -> v@ /internal to a component/ satisfies
@phaseOf v == (phaseOf u + 1) `mod` d@. Therefore, grouping a component's
vertices by phase yields its cyclic classes. An edge leaving a component
relates two independent phasings and carries no such relation. A component of
period @d@ has phases in @{0 .. d-1}@; a vertex whose component has no cycles
has phase @0@.

Complexity: the first query takes @O((V + E) log V)@ time and @O(V)@ cache
space; subsequent queries take @O(1)@ time and @O(1)@ space.
-}
phaseOf :: Graph -> Int -> Int
phaseOf graph vertex
    | vertex < 0 || vertex >= graphDim graph =
        error "Dtmc.Transition.Matrix.Internal.Graph.phaseOf: vertex out of bounds"
    | otherwise = graphPhaseOf graph Unboxed.! vertex

-- The sole caller supplies strongly connected components, but retain a
-- defensive fallback for an incomplete BFS: avoid missing-level lookups and
-- assign no period and phase 0. Otherwise one BFS supplies both cached tables.
componentPhasing :: DG.Graph -> [Int] -> (Maybe Natural, [(Int, Int)])
componentPhasing _ [] = (Nothing, [])
componentPhasing successors component@(root : _)
    | not reachedAll = (Nothing, [(vertex, 0) | vertex <- component])
    | period == 0 = (Nothing, [(vertex, 0) | vertex <- component])
    | otherwise =
        ( Just (fromIntegral period)
        , [(vertex, (levels IntMap.! vertex) `mod` period) | vertex <- component]
        )
  where
    member = IntSet.fromList component
    levels = bfsLevels successors member root
    reachedAll = IntMap.size levels == IntSet.size member

    period =
        List.foldl' accumulateVertex 0 (IntSet.toList member)

    accumulateVertex currentGcd from =
        List.foldl'
            (accumulateEdge (levels IntMap.! from))
            currentGcd
            (successors Array.! from)

    accumulateEdge fromLevel currentGcd to
        | not (IntSet.member to member) = currentGcd
        | otherwise =
            gcd currentGcd (abs (fromLevel + 1 - levels IntMap.! to))

-- Breadth-first levels within one component. Vertices are inserted into the
-- level map when enqueued, so each is enqueued exactly once.
bfsLevels :: DG.Graph -> IntSet.IntSet -> Int -> IntMap.IntMap Int
bfsLevels successors member root =
    search (Sequence.singleton root) (IntMap.singleton root 0)
  where
    search queue levels =
        case Sequence.viewl queue of
            Sequence.EmptyL -> levels
            from Sequence.:< rest ->
                search queue' levels'
              where
                fromLevel = levels IntMap.! from
                (queue', levels') =
                    List.foldl'
                        (discover fromLevel)
                        (rest, levels)
                        (successors Array.! from)

    discover fromLevel state@(queue, levels) candidate
        | not (IntSet.member candidate member) = state
        | IntMap.member candidate levels = state
        | otherwise =
            ( queue Sequence.|> candidate
            , IntMap.insert candidate (fromLevel + 1) levels
            )
