-- |
-- Module      : Dtmc.Internal.Graph
-- Description : Support-graph combinatorics: reachability, components, recurrence, period, phase.
--
-- A small DTMC-specific layer over "Data.Graph". It knows nothing about
-- probabilities: vertices are the integers @{0 .. n-1}@ and edges are the
-- positive entries of the transition matrix's support.
--
-- The graph is stored as adjacency lists in both directions. Keeping the
-- transpose makes forward and reverse traversals proportional to the graph
-- actually visited instead of requiring matrix row or column scans.
module Dtmc.Internal.Graph (
    Graph,
    graphDim,
    fromAdjacency,
    hasEdge,
    reachable,
    reachesAny,
    backwardReachable,
    components,
    componentOf,
    isClosed,
    inClosedComponent,
    componentPeriod,
    periodOf,
    phaseOf,
) where

import qualified Data.Array as Array
import qualified Data.Array.Unboxed as Unboxed
import qualified Data.Graph as DG
import qualified Data.IntMap.Strict as IntMap
import qualified Data.IntSet as IntSet
import qualified Data.List as List
import qualified Data.Sequence as Sequence
import Data.Tree (Tree, flatten)
import Numeric.Natural (Natural)

-- | An immutable directed graph.
--
-- 'graphSuccessors' contains outgoing neighbours. 'graphPredecessors' is the
-- transposed graph and therefore contains incoming neighbours. Both describe
-- the same logical edge set.
--
-- Strongly connected components, the component lookup table, the
-- closed-component table, and the per-vertex period/phase tables are lazy
-- derived fields. The
-- component structure comes from 'DG.scc'; the period/phase tables come from one
-- BFS per component. All are computed on first use.
data Graph = Graph
    { graphDim :: Int
    -- ^ Number of vertices @V@.
    , graphSuccessors :: DG.Graph
    -- ^ Original graph: the row for @u@ contains every @v@ with @u -> v@.
    , graphPredecessors :: DG.Graph
    -- ^ Transpose: the row for @v@ contains every @u@ with @u -> v@.
    , graphSccs :: [[Int]]
    -- ^ Normalised strongly connected components.
    , graphComponentOf :: Array.Array Int [Int]
    -- ^ Constant-time vertex-to-component lookup after SCC construction.
    , graphClosedComponentTable :: Unboxed.UArray Int Bool
    -- ^ Per-vertex closedness of its component: @True@ iff the vertex's
    -- strongly connected component is a sink of the condensation (no edge
    -- leaves it). Settled in one pass over all edges.
    , graphPeriodOf :: Array.Array Int (Maybe Natural)
    -- ^ Per-vertex period of its strongly connected component (@Nothing@ when
    -- the component has no cycles). Filled by one BFS per component.
    , graphPhaseOf :: Unboxed.UArray Int Int
    -- ^ Per-vertex phase within its component: the BFS level from the
    -- component's least vertex, modulo the period. Every edge /within a
    -- component/ advances the phase by one (modulo that period); edges leaving
    -- a component relate phases across different components and obey no such
    -- rule. Shares its BFS with 'graphPeriodOf'.
    }

-- | Build a graph from its dimension and a complete Boolean adjacency
-- association list. Only entries whose value is 'True' become edges.
--
-- The input format contains @V^2@ entries, so reading it necessarily takes
-- @O(V^2)@ even when the resulting support graph is sparse. Once built, the
-- graph itself occupies @O(V + E)@ space.
fromAdjacency :: Int -> [((Int, Int), Bool)] -> Graph
fromAdjacency dim entries
    | dim < 0 = error "Dtmc.Internal.Graph.fromAdjacency: negative dimension"
    | otherwise =
        Graph
            { graphDim = dim
            , graphSuccessors = successors
            , graphPredecessors = DG.transposeG successors
            , graphSccs = sccs
            , graphComponentOf = componentTable
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

    -- Component id per vertex, then the closed-component table. A component is
    -- closed (its states recurrent) iff no edge leaves it -- i.e. it is a sink
    -- of the condensation. One pass over all edges settles every component: an
    -- edge whose endpoints lie in different components marks the source
    -- component open.
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

-- | Direct one-step edge test.
--
-- Time: @O(outDegree(u))@ because a 'Data.Graph' row is a list. Algorithms
-- should normally enumerate the row instead of repeatedly calling 'hasEdge'.
hasEdge :: Graph -> Int -> Int -> Bool
hasEdge graph from to = to `elem` (graphSuccessors graph Array.! from)

-- | Whether @to@ is reachable from @from@ in zero or more steps.
--
-- This delegates to 'DG.path' and performs a graph search rather than
-- retaining a quadratic transitive closure.
--
-- Time: @O(V + E)@ worst case per query. Space: @O(V)@ traversal state.
reachable :: Graph -> Int -> Int -> Bool
reachable graph = DG.path (graphSuccessors graph)

-- | Whether @from@ can reach any of the supplied targets in zero or more
-- steps.
--
-- Targets are materialised as a Boolean membership array. The lazy reachable
-- stream is then consumed until it encounters a target, so the traversal can
-- terminate early.
--
-- Time: @O(V + T + E_r)@ worst case, where @T@ is the number of supplied
-- targets and @E_r@ is the portion of the graph examined before termination.
-- Space: @O(V)@.
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

-- | Reverse reachability from a seed set inside the subgraph induced by
-- @allowed@.
--
-- The allowed predicate is evaluated once per vertex. The transpose is
-- filtered to the induced subgraph, after which 'DG.dfs' performs one
-- multi-source traversal. A Boolean result mask restores ascending output
-- order without an @O(R log R)@ comparison sort.
--
-- Time: @O(V + E)@, excluding the cost of the @V@ predicate calls.
-- Space: @O(V + E)@ for the filtered adjacency lists and traversal state.
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

-- | Strongly connected components. Vertices within each component are in
-- ascending order, and components are ordered by their least vertex.
--
-- The 'DG.scc' graph phase is @O(V + E)@. Normalising the stable public output
-- order adds @O(V log V)@ worst-case sorting work.
components :: Graph -> [[Int]]
components = graphSccs

normaliseComponents :: [Tree Int] -> [[Int]]
normaliseComponents =
    List.sortOn componentKey . map (List.sort . flatten)
  where
    componentKey [] = -1
    componentKey (first : _) = first

-- | The strongly connected component containing a vertex.
--
-- Time: @O(1)@ after the one-off SCC and lookup-table construction.
componentOf :: Graph -> Int -> [Int]
componentOf graph vertex
    | vertex < 0 || vertex >= graphDim graph =
        error "Dtmc.Internal.Graph.componentOf: vertex out of bounds"
    | otherwise = graphComponentOf graph Array.! vertex

-- | Whether a vertex set is closed: no direct edge leaves it.
--
-- This examines only actual outgoing edges of vertices in the set.
--
-- Time: @O(V + E_C)@, where @E_C@ is the total out-degree of the set's vertices
-- (all their outgoing edges, not only the ones that leave the set).
-- Space: @O(V)@ for constant-time set membership.
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

-- | Whether the vertex lies in a closed strongly connected component -- one
-- that is a sink of the condensation, with no edge leaving it. (In finite-chain
-- terms this is exactly recurrence, but that reading belongs to the
-- Markov-vocabulary layer, not to this graph module.)
--
-- This is the specialised, precomputed form of
-- @'isClosed' g ('componentOf' g v)@: the open/closed status of every component
-- is settled once in a single @O(V + E)@ pass over all edges and cached, so each
-- query is a constant-time array read.
--
-- Time: @O(1)@ after the one-off pass. Space: @O(1)@ per query.
inClosedComponent :: Graph -> Int -> Bool
inClosedComponent graph vertex
    | vertex < 0 || vertex >= graphDim graph =
        error "Dtmc.Internal.Graph.inClosedComponent: vertex out of bounds"
    | otherwise = graphClosedComponentTable graph Unboxed.! vertex

-- | Period of a strongly connected component: the gcd of the lengths of all
-- its closed walks. 'Nothing' denotes an empty component or a singleton with
-- no self-loop.
--
-- The input is expected to be a genuine strongly connected component; the value
-- returned is the period of the component containing its first vertex, read from
-- the precomputed 'graphPeriodOf' table.
--
-- Time: @O(1)@ after the one-off @O(V + E)@ phasing pass. Space: @O(1)@.
componentPeriod :: Graph -> [Int] -> Maybe Natural
componentPeriod _ [] = Nothing
componentPeriod graph (root : _) = periodOf graph root

-- | Period of the strongly connected component containing the vertex
-- (@Nothing@ when that component has no cycles), read from the precomputed
-- table.
--
-- Time: @O(1)@ after the one-off phasing pass. Space: @O(1)@ per query.
periodOf :: Graph -> Int -> Maybe Natural
periodOf graph vertex
    | vertex < 0 || vertex >= graphDim graph =
        error "Dtmc.Internal.Graph.periodOf: vertex out of bounds"
    | otherwise = graphPeriodOf graph Array.! vertex

-- | Phase of the vertex within its strongly connected component: its BFS level
-- from the component's least vertex, modulo the component's period @d@. Every
-- edge @u -> v@ /internal to a component/ satisfies
-- @phaseOf v == (phaseOf u + 1) `mod` d@ -- so grouping a component's vertices
-- by phase yields its cyclic (periodicity) classes -- whereas an edge leaving a
-- component relates two independent phasings and carries no such relation. A
-- vertex whose component has no cycles has phase @0@.
--
-- Time: @O(1)@ after the one-off phasing pass. Space: @O(1)@ per query.
phaseOf :: Graph -> Int -> Int
phaseOf graph vertex
    | vertex < 0 || vertex >= graphDim graph =
        error "Dtmc.Internal.Graph.phaseOf: vertex out of bounds"
    | otherwise = graphPhaseOf graph Unboxed.! vertex

-- Period and phases of one component from a single BFS. @Nothing@ period for a
-- component with no cycles (a singleton without a self-loop, or -- defensively
-- -- a set the BFS fails to cover), in which case every supplied vertex is
-- given phase 0. Otherwise the period is the gcd of the edge-level
-- discrepancies and each vertex's phase is its BFS level modulo that period.
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
