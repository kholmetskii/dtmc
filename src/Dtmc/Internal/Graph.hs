-- |
-- Module      : Dtmc.Internal.Graph
-- Description : Support-graph combinatorics: reachability, components, period.
--
-- A small, purpose-built directed-graph layer over a boolean adjacency relation.
-- It knows nothing about probabilities or transition matrices: it is exact
-- combinatorics on the vertex set @{0 .. n-1}@. "Dtmc.Classification" builds a
-- 'Graph' from the support of a chain (the edges @i -> j@ with @P(i,j) > 0@) and
-- renames these graph facts into Markov-chain vocabulary -- communicating
-- classes are the strongly connected components, and a class's period is a
-- property of that component.
--
-- This is deliberately not a general graph library: it exposes only the few
-- operations the qualitative theory needs. The adjacency and its
-- reflexive-transitive closure are computed once in 'fromAdjacency' (the closure
-- by Floyd-Warshall over an immutable 'Array', @O(n^3)@ with @O(1)@ lookups),
-- and every query below reads from that.
module Dtmc.Internal.Graph (
    Graph,
    graphDim,
    fromAdjacency,
    edge,
    reachable,
    components,
    componentOf,
    closed,
    componentPeriod,
) where

import Data.Array (
    Array,
    array,
    (!),
 )
import Data.Maybe (
    fromMaybe,
    isNothing,
 )
import Numeric.Natural (Natural)

-- | A directed graph on @{0 .. graphDim-1}@, holding both the direct adjacency
-- relation and its reflexive-transitive closure. Fields are lazy, so a query
-- that needs only adjacency (e.g. 'edge') never forces the closure.
data Graph = Graph
    { graphDim :: Int
    -- ^ Number of vertices.
    , graphAdjacency :: Array (Int, Int) Bool
    , graphReach :: Array (Int, Int) Bool
    }

-- | Build a graph from its dimension and the adjacency relation as an
-- association list @((i,j), i -> j)@. The reflexive-transitive closure is
-- computed here, once.
fromAdjacency :: Int -> [((Int, Int), Bool)] -> Graph
fromAdjacency dim entries =
    Graph
        { graphDim = dim
        , graphAdjacency = adjacency
        , graphReach = reachClosure dim adjacency
        }
  where
    adjacency = array ((0, 0), (dim - 1, dim - 1)) entries

-- Reflexive-transitive closure of a boolean adjacency array by Floyd-Warshall:
-- start from adjacency with the diagonal forced true, then for each intermediate
-- vertex @k@ add @i -> j@ whenever @i -> k@ and @k -> j@. Entry @(i,j)@ ends true
-- iff @j@ is reachable from @i@ in zero or more steps.
reachClosure :: Int -> Array (Int, Int) Bool -> Array (Int, Int) Bool
reachClosure dim adjacency =
    foldl' pass reflexive [0 .. dim - 1]
  where
    bounds' = ((0, 0), (dim - 1, dim - 1))
    indices = [0 .. dim - 1]
    reflexive =
        array
            bounds'
            [ ((i, j), i == j || adjacency ! (i, j))
            | i <- indices
            , j <- indices
            ]
    pass reach k =
        array
            bounds'
            [ ((i, j), reach ! (i, j) || (reach ! (i, k) && reach ! (k, j)))
            | i <- indices
            , j <- indices
            ]

-- | Direct one-step edge: @True@ iff the adjacency relation holds for @(i,j)@.
edge :: Graph -> Int -> Int -> Bool
edge g u v = graphAdjacency g ! (u, v)

-- | Reachability @i -> j@ in zero or more steps.
reachable :: Graph -> Int -> Int -> Bool
reachable g u v = graphReach g ! (u, v)

-- | Strongly connected components (mutual-reachability classes), each in
-- ascending-vertex order, and themselves ordered by least member.
components :: Graph -> [[Int]]
components g = go [0 .. graphDim g - 1]
  where
    mutual i j = reachable g i j && reachable g j i
    go [] = []
    go (x : xs) =
        (x : filter (mutual x) xs) : go (filter (not . mutual x) xs)

-- | The component containing a given vertex, in ascending-vertex order (so its
-- head is the least member -- matching the corresponding element of
-- 'components').
componentOf :: Graph -> Int -> [Int]
componentOf g v =
    filter (\u -> reachable g v u && reachable g u v) [0 .. graphDim g - 1]

-- | Whether a vertex set is closed: no edge leaves it. A closed component is a
-- recurrent/absorbing set of the chain.
closed :: Graph -> [Int] -> Bool
closed g vertices =
    and
        [ not (edge g u v)
        | u <- vertices
        , v <- [0 .. graphDim g - 1]
        , v `notElem` vertices
        ]

-- | Period of a component: BFS-label it from its head, then take the gcd of
-- @level(u) + 1 - level(v)@ over every internal edge @u -> v@. A gcd of @0@ (no
-- edges/cycles) means the period is undefined.
componentPeriod :: Graph -> [Int] -> Maybe Natural
componentPeriod _ [] = Nothing
componentPeriod g component@(root : _) =
    if d == 0 then Nothing else Just (fromIntegral d)
  where
    dist = bfsWithin g component root
    lvl u = fromMaybe 0 (lookup u dist)
    d =
        foldl'
            gcd
            0
            [ abs (lvl u + 1 - lvl v)
            | u <- component
            , v <- component
            , edge g u v
            ]

-- Breadth-first level assignment restricted to a single component, returning
-- each reached vertex paired with its distance from @root@.
bfsWithin :: Graph -> [Int] -> Int -> [(Int, Int)]
bfsWithin g component root = go [root] [(root, 0)]
  where
    go [] dist = dist
    go (u : queue) dist =
        go (queue ++ fresh) (dist ++ [(v, du + 1) | v <- fresh])
      where
        du = fromMaybe 0 (lookup u dist)
        fresh =
            [ v
            | v <- component
            , edge g u v
            , isNothing (lookup v dist)
            , v `notElem` queue
            ]
