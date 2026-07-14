-- |
-- Module      : Dtmc.Classification
-- Description : Communicating classes, irreducibility, and periodicity.
--
-- Qualitative structure of a chain derived purely from the support graph of @P@
-- (a directed edge @i -> j@ whenever @P(i,j) > 0@). Because it depends only on
-- which entries are positive, everything here is exact combinatorics, independent
-- of the actual probabilities: reachability, the partition into communicating
-- classes, irreducibility, and the period of each state/class.
module Dtmc.Classification (
    -- * Support graph
    supportEdge,
    accessible,
    communicates,

    -- * Communicating classes
    communicatingClasses,
    irreducible,

    -- * Periodicity
    period,
    aperiodic,

    -- * Classification summary
    CommClass (..),
    Classification,
    classesOf,
    classify,

    -- * Irreducibility witness
    Irreducible,
    witnessIrreducible,
    irreducibleMatrix,
) where

import Data.Array (
    Array,
    array,
    (!),
 )
import Data.Finite (
    Finite,
    finite,
    getFinite,
 )
import Data.Maybe (
    fromMaybe,
    isNothing,
 )
import Dtmc.Internal.Types (TransitionMatrix, unTransitionMatrix)
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (Natural)

-- The support graph of a chain, built once and shared by every query below.
-- @graphAdjacency@ is the direct support relation (@P(i,j) > 0@);
-- @graphReach@ is its reflexive-transitive closure. Fields are lazy, so a query
-- that needs only adjacency (e.g. 'supportEdge') never forces the closure.
data Graph = Graph
    { graphDim :: Int
    , graphAdjacency :: Array (Int, Int) Bool
    , graphReach :: Array (Int, Int) Bool
    }

-- Build the support graph and its closure from a transition matrix. The
-- adjacency array records @P(i,j) > 0@; 'reachClosure' then adds all indirect
-- paths (and the reflexive diagonal).
buildGraph :: (KnownNat n) => TransitionMatrix n -> Graph
buildGraph p =
    Graph
        { graphDim = dim
        , graphAdjacency = adjacency
        , graphReach = reachClosure dim adjacency
        }
  where
    rows = LA.toLists (S.extract (unTransitionMatrix p))
    dim = length rows
    adjacency =
        array
            ((0, 0), (dim - 1, dim - 1))
            [ ((i, j), entry > 0)
            | (i, row) <- zip [0 ..] rows
            , (j, entry) <- zip [0 ..] row
            ]

-- Reflexive-transitive closure of a boolean adjacency array by Floyd-Warshall:
-- start from the adjacency relation with the diagonal forced true, then for each
-- intermediate vertex @k@ add the edge @i -> j@ whenever @i -> k@ and @k -> j@.
-- Entry @(i,j)@ ends true iff @j@ is reachable from @i@ in zero or more steps.
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

-- Convert a raw @Int@ state index into the bounded 'Finite' @n@ index.
toFinite :: (KnownNat n) => Int -> Finite n
toFinite = finite . fromIntegral

-- Convert a bounded 'Finite' @n@ index back to a raw @Int@.
toIndex :: Finite n -> Int
toIndex = fromIntegral . getFinite

-- Partition state indices into communicating classes given a built graph:
-- greedily group each state with all later states it communicates with. States
-- (and classes) come out in ascending-index order.
rawClasses :: Graph -> [[Int]]
rawClasses g = go [0 .. graphDim g - 1]
  where
    reach = graphReach g
    comm i j = reach ! (i, j) && reach ! (j, i)
    go [] = []
    go (x : xs) =
        (x : filter (comm x) xs) : go (filter (not . comm x) xs)

-- | Direct one-step reachability: @True@ iff @P(i,j) > 0@.
supportEdge :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
supportEdge p i j = graphAdjacency (buildGraph p) ! (toIndex i, toIndex j)

-- | Accessibility @i -> j@: @j@ is reachable from @i@ in zero or more steps.
accessible :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
accessible p i j = graphReach (buildGraph p) ! (toIndex i, toIndex j)

-- | Communication @i <-> j@: @i@ and @j@ are mutually accessible. This is the
-- equivalence relation whose classes are the communicating classes.
communicates :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
communicates p i j =
    reach ! (a, b) && reach ! (b, a)
  where
    reach = graphReach (buildGraph p)
    a = toIndex i
    b = toIndex j

-- | The communicating classes of the chain, each as a list of states.
communicatingClasses :: (KnownNat n) => TransitionMatrix n -> [[Finite n]]
communicatingClasses p =
    map (map toFinite) (rawClasses (buildGraph p))

-- | Whether the chain is irreducible: all states form a single (non-empty)
-- communicating class, so every state is reachable from every other.
irreducible :: (KnownNat n) => TransitionMatrix n -> Bool
irreducible p =
    case rawClasses (buildGraph p) of
        [c] -> not (null c)
        _ -> False

-- | Period of state @i@: the gcd of the lengths of all closed walks through
-- @i@, computed within @i@'s communicating class. 'Nothing' when the class has
-- no cycles (the period is undefined).
period :: (KnownNat n) => TransitionMatrix n -> Finite n -> Maybe Natural
period p i =
    periodOfClass (edgeOf g) klass
  where
    g = buildGraph p
    reach = graphReach g
    a = toIndex i
    klass = filter (\j -> reach ! (a, j) && reach ! (j, a)) [0 .. graphDim g - 1]

-- Read the support (adjacency) relation of a graph as a plain edge predicate,
-- for the period BFS below.
edgeOf :: Graph -> Int -> Int -> Bool
edgeOf g u v = graphAdjacency g ! (u, v)

-- Period of a communicating class: BFS-label the class from a root, then take
-- the gcd of @level(u) + 1 - level(v)@ over every intra-class edge @u -> v@.
-- A gcd of @0@ (no edges/cycles) means the period is undefined.
periodOfClass :: (Int -> Int -> Bool) -> [Int] -> Maybe Natural
periodOfClass _ [] = Nothing
periodOfClass edge klass@(root : _) =
    if d == 0 then Nothing else Just (fromIntegral d)
  where
    dist = bfsWithin edge klass root
    lvl u = fromMaybe 0 (lookup u dist)
    d =
        foldl'
            gcd
            0
            [ abs (lvl u + 1 - lvl v)
            | u <- klass
            , v <- klass
            , edge u v
            ]

-- Breadth-first level assignment restricted to a single class, returning each
-- reached vertex paired with its distance from @root@.
bfsWithin :: (Int -> Int -> Bool) -> [Int] -> Int -> [(Int, Int)]
bfsWithin edge klass root = go [root] [(root, 0)]
  where
    go [] dist = dist
    go (u : queue) dist =
        go (queue ++ fresh) (dist ++ [(v, du + 1) | v <- fresh])
      where
        du = fromMaybe 0 (lookup u dist)
        fresh =
            [ v
            | v <- klass
            , edge u v
            , isNothing (lookup v dist)
            , v `notElem` queue
            ]

-- | Whether the chain is aperiodic: every communicating class has period @1@
-- (and there is at least one state). Derived from 'classify' so the support
-- graph and closure are built only once, rather than per state.
aperiodic :: (KnownNat n) => TransitionMatrix n -> Bool
aperiodic p =
    not (null classes) && all ((== Just 1) . classPeriod) classes
  where
    classes = classesOf (classify p)

-- | A summary of one communicating class: its member states, its 'period'
-- (@Nothing@ if undefined), and whether it is 'classClosed' -- i.e. no edge
-- leaves the class, making it a recurrent/absorbing set.
data CommClass n = CommClass
    { classMembers :: [Finite n]
    , classPeriod :: Maybe Natural
    , classClosed :: Bool
    }

deriving instance (KnownNat n) => Eq (CommClass n)

deriving instance (KnownNat n) => Show (CommClass n)

-- | A full decomposition of a chain into its communicating classes.
newtype Classification n = Classification [CommClass n]

type role Classification nominal

deriving instance (KnownNat n) => Eq (Classification n)

deriving instance (KnownNat n) => Show (Classification n)

-- | Extract the list of communicating classes from a 'Classification'.
classesOf :: Classification n -> [CommClass n]
classesOf (Classification cs) = cs

-- | Decompose a chain into communicating classes, annotating each with its
-- period and whether it is closed.
classify :: (KnownNat n) => TransitionMatrix n -> Classification n
classify p =
    Classification
        [ CommClass
            { classMembers = map toFinite c
            , classPeriod = periodOfClass edge c
            , classClosed = isClosed c
            }
        | c <- classes
        ]
  where
    g = buildGraph p
    dim = graphDim g
    edge = edgeOf g
    classes = rawClasses g
    isClosed c =
        and
            [ not (edge u v)
            | u <- c
            , v <- [0 .. dim - 1]
            , v `notElem` c
            ]

-- | A transition matrix carrying a proof that it is irreducible, obtainable
-- only via 'witnessIrreducible'. Lets downstream code demand irreducibility in
-- a type rather than re-checking it.
newtype Irreducible n = Irreducible (TransitionMatrix n)

type role Irreducible nominal

deriving instance (KnownNat n) => Show (Irreducible n)

-- | Certify irreducibility: @Just@ the wrapped matrix when 'irreducible' holds,
-- @Nothing@ otherwise.
witnessIrreducible :: (KnownNat n) => TransitionMatrix n -> Maybe (Irreducible n)
witnessIrreducible p
    | irreducible p = Just (Irreducible p)
    | otherwise = Nothing

-- | Recover the underlying transition matrix from an 'Irreducible' witness.
irreducibleMatrix :: Irreducible n -> TransitionMatrix n
irreducibleMatrix (Irreducible p) = p
