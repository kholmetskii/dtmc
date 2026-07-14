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

import Data.Finite (
    Finite,
    finite,
    finites,
    getFinite,
 )
import Data.Maybe (
    fromMaybe, isNothing,
 )
import Dtmc.Internal.Types ( TransitionMatrix, unTransitionMatrix )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (Natural)

-- Adjacency of the support graph: entry @(i,j)@ is true iff @P(i,j) > 0@.
supportMatrix :: (KnownNat n) => TransitionMatrix n -> [[Bool]]
supportMatrix p =
    map (map (> 0)) (LA.toLists (S.extract (unTransitionMatrix p)))

-- Index a boolean adjacency/reachability matrix as @m[i][j]@.
at :: [[Bool]] -> Int -> Int -> Bool
at m i j = (m !! i) !! j

-- Convert a raw @Int@ state index into the bounded 'Finite' @n@ index.
toFinite :: (KnownNat n) => Int -> Finite n
toFinite = finite . fromIntegral

-- Convert a bounded 'Finite' @n@ index back to a raw @Int@.
toIndex :: Finite n -> Int
toIndex = fromIntegral . getFinite

-- Reflexive-transitive closure of the support graph by @n@ rounds of
-- Boolean relation composition: entry @(i,j)@ becomes true iff @j@ is reachable
-- from @i@ in zero or more steps.
reachClosure :: [[Bool]] -> [[Bool]]
reachClosure s = iterate stepClosure start !! dim
  where
    dim = length s
    idxs = [0 .. dim - 1]
    start = [[i == j || at s i j | j <- idxs] | i <- idxs]
    stepClosure r =
        [ [at r i j || or [at r i k && at s k j | k <- idxs] | j <- idxs]
        | i <- idxs
        ]

-- | Direct one-step reachability: @True@ iff @P(i,j) > 0@.
supportEdge :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
supportEdge p i j = at (supportMatrix p) (toIndex i) (toIndex j)

-- | Accessibility @i -> j@: @j@ is reachable from @i@ in zero or more steps.
accessible :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
accessible p i j = at (reachClosure (supportMatrix p)) (toIndex i) (toIndex j)

-- | Communication @i <-> j@: @i@ and @j@ are mutually accessible. This is the
-- equivalence relation whose classes are the communicating classes.
communicates :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
communicates p i j =
    at r a b && at r b a
  where
    r = reachClosure (supportMatrix p)
    a = toIndex i
    b = toIndex j

-- Partition state indices into communicating classes given a reachability
-- matrix: greedily group each state with all others it communicates with.
rawClasses :: [[Bool]] -> [[Int]]
rawClasses r = go [0 .. length r - 1]
  where
    comm i j = at r i j && at r j i
    go [] = []
    go (x : xs) =
        (x : filter (comm x) xs) : go (filter (not . comm x) xs)

-- | The communicating classes of the chain, each as a list of states.
communicatingClasses :: (KnownNat n) => TransitionMatrix n -> [[Finite n]]
communicatingClasses p =
    map (map toFinite) (rawClasses (reachClosure (supportMatrix p)))

-- | Whether the chain is irreducible: all states form a single (non-empty)
-- communicating class, so every state is reachable from every other.
irreducible :: (KnownNat n) => TransitionMatrix n -> Bool
irreducible p =
    case rawClasses (reachClosure (supportMatrix p)) of
        [c] -> not (null c)
        _ -> False


-- | Period of state @i@: the gcd of the lengths of all closed walks through
-- @i@, computed within @i@'s communicating class. 'Nothing' when the class has
-- no cycles (the period is undefined).
period :: (KnownNat n) => TransitionMatrix n -> Finite n -> Maybe Natural
period p i =
    periodOfClass s klass
  where
    s = supportMatrix p
    r = reachClosure s
    a = toIndex i
    klass = filter (\j -> at r a j && at r j a) [0 .. length s - 1]

-- Period of a communicating class: BFS-label the class from a root, then take
-- the gcd of @level(u) + 1 - level(v)@ over every intra-class edge @u -> v@.
-- A gcd of @0@ (no edges/cycles) means the period is undefined.
periodOfClass :: [[Bool]] -> [Int] -> Maybe Natural
periodOfClass _ [] = Nothing
periodOfClass s klass@(root : _) =
    if d == 0 then Nothing else Just (fromIntegral d)
  where
    dist = bfsWithin s klass root
    lvl u = fromMaybe 0 (lookup u dist)
    d =
        foldl'
            gcd
            0
            [ abs (lvl u + 1 - lvl v)
            | u <- klass
            , v <- klass
            , at s u v
            ]

-- Breadth-first level assignment restricted to a single class, returning each
-- reached vertex paired with its distance from @root@.
bfsWithin :: [[Bool]] -> [Int] -> Int -> [(Int, Int)]
bfsWithin s klass root = go [root] [(root, 0)]
  where
    go [] dist = dist
    go (u : queue) dist =
        go (queue ++ fresh) (dist ++ [(v, du + 1) | v <- fresh])
      where
        du = fromMaybe 0 (lookup u dist)
        fresh =
            [ v
            | v <- klass
            , at s u v
            , isNothing (lookup v dist)
            , v `notElem` queue
            ]

-- | Whether the chain is aperiodic: every state has period @1@ (and there is
-- at least one state).
aperiodic :: forall n. (KnownNat n) => TransitionMatrix n -> Bool
aperiodic p =
    not (null states) && all (\i -> period p i == Just 1) states
  where
    states = finites :: [Finite n]

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
            , classPeriod = periodOfClass s c
            , classClosed = isClosed c
            }
        | c <- classes
        ]
  where
    s = supportMatrix p
    dim = length s
    classes = rawClasses (reachClosure s)
    isClosed c =
        and
            [ not (at s u v)
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
