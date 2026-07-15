-- |
-- Module      : Dtmc.Classification
-- Description : Communicating classes, irreducibility, and periodicity.
--
-- Qualitative structure of a chain derived purely from the support graph of @P@
-- (a directed edge @i -> j@ whenever @P(i,j) > 0@). Because it depends only on
-- which entries are positive, everything here is exact combinatorics, independent
-- of the actual probabilities: reachability, the partition into communicating
-- classes, irreducibility, and the period of each state/class.
--
-- The graph combinatorics live in "Dtmc.Internal.Graph"; this module is the thin
-- bridge that reads a chain's support into a 'Graph' (once, via 'supportGraph')
-- and renames the graph facts into Markov-chain vocabulary: communicating
-- classes are the strongly connected components, and a class's period is the
-- period of that component.
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
    getFinite,
 )
import Dtmc.Internal.Graph (
    Graph,
    closed,
    componentOf,
    componentPeriod,
    components,
    edge,
    fromAdjacency,
    reachable,
 )
import Dtmc.Internal.Types (TransitionMatrix, unTransitionMatrix)
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (Natural)

-- The one bridge from probabilities to combinatorics: build the support graph
-- of @P@ (edge @i -> j@ iff @P(i,j) > 0@). Every function below goes through
-- this and then speaks only in graph terms.
supportGraph :: (KnownNat n) => TransitionMatrix n -> Graph
supportGraph p =
    fromAdjacency
        dim
        [ ((i, j), entry > 0)
        | (i, row) <- zip [0 ..] rows
        , (j, entry) <- zip [0 ..] row
        ]
  where
    rows = LA.toLists (S.extract (unTransitionMatrix p))
    dim = length rows

-- Convert a raw @Int@ state index into the bounded 'Finite' @n@ index.
toFinite :: (KnownNat n) => Int -> Finite n
toFinite = finite . fromIntegral

-- Convert a bounded 'Finite' @n@ index back to a raw @Int@.
toIndex :: Finite n -> Int
toIndex = fromIntegral . getFinite

-- | Direct one-step reachability: @True@ iff @P(i,j) > 0@.
supportEdge :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
supportEdge p i j = edge (supportGraph p) (toIndex i) (toIndex j)

-- | Accessibility @i -> j@: @j@ is reachable from @i@ in zero or more steps.
accessible :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
accessible p i j = reachable (supportGraph p) (toIndex i) (toIndex j)

-- | Communication @i <-> j@: @i@ and @j@ are mutually accessible. This is the
-- equivalence relation whose classes are the communicating classes.
communicates :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
communicates p i j =
    reachable g a b && reachable g b a
  where
    g = supportGraph p
    a = toIndex i
    b = toIndex j

-- | The communicating classes of the chain, each as a list of states.
communicatingClasses :: (KnownNat n) => TransitionMatrix n -> [[Finite n]]
communicatingClasses p =
    map (map toFinite) (components (supportGraph p))

-- | Whether the chain is irreducible: all states form a single (non-empty)
-- communicating class, so every state is reachable from every other.
irreducible :: (KnownNat n) => TransitionMatrix n -> Bool
irreducible p =
    case components (supportGraph p) of
        [c] -> not (null c)
        _ -> False

-- | Period of state @i@: the gcd of the lengths of all closed walks through
-- @i@, computed within @i@'s communicating class. 'Nothing' when the class has
-- no cycles (the period is undefined).
period :: (KnownNat n) => TransitionMatrix n -> Finite n -> Maybe Natural
period p i =
    componentPeriod g (componentOf g (toIndex i))
  where
    g = supportGraph p

-- | Whether the chain is aperiodic: every communicating class has period @1@
-- (and there is at least one state). A class with undefined period -- a
-- cycle-free transient class, where 'period' is 'Nothing' -- does not have
-- period one, so its presence makes the chain non-aperiodic under this
-- definition. Derived from 'classify' so the support graph and closure are
-- built only once, rather than per state.
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
            , classPeriod = componentPeriod g c
            , classClosed = closed g c
            }
        | c <- components g
        ]
  where
    g = supportGraph p

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
