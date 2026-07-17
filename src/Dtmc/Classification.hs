-- |
-- Module      : Dtmc.Classification
-- Description : Communicating classes, irreducibility, periodicity, and recurrence.
--
-- Qualitative structure of a chain derived purely from the support graph of @P@
-- (a directed edge @i -> j@ whenever @P(i,j) > 0@). Because it depends only on
-- which entries are positive, everything here is exact combinatorics, independent
-- of the actual probabilities: reachability, the partition into communicating
-- classes, irreducibility, the period of each state/class, and the
-- recurrence/transience of each state (for a finite chain, recurrent means
-- exactly that the communicating class is closed).
--
-- The graph combinatorics live in "Dtmc.Internal.Graph"; this module is the thin
-- bridge that reads a chain's support into a 'Graph' and renames the graph facts
-- into Markov-chain vocabulary: communicating classes are the strongly connected
-- components, and a class's period is the period of that component.
--
-- The API has two tiers. 'supportGraphOf' builds a 'SupportGraph' -- the
-- support relation with its reachability closure precomputed, an @O(n^3)@
-- construction -- and the @..In@ variants answer any number of queries against
-- it cheaply. The one-shot functions ('accessible', 'communicates', 'period',
-- ...) are conveniences that rebuild the graph on every call: fine for a
-- single question, wasteful in a loop.
module Dtmc.Classification (
    -- * Precomputed support graph
    SupportGraph,
    supportGraphOf,
    supportEdgeIn,
    accessibleIn,
    communicatesIn,
    communicatingClassesIn,
    irreducibleIn,
    periodIn,
    aperiodicIn,
    recurrentStateIn,
    transientStateIn,
    recurrentStatesIn,
    transientStatesIn,
    classifyIn,

    -- * One-shot support-graph queries
    supportEdge,
    accessible,
    communicates,

    -- * Communicating classes
    communicatingClasses,
    irreducible,

    -- * Periodicity
    period,
    aperiodic,

    -- * Recurrence and transience
    recurrentState,
    transientState,
    recurrentStates,
    transientStates,

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
    Nat,
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

-- | The support graph of a chain with its reachability closure precomputed.
-- Build it once with 'supportGraphOf' (@O(n^3)@), then run any number of
-- queries against it via the @..In@ functions. The phantom @n@ ties the
-- 'Finite' state indices to the graph's dimension, just as for
-- 'Dtmc.Distribution.Distribution'.
newtype SupportGraph (n :: Nat) = SupportGraph Graph

type role SupportGraph nominal

-- | Build the 'SupportGraph' of a chain: the single @O(n^3)@ closure
-- computation that every query in this module ultimately reads from.
supportGraphOf :: (KnownNat n) => TransitionMatrix n -> SupportGraph n
supportGraphOf = SupportGraph . supportGraph

-- Convert a raw @Int@ state index into the bounded 'Finite' @n@ index.
toFinite :: (KnownNat n) => Int -> Finite n
toFinite = finite . fromIntegral

-- Convert a bounded 'Finite' @n@ index back to a raw @Int@.
toIndex :: Finite n -> Int
toIndex = fromIntegral . getFinite

-- | Direct one-step reachability in a prebuilt graph: @True@ iff @P(i,j) > 0@.
supportEdgeIn :: SupportGraph n -> Finite n -> Finite n -> Bool
supportEdgeIn (SupportGraph g) i j = edge g (toIndex i) (toIndex j)

-- | Accessibility @i -> j@ in a prebuilt graph: @j@ is reachable from @i@ in
-- zero or more steps.
accessibleIn :: SupportGraph n -> Finite n -> Finite n -> Bool
accessibleIn (SupportGraph g) i j = reachable g (toIndex i) (toIndex j)

-- | Communication @i <-> j@ in a prebuilt graph: mutual accessibility.
communicatesIn :: SupportGraph n -> Finite n -> Finite n -> Bool
communicatesIn (SupportGraph g) i j =
    reachable g a b && reachable g b a
  where
    a = toIndex i
    b = toIndex j

-- | The communicating classes of a prebuilt graph, each as a list of states.
communicatingClassesIn :: (KnownNat n) => SupportGraph n -> [[Finite n]]
communicatingClassesIn (SupportGraph g) =
    map (map toFinite) (components g)

-- | Irreducibility of a prebuilt graph: all states form a single (non-empty)
-- communicating class.
irreducibleIn :: SupportGraph n -> Bool
irreducibleIn (SupportGraph g) =
    case components g of
        [c] -> not (null c)
        _ -> False

-- | Period of state @i@ in a prebuilt graph (see 'period').
periodIn :: SupportGraph n -> Finite n -> Maybe Natural
periodIn (SupportGraph g) i =
    componentPeriod g (componentOf g (toIndex i))

-- | Aperiodicity of a prebuilt graph (see 'aperiodic').
aperiodicIn :: SupportGraph n -> Bool
aperiodicIn (SupportGraph g) =
    not (null cs) && all ((== Just 1) . componentPeriod g) cs
  where
    cs = components g

-- | Recurrence of state @i@ in a prebuilt graph (see 'recurrentState').
recurrentStateIn :: SupportGraph n -> Finite n -> Bool
recurrentStateIn (SupportGraph g) i =
    closed g (componentOf g (toIndex i))

-- | Transience of state @i@ in a prebuilt graph (see 'transientState').
transientStateIn :: SupportGraph n -> Finite n -> Bool
transientStateIn sg i = not (recurrentStateIn sg i)

-- | All recurrent states of a prebuilt graph (see 'recurrentStates').
recurrentStatesIn :: (KnownNat n) => SupportGraph n -> [Finite n]
recurrentStatesIn (SupportGraph g) =
    concatMap (map toFinite) (filter (closed g) (components g))

-- | All transient states of a prebuilt graph (see 'transientStates').
transientStatesIn :: (KnownNat n) => SupportGraph n -> [Finite n]
transientStatesIn (SupportGraph g) =
    concatMap (map toFinite) (filter (not . closed g) (components g))

-- | Direct one-step reachability: @True@ iff @P(i,j) > 0@. Rebuilds the
-- support graph on every call; for repeated queries build one 'SupportGraph'
-- with 'supportGraphOf' and use 'supportEdgeIn'.
supportEdge :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
supportEdge p = supportEdgeIn (supportGraphOf p)

-- | Accessibility @i -> j@: @j@ is reachable from @i@ in zero or more steps.
-- One-shot; see 'accessibleIn' for the amortised variant.
accessible :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
accessible p = accessibleIn (supportGraphOf p)

-- | Communication @i <-> j@: @i@ and @j@ are mutually accessible. This is the
-- equivalence relation whose classes are the communicating classes.
-- One-shot; see 'communicatesIn' for the amortised variant.
communicates :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
communicates p = communicatesIn (supportGraphOf p)

-- | The communicating classes of the chain, each as a list of states.
-- One-shot; see 'communicatingClassesIn' for the amortised variant.
communicatingClasses :: (KnownNat n) => TransitionMatrix n -> [[Finite n]]
communicatingClasses = communicatingClassesIn . supportGraphOf

-- | Whether the chain is irreducible: all states form a single (non-empty)
-- communicating class, so every state is reachable from every other.
-- One-shot; see 'irreducibleIn' for the amortised variant.
irreducible :: (KnownNat n) => TransitionMatrix n -> Bool
irreducible = irreducibleIn . supportGraphOf

-- | Period of state @i@: the gcd of the lengths of all closed walks through
-- @i@, computed within @i@'s communicating class. 'Nothing' when the class has
-- no cycles (the period is undefined).
-- One-shot; see 'periodIn' for the amortised variant.
period :: (KnownNat n) => TransitionMatrix n -> Finite n -> Maybe Natural
period p = periodIn (supportGraphOf p)

-- | Whether the chain is aperiodic: every communicating class has period @1@
-- (and there is at least one state). A class with undefined period -- a
-- cycle-free transient class, where 'period' is 'Nothing' -- does not have
-- period one, so its presence makes the chain non-aperiodic under this
-- definition. One-shot; see 'aperiodicIn' for the amortised variant.
aperiodic :: (KnownNat n) => TransitionMatrix n -> Bool
aperiodic = aperiodicIn . supportGraphOf

-- | Whether state @i@ is recurrent: started at @i@, the chain returns to @i@
-- with probability one. For a /finite/ chain this is purely combinatorial:
-- A state of a finite chain is recurrent iff its communicating
-- class is closed.
--
-- The equivalence is specific to finite chains: on an infinite state space a
-- closed class may be transient (e.g. the asymmetric random walk on the
-- integers), so this must not be read as a statement about infinite chains
-- truncated to finite matrices.
-- One-shot; see 'recurrentStateIn' for the amortised variant.
recurrentState :: (KnownNat n) => TransitionMatrix n -> Finite n -> Bool
recurrentState p = recurrentStateIn (supportGraphOf p)

-- | Whether state @i@ is transient: positive probability of never returning;
-- the negation of 'recurrentState'. One-shot; see 'transientStateIn' for the
-- amortised variant.
transientState :: (KnownNat n) => TransitionMatrix n -> Finite n -> Bool
transientState p = transientStateIn (supportGraphOf p)

-- | All recurrent states: the members of the closed communicating classes, in
-- the order 'classify' lists them. Never empty: the classes of a finite chain
-- form an acyclic reachability digraph (a cycle of classes would merge them
-- into one class), so some class has no outgoing edges, and a sink class is
-- closed. One-shot; see 'recurrentStatesIn' for the amortised variant.
recurrentStates :: (KnownNat n) => TransitionMatrix n -> [Finite n]
recurrentStates = recurrentStatesIn . supportGraphOf

-- | All transient states: the members of the non-closed communicating
-- classes, in the order 'classify' lists them. Empty iff every class is
-- closed (in particular for any irreducible chain). One-shot; see
-- 'transientStatesIn' for the amortised variant.
transientStates :: (KnownNat n) => TransitionMatrix n -> [Finite n]
transientStates = transientStatesIn . supportGraphOf

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

-- | Decompose a prebuilt graph into communicating classes, annotating each
-- with its period and whether it is closed.
classifyIn :: (KnownNat n) => SupportGraph n -> Classification n
classifyIn (SupportGraph g) =
    Classification
        [ CommClass
            { classMembers = map toFinite c
            , classPeriod = componentPeriod g c
            , classClosed = closed g c
            }
        | c <- components g
        ]

-- | Decompose a chain into communicating classes, annotating each with its
-- period and whether it is closed. One-shot; see 'classifyIn' for the
-- amortised variant.
classify :: (KnownNat n) => TransitionMatrix n -> Classification n
classify = classifyIn . supportGraphOf

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
