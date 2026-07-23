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
-- bridge that reads a chain's support into a graph and renames the graph facts
-- into Markov-chain vocabulary: communicating classes are the strongly connected
-- components, and a class's period is the period of that component.
--
-- Every function takes a 'TransitionMatrix' and reads its support graph, which
-- the matrix carries as a lazy field ('Dtmc.Internal.Types.tmSupport') built
-- once on first use. Repeated queries on the /same/ matrix value therefore share
-- that single build automatically -- there is no separate prebuilt-graph object
-- to construct and thread.
module Dtmc.Classification (
    -- * Reachability
    supportEdge,
    accessible,
    reachesAny,
    backwardReachable,
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
    CommClass (..),
    Classification,
    classesOf,
    isIrreducible,
    isAperiodic,
    isErgodic,
    chainPeriod,
    recurrentStatesOf,
    transientStatesOf,
    absorbingStates,
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
import Dtmc.Internal.Graph qualified as G
import Dtmc.Internal.Types (TransitionMatrix, tmSupport)
import GHC.TypeNats (KnownNat)
import Numeric.Natural (Natural)

-- Convert a raw @Int@ state index into the bounded 'Finite' @n@ index.
toFinite :: (KnownNat n) => Int -> Finite n
toFinite = finite . fromIntegral

-- Convert a bounded 'Finite' @n@ index back to a raw @Int@.
toIndex :: Finite n -> Int
toIndex = fromIntegral . getFinite

-- | Direct one-step reachability: @True@ iff @P(i,j) > 0@.
supportEdge :: TransitionMatrix n -> Finite n -> Finite n -> Bool
supportEdge p i j = G.hasEdge (tmSupport p) (toIndex i) (toIndex j)

-- | Accessibility @i -> j@: @j@ is reachable from @i@ in zero or more steps.
accessible :: TransitionMatrix n -> Finite n -> Finite n -> Bool
accessible p i j = G.reachable (tmSupport p) (toIndex i) (toIndex j)

-- | Whether state @i@ can reach /any/ of the given targets in zero or more
-- steps. Expands @i@'s reachable set once and tests every target against it,
-- rather than searching afresh per target as @'any' ('accessible' p i)@ would.
reachesAny :: TransitionMatrix n -> Finite n -> [Finite n] -> Bool
reachesAny p i targets =
    G.reachesAny (tmSupport p) (toIndex i) (map toIndex targets)

-- | Reverse reachability within an induced subgraph: the states from which some
-- @seed@ is reachable along a directed support path whose states /all/ satisfy
-- @allowed@ (seeds included, when allowed). Unlike 'accessible', which ranges
-- over the whole graph, the path here must stay inside @allowed@.
backwardReachable ::
    (KnownNat n) =>
    TransitionMatrix n ->
    (Finite n -> Bool) ->
    [Finite n] ->
    [Finite n]
backwardReachable p allowed seeds =
    map toFinite (G.backwardReachable (tmSupport p) (allowed . toFinite) (map toIndex seeds))

-- | Communication @i <-> j@: @i@ and @j@ are mutually accessible. This is the
-- equivalence relation whose classes are the communicating classes.
communicates :: TransitionMatrix n -> Finite n -> Finite n -> Bool
communicates p i j =
    G.reachable g a b && G.reachable g b a
  where
    g = tmSupport p
    a = toIndex i
    b = toIndex j

-- | The communicating classes of the chain, each as a list of states: the
-- strongly connected components of the support graph, ordered by least member.
communicatingClasses :: (KnownNat n) => TransitionMatrix n -> [[Finite n]]
communicatingClasses p = map (map toFinite) (G.components (tmSupport p))

-- | Whether the chain is irreducible: all states form a single (non-empty)
-- communicating class, so every state is reachable from every other.
irreducible :: TransitionMatrix n -> Bool
irreducible p =
    case G.components (tmSupport p) of
        [c] -> not (null c)
        _ -> False

-- | Period of state @i@: the gcd of the lengths of all closed walks through
-- @i@, computed within @i@'s communicating class. 'Nothing' when the class has
-- no cycles (a single state with no self-loop), where the period is undefined.
period :: TransitionMatrix n -> Finite n -> Maybe Natural
period p i = G.periodOf (tmSupport p) (toIndex i)

-- | Whether the chain is aperiodic: every communicating class has period @1@
-- (and there is at least one state). A class with undefined period -- a
-- cycle-free transient class, where 'period' is 'Nothing' -- does not have
-- period one, so its presence makes the chain non-aperiodic under this
-- definition.
aperiodic :: TransitionMatrix n -> Bool
aperiodic p =
    not (null cs) && all ((== Just 1) . G.componentPeriod g) cs
  where
    g = tmSupport p
    cs = G.components g

-- | The cyclic (periodicity) classes of an /irreducible/ chain, in cyclic order
-- @C_0, C_1, ..., C_{d-1}@ where @d@ is the period: from any state in @C_r@,
-- every one-step move lands in @C_{(r+1) mod d}@. 'Nothing' when the chain is
-- not irreducible or its period is undefined (a single state with no self-loop).
--
-- States are grouped by their phase within the (single) communicating class, so
-- this is @O(n)@ on top of the shared support graph.
cyclicClasses :: (KnownNat n) => TransitionMatrix n -> Maybe [[Finite n]]
cyclicClasses p
    | not (irreducible p) = Nothing
    | otherwise =
        case G.periodOf g 0 of
            Nothing -> Nothing
            Just d ->
                Just
                    [ [toFinite v | v <- [0 .. G.graphDim g - 1], G.phaseOf g v == r]
                    | r <- [0 .. fromIntegral d - 1]
                    ]
  where
    g = tmSupport p

-- | Whether state @i@ is recurrent: started at @i@, the chain returns to @i@
-- with probability one. For a /finite/ chain this is purely combinatorial: a
-- state of a finite chain is recurrent iff its communicating class is closed.
--
-- The equivalence is specific to finite chains: on an infinite state space a
-- closed class may be transient (e.g. the asymmetric random walk on the
-- integers), so this must not be read as a statement about infinite chains
-- truncated to finite matrices.
recurrentState :: TransitionMatrix n -> Finite n -> Bool
recurrentState p i = G.inClosedComponent (tmSupport p) (toIndex i)

-- | Whether state @i@ is transient: positive probability of never returning;
-- the negation of 'recurrentState'.
transientState :: TransitionMatrix n -> Finite n -> Bool
transientState p i = not (recurrentState p i)

-- | All recurrent states: the members of the closed communicating classes, in
-- the order 'classify' lists them. Never empty: the classes of a finite chain
-- form an acyclic reachability digraph (a cycle of classes would merge them into
-- one class), so some class has no outgoing edges, and a sink class is closed.
recurrentStates :: (KnownNat n) => TransitionMatrix n -> [Finite n]
recurrentStates p =
    concatMap (map toFinite) (filter (G.isClosed g) (G.components g))
  where
    g = tmSupport p

-- | All transient states: the members of the non-closed communicating classes,
-- in the order 'classify' lists them. Empty iff every class is closed (in
-- particular for any irreducible chain).
transientStates :: (KnownNat n) => TransitionMatrix n -> [Finite n]
transientStates p =
    concatMap (map toFinite) (filter (not . G.isClosed g) (G.components g))
  where
    g = tmSupport p

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

-- | A full qualitative report of a chain: its communicating classes together
-- with the chain-level facts derived from them. Built only by 'classify' (the
-- constructor is hidden), so a 'Classification' is always internally consistent.
-- The fact fields are plain projections, so reading them carries no 'KnownNat'
-- constraint -- all the type-level work happens once, in 'classify'.
data Classification n = Classification
    { classesOf :: [CommClass n]
    -- ^ The communicating classes, ordered by least member.
    , isIrreducible :: Bool
    -- ^ Whether the states form a single (non-empty) communicating class.
    , isAperiodic :: Bool
    -- ^ Whether every class has period @1@ (and there is at least one class).
    , isErgodic :: Bool
    -- ^ Whether the chain is ergodic -- irreducible and aperiodic -- the
    -- condition under which it converges to a unique stationary distribution.
    , chainPeriod :: Maybe Natural
    -- ^ The period of the chain when it is irreducible (@Just d@); @Nothing@ for
    -- a reducible chain (where the period is a per-class notion) or when the
    -- single class has no cycles.
    , recurrentStatesOf :: [Finite n]
    -- ^ States lying in closed classes -- recurrent, in the finite-chain sense.
    , transientStatesOf :: [Finite n]
    -- ^ States lying in non-closed classes -- transient.
    , absorbingStates :: [Finite n]
    -- ^ Absorbing states: those forming a closed singleton class, i.e. states
    -- @i@ with @P(i,i) = 1@.
    }

type role Classification nominal

deriving instance (KnownNat n) => Eq (Classification n)

deriving instance (KnownNat n) => Show (Classification n)

-- | Decompose a chain into its communicating classes and summarise the
-- chain-level structure. Every field is computed from a single shared support
-- graph, so this one call answers "classes, irreducibility, aperiodicity,
-- recurrent and transient states" together.
classify :: (KnownNat n) => TransitionMatrix n -> Classification n
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
            { classMembers = map toFinite c
            , classPeriod = G.periodOf g v
            , classClosed = G.inClosedComponent g v
            }
        | c@(v : _) <- G.components g
        ]
    irreducible' = case cs of
        [_] -> True
        _ -> False
    aperiodic' = not (null cs) && all ((== Just 1) . classPeriod) cs
    chainPeriod' = case cs of
        [c] -> classPeriod c
        _ -> Nothing

-- | A transition matrix carrying a proof that it is irreducible, obtainable
-- only via 'witnessIrreducible'. Lets downstream code demand irreducibility in
-- a type rather than re-checking it.
newtype Irreducible n = Irreducible (TransitionMatrix n)

type role Irreducible nominal

deriving instance (KnownNat n) => Show (Irreducible n)

-- | Certify irreducibility: @Just@ the wrapped matrix when 'irreducible' holds,
-- @Nothing@ otherwise.
witnessIrreducible :: TransitionMatrix n -> Maybe (Irreducible n)
witnessIrreducible p
    | irreducible p = Just (Irreducible p)
    | otherwise = Nothing

-- | Recover the underlying transition matrix from an 'Irreducible' witness.
irreducibleMatrix :: Irreducible n -> TransitionMatrix n
irreducibleMatrix (Irreducible p) = p
