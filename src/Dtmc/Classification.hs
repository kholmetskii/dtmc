{- |
Module      : Dtmc.Classification
Description : Communicating classes, irreducibility, periodicity, and recurrence.

Qualitative DTMC properties derived from the support graph of @P@: there is an
edge @i -> j@ exactly when the stored @P(i,j) > 0@. The comparison has no
tolerance, so a tiny positive value is a transition while zero or a negative
value is not. Results depend on which entries are positive, not their
magnitudes. Recurrence statements assume a finite, valid transition matrix.

For @n@ states and @E@ support edges, the first graph query scans @n^2@
entries and may spend @O(n log n + E)@ building the requested cached graph
facts. Queries on the same matrix share those lazy caches. A @0 x 0@ matrix
has no communicating classes and is neither irreducible nor aperiodic here.
-}
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
    unIrreducible,
) where

import Data.Array qualified as Array
import Data.Finite (
    Finite,
    finite,
    getFinite,
 )
import Dtmc.Classification.Internal (
    Classification (..),
    CommClass (..),
    Irreducible (Irreducible),
    unIrreducible,
 )
import Dtmc.Internal.Graph qualified as G
import Dtmc.Transition.Matrix.Internal (TransitionMatrix, tmSupport)
import GHC.TypeNats (KnownNat)
import Numeric.Natural (Natural)

toFinite :: (KnownNat n) => Int -> Finite n
toFinite = finite . fromIntegral

toIndex :: Finite n -> Int
toIndex = fromIntegral . getFinite

{- | Whether @P(i,j) > 0@: a direct transition in the support graph. No
tolerance is applied.

Time: @O(outDegree(i))@ after the support graph is built.
-}
supportEdge :: TransitionMatrix n -> Finite n -> Finite n -> Bool
supportEdge p i j = G.hasEdge (tmSupport p) (toIndex i) (toIndex j)

{- | Whether @j@ is reachable from @i@ in zero or more transitions. Hence every
state is reachable from itself, even without a self-transition.

Time: @O(n + E)@; traversal space: @O(n)@.
-}
accessible :: TransitionMatrix n -> Finite n -> Finite n -> Bool
accessible p i j = G.reachable (tmSupport p) (toIndex i) (toIndex j)

{- | Whether any target is reachable from @i@ in zero or more transitions.
An empty target list gives 'False'; including @i@ gives 'True'.

The graph is traversed once rather than once per target. Worst-case time:
@O(n + E + t)@ for @t@ supplied targets; space: @O(n)@.
-}
reachesAny :: TransitionMatrix n -> Finite n -> [Finite n] -> Bool
reachesAny p i targets =
    G.reachesAny (tmSupport p) (toIndex i) (map toIndex targets)

{- | States from which an allowed seed is reachable along a support path
containing only states accepted by @allowed@. Disallowed seeds are ignored;
the result is duplicate-free and ordered by state index.

Time: @O(n + E + s)@ plus @n@ predicate evaluations for @s@ seeds. Temporary
space: @O(n + E)@.
-}
backwardReachable ::
    (KnownNat n) =>
    TransitionMatrix n ->
    (Finite n -> Bool) ->
    [Finite n] ->
    [Finite n]
backwardReachable p allowed seeds =
    map toFinite (G.backwardReachable (tmSupport p) (allowed . toFinite) (map toIndex seeds))

{- | Whether @i@ and @j@ communicate: each state is reachable from the other.
This is an equivalence relation on the state space.

Time: @O(1)@ after strongly connected components are cached.
-}
communicates :: TransitionMatrix n -> Finite n -> Finite n -> Bool
communicates p i j =
    G.sameComponent (tmSupport p) (toIndex i) (toIndex j)

{- | The communicating classes, equivalent to the strongly connected
components of the support graph. States within a class and the classes by
least member are both in ascending order.

Result construction takes @O(n)@ after components are cached.
-}
communicatingClasses :: (KnownNat n) => TransitionMatrix n -> [[Finite n]]
communicatingClasses p = map (map toFinite) (G.components (tmSupport p))

{- | Whether every state communicates with every other state. The empty chain
is not irreducible.

Time: @O(1)@ after communicating classes are cached.
-}
irreducible :: TransitionMatrix n -> Bool
irreducible p =
    case G.components (tmSupport p) of
        [c] -> not (null c)
        _ -> False

{- | The period of @i@:
@gcd { k >= 1 | (P^k)(i,i) > 0 }@. Returns 'Nothing' when @i@ has no
positive-length return path, necessarily a singleton class without a
self-transition.

Time: @O(1)@ after periods are cached.
-}
period :: TransitionMatrix n -> Finite n -> Maybe Natural
period p i = G.periodOf (tmSupport p) (toIndex i)

{- | Whether every communicating class has period @1@. The empty chain and a
chain containing a class with undefined period are not aperiodic under this
definition.

Time: @O(c)@ for @c@ cached communicating classes.
-}
aperiodic :: TransitionMatrix n -> Bool
aperiodic p =
    not (null cs) && all ((== Just 1) . G.componentPeriod g) cs
  where
    g = tmSupport p
    cs = G.components g

{- | Partition an irreducible chain of period @d@ into cyclic classes
@C_0, ..., C_(d-1)@. Every transition from @C_r@ enters
@C_((r+1) mod d)@; @C_0@ contains the least state. Returns 'Nothing' for a
reducible chain or an undefined period.

Time and result space: @O(n)@ after graph facts are cached.
-}
cyclicClasses :: (KnownNat n) => TransitionMatrix n -> Maybe [[Finite n]]
cyclicClasses p
    | not (irreducible p) = Nothing
    | otherwise =
        case G.periodOf g 0 of
            Nothing -> Nothing
            Just d ->
                let dInt = fromIntegral d
                    -- One pass buckets every vertex by its phase (@O(V + d)@),
                    -- instead of scanning all vertices once per phase (@O(dV)@).
                    buckets =
                        Array.accumArray
                            (flip (:))
                            []
                            (0, dInt - 1)
                            [(G.phaseOf g v, toFinite v) | v <- [0 .. G.graphDim g - 1]]
                 in Just [reverse (buckets Array.! r) | r <- [0 .. dInt - 1]]
  where
    g = tmSupport p

{- | Whether the chain returns to @i@ with probability one when started there.
For a finite DTMC this holds exactly when @i@ belongs to a closed communicating
class.

Time: @O(1)@ after class closedness is cached.
-}
recurrentState :: TransitionMatrix n -> Finite n -> Bool
recurrentState p i = G.inClosedComponent (tmSupport p) (toIndex i)

{- | Whether @i@ is transient: the probability of returning from @i@ is less
than one. This is the negation of 'recurrentState' for a finite DTMC.

Time: @O(1)@ after class closedness is cached.
-}
transientState :: TransitionMatrix n -> Finite n -> Bool
transientState p i = not (recurrentState p i)

{- | Members of the closed communicating classes, ordered by class and state
index. Every non-empty finite DTMC has at least one; the empty chain returns
the empty list.

Time and result space: @O(n)@ after class closedness is cached.
-}
recurrentStates :: (KnownNat n) => TransitionMatrix n -> [Finite n]
recurrentStates p =
    concatMap (map toFinite) closedComponents
  where
    g = tmSupport p
    -- Closedness is read per component in O(1) from the cached table (via the
    -- representative vertex @v@), rather than recomputed by 'G.isClosed'.
    closedComponents =
        [component | component@(v : _) <- G.components g, G.inClosedComponent g v]

{- | Members of the non-closed communicating classes, ordered by class and
state index. The result is empty exactly when every class is closed.

Time and result space: @O(n)@ after class closedness is cached.
-}
transientStates :: (KnownNat n) => TransitionMatrix n -> [Finite n]
transientStates p =
    concatMap (map toFinite) openComponents
  where
    g = tmSupport p
    openComponents =
        [ component
        | component@(v : _) <- G.components g
        , not (G.inClosedComponent g v)
        ]

{- | Build all exported class, period, recurrence, absorbing-state,
irreducibility, and aperiodicity summaries from one shared support graph.

On an unforced matrix, time is @O(n^2 + n log n + E)@ and cached graph plus
report space is @O(n + E)@.
-}
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

{- | Return a witness exactly when the matrix is irreducible. The matrix is
wrapped unchanged.

Time: @O(1)@ after communicating classes are cached.
-}
witnessIrreducible :: TransitionMatrix n -> Maybe (Irreducible n)
witnessIrreducible p
    | irreducible p = Just (Irreducible p)
    | otherwise = Nothing
