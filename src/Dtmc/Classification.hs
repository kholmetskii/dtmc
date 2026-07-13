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

supportMatrix :: (KnownNat n) => TransitionMatrix n -> [[Bool]]
supportMatrix p =
    map (map (> 0)) (LA.toLists (S.extract (unTransitionMatrix p)))

at :: [[Bool]] -> Int -> Int -> Bool
at m i j = (m !! i) !! j

toFinite :: (KnownNat n) => Int -> Finite n
toFinite = finite . fromIntegral

toIndex :: Finite n -> Int
toIndex = fromIntegral . getFinite

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

supportEdge :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
supportEdge p i j = at (supportMatrix p) (toIndex i) (toIndex j)

accessible :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
accessible p i j = at (reachClosure (supportMatrix p)) (toIndex i) (toIndex j)

communicates :: (KnownNat n) => TransitionMatrix n -> Finite n -> Finite n -> Bool
communicates p i j =
    at r a b && at r b a
  where
    r = reachClosure (supportMatrix p)
    a = toIndex i
    b = toIndex j

rawClasses :: [[Bool]] -> [[Int]]
rawClasses r = go [0 .. length r - 1]
  where
    comm i j = at r i j && at r j i
    go [] = []
    go (x : xs) =
        (x : filter (comm x) xs) : go (filter (not . comm x) xs)

communicatingClasses :: (KnownNat n) => TransitionMatrix n -> [[Finite n]]
communicatingClasses p =
    map (map toFinite) (rawClasses (reachClosure (supportMatrix p)))

irreducible :: (KnownNat n) => TransitionMatrix n -> Bool
irreducible p =
    case rawClasses (reachClosure (supportMatrix p)) of
        [c] -> not (null c)
        _ -> False


period :: (KnownNat n) => TransitionMatrix n -> Finite n -> Maybe Natural
period p i =
    periodOfClass s klass
  where
    s = supportMatrix p
    r = reachClosure s
    a = toIndex i
    klass = filter (\j -> at r a j && at r j a) [0 .. length s - 1]

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

aperiodic :: forall n. (KnownNat n) => TransitionMatrix n -> Bool
aperiodic p =
    not (null states) && all (\i -> period p i == Just 1) states
  where
    states = finites :: [Finite n]

data CommClass n = CommClass
    { classMembers :: [Finite n]
    , classPeriod :: Maybe Natural
    , classClosed :: Bool
    }

deriving instance (KnownNat n) => Eq (CommClass n)

deriving instance (KnownNat n) => Show (CommClass n)

newtype Classification n = Classification [CommClass n]

type role Classification nominal

deriving instance (KnownNat n) => Eq (Classification n)

deriving instance (KnownNat n) => Show (Classification n)

classesOf :: Classification n -> [CommClass n]
classesOf (Classification cs) = cs

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

newtype Irreducible n = Irreducible (TransitionMatrix n)

type role Irreducible nominal

deriving instance (KnownNat n) => Show (Irreducible n)

witnessIrreducible :: (KnownNat n) => TransitionMatrix n -> Maybe (Irreducible n)
witnessIrreducible p
    | irreducible p = Just (Irreducible p)
    | otherwise = Nothing

irreducibleMatrix :: Irreducible n -> TransitionMatrix n
irreducibleMatrix (Irreducible p) = p
