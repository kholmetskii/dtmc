{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.ClassificationSpec (
    spec,
) where

import Data.Finite (
    Finite,
    finites,
    getFinite,
 )
import Data.List (
    sort,
 )
import Dtmc.Analysis.Classification (
    CommClass (..),
    absorbingStates,
    accessible,
    aperiodic,
    chainPeriod,
    communicates,
    communicatingClasses,
    cyclicClasses,
    ergodic,
    irreducible,
    period,
    reachesAny,
    recurrentState,
    recurrentStates,
    supportEdge,
    transientState,
    transientStates,
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.TestSupport (
    chunksOf,
    genTransitionRows,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    TransitionMatrixError,
    fromRows,
    toRows,
 )
import GHC.Generics (
    Generic,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.Natural (Natural)
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
 )
import Test.Hspec.QuickCheck (
    prop,
 )
import Test.QuickCheck (
    Property,
    conjoin,
    counterexample,
    forAll,
    property,
    (===),
 )

data NamedClassState = ClassA | ClassB | ClassC
    deriving (Eq, Ord, Show, Generic)

instance FiniteState NamedClassState

checked :: (Show e) => Either e a -> a
checked = either (error . show) id

threeCycle :: TransitionMatrix (Finite 3)
threeCycle =
    checked $
        fromRows
            ( chunksOf
                3
                [ 0
                , 1
                , 0
                , 0
                , 0
                , 1
                , 1
                , 0
                , 0
                ]
            )

namedThreeCycle :: TransitionMatrix NamedClassState
namedThreeCycle =
    checked $
        fromRows @NamedClassState
            (chunksOf 3 [0, 1, 0, 0, 0, 1, 1, 0, 0])

selfLoopTwo :: TransitionMatrix (Finite 2)
selfLoopTwo =
    checked $
        fromRows
            ( chunksOf
                2
                [ 0.5
                , 0.5
                , 1.0
                , 0.0
                ]
            )

bipartiteTwo :: TransitionMatrix (Finite 2)
bipartiteTwo =
    checked $
        fromRows
            ( chunksOf
                2
                [ 0
                , 1
                , 1
                , 0
                ]
            )

sevenState :: TransitionMatrix (Finite 7)
sevenState =
    checked $
        fromRows
            ( chunksOf
                7
                [ 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0.4
                , 0
                , 0.6
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0.3
                , 0
                , 0.7
                , 0
                , 0
                , 0
                , 0
                , 0.3
                , 0.4
                , 0
                , 0.3
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0.2
                , 0
                , 0.8
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1
                ]
            )

identityThree :: TransitionMatrix (Finite 3)
identityThree =
    checked $
        fromRows
            ( chunksOf
                3
                [ 1
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 1
                ]
            )

-- Exercise 3.2.2: irreducible, period 2, cyclic classes {A,B} and {C,D}.
fourStateCyclic :: TransitionMatrix (Finite 4)
fourStateCyclic =
    checked $
        fromRows
            ( chunksOf
                4
                [ 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 1
                , 0.5
                , 0.5
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                ]
            )

matrixSupport :: TransitionMatrix (Finite n) -> [[Bool]]
matrixSupport = map (map (> 0)) . toRows

boolMul :: [[Bool]] -> [[Bool]] -> [[Bool]]
boolMul a b =
    [ [or [ai && (b !! k !! j) | (k, ai) <- zip [0 ..] row] | j <- idxs]
    | row <- a
    ]
  where
    idxs = [0 .. length a - 1]

boolIdentity :: Int -> [[Bool]]
boolIdentity dim = [[i == j | j <- [0 .. dim - 1]] | i <- [0 .. dim - 1]]

referencePeriod :: [[Bool]] -> Int -> Maybe Natural
referencePeriod s i =
    case returns of
        [] -> Nothing
        _ -> Just (fromIntegral (foldr1 gcd returns))
  where
    dim = length s
    bound = 4 * dim * dim + 1
    powers = take bound (drop 1 (iterate (boolMul s) (boolIdentity dim)))
    returns = [k | (k, m) <- zip [1 :: Int ..] powers, (m !! i) !! i]

classesAsInts :: TransitionMatrix (Finite n) -> [[Integer]]
classesAsInts = map (map getFinite . classMembers) . communicatingClasses

cyclicClassesAsInts :: (KnownNat n) => TransitionMatrix (Finite n) -> Maybe [[Integer]]
cyclicClassesAsInts = fmap (map (map getFinite)) . cyclicClasses

sortUnique :: (Ord a) => [a] -> [a]
sortUnique = foldr insert []
  where
    insert x [] = [x]
    insert x (y : ys)
        | x < y = x : y : ys
        | x == y = y : ys
        | otherwise = y : insert x ys

periodMatchesReference :: (KnownNat n) => TransitionMatrix (Finite n) -> [Finite n] -> Property
periodMatchesReference p states =
    conjoin
        [ period p i === referencePeriod s (fromIntegral (getFinite i))
        | i <- states
        ]
  where
    s = matrixSupport p

spec :: Spec
spec = do
    describe "communication is an equivalence relation" $ do
        prop "is reflexive, symmetric, and transitive on random support graphs" $
            forAll (genTransitionRows 4) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 4)) of
                    Right p ->
                        let states = finites :: [Finite 4]
                         in conjoin
                                [ conjoin
                                    [ counterexample "reflexivity" (communicates p i i)
                                    | i <- states
                                    ]
                                , conjoin
                                    [ counterexample "symmetry" $
                                        communicates p i j === communicates p j i
                                    | i <- states
                                    , j <- states
                                    ]
                                , conjoin
                                    [ counterexample "transitivity" $
                                        not (communicates p i j && communicates p j k)
                                            || communicates p i k
                                    | i <- states
                                    , j <- states
                                    , k <- states
                                    ]
                                ]
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

        prop "accessibility is reflexive" $
            forAll (genTransitionRows 4) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 4)) of
                    Right p ->
                        conjoin
                            [ property (accessible p i i)
                            | i <- finites :: [Finite 4]
                            ]
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

    describe "reachesAny" $ do
        it "finds a reachable target" $
            reachesAny threeCycle 0 [2] `shouldBe` True

        it "returns false for an empty target set" $
            reachesAny threeCycle 0 [] `shouldBe` False

        it "uses zero-step reachability" $
            reachesAny identityThree 1 [1] `shouldBe` True

    describe "period" $ do
        it "is 3 for every state of the three-cycle" $
            map (period threeCycle) (finites :: [Finite 3])
                `shouldBe` [Just 3, Just 3, Just 3]

        it "is 1 for the self-loop chain (aperiodic)" $ do
            map (period selfLoopTwo) (finites :: [Finite 2])
                `shouldBe` [Just 1, Just 1]
            aperiodic selfLoopTwo `shouldBe` True

        it "is 2 for the bipartite swap (periodic)" $ do
            map (period bipartiteTwo) (finites :: [Finite 2])
                `shouldBe` [Just 2, Just 2]
            aperiodic bipartiteTwo `shouldBe` False

        it "matches the hand-computed periods of the seven-state chain" $
            map (period sevenState) (finites :: [Finite 7])
                `shouldBe` [Just 2, Just 2, Just 1, Just 1, Just 1, Just 1, Just 1]

        prop "agrees with the gcd of return-time lengths (random @4)" $
            forAll (genTransitionRows 4) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 4)) of
                    Right p -> periodMatchesReference p (finites :: [Finite 4])
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

        prop "agrees with the gcd of return-time lengths (random @3)" $
            forAll (genTransitionRows 3) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 3)) of
                    Right p -> periodMatchesReference p (finites :: [Finite 3])
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

    describe "communicatingClasses" $ do
        it "splits the seven-state chain into {A,B}, {C,D,E,F}, {G}" $
            classesAsInts sevenState `shouldBe` [[0, 1], [2, 3, 4, 5], [6]]

        it "returns a single class for the irreducible three-cycle" $
            classesAsInts threeCycle `shouldBe` [[0, 1, 2]]

        prop "the classes partition the state space (random @4)" $
            forAll (genTransitionRows 4) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 4)) of
                    Right p ->
                        property (sortUnique (concat (classesAsInts p)) == [0 .. 3])
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

        prop "communication agrees with the class partition (random @4)" $
            forAll (genTransitionRows 4) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 4)) of
                    Right p ->
                        let states = finites :: [Finite 4]
                            classIx = map classMembers (communicatingClasses p)
                            sameClass i j = or [i `elem` c && j `elem` c | c <- classIx]
                         in conjoin
                                [ counterexample (show (i, j)) $
                                    communicates p i j === sameClass i j
                                | i <- states
                                , j <- states
                                ]
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

    describe "cyclicClasses" $ do
        it "splits the period-2 four-state chain into {A,B} and {C,D}" $
            cyclicClassesAsInts fourStateCyclic `shouldBe` Just [[0, 1], [2, 3]]

        it "splits the three-cycle into three singletons" $
            cyclicClassesAsInts threeCycle `shouldBe` Just [[0], [1], [2]]

        it "is Nothing for the reducible seven-state chain" $
            cyclicClassesAsInts sevenState `shouldBe` Nothing

        prop "classes partition the states and advance one step (random @4)" $
            forAll (genTransitionRows 4) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 4)) of
                    Right p ->
                        case cyclicClasses p of
                            Nothing -> property True
                            Just cs ->
                                let d = length cs
                                    states = finites :: [Finite 4]
                                 in conjoin
                                        [ counterexample "partition" (sort (concat cs) === states)
                                        , conjoin
                                            [ counterexample (show (i, j)) $
                                                property (j `elem` (cs !! ((r + 1) `mod` d)))
                                            | (r, c) <- zip [0 ..] cs
                                            , i <- c
                                            , j <- states
                                            , supportEdge p i j
                                            ]
                                        ]
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

    describe "irreducible" $ do
        it "holds for the three-cycle and swap, fails for the seven-state chain" $ do
            irreducible threeCycle `shouldBe` True
            irreducible bipartiteTwo `shouldBe` True
            irreducible sevenState `shouldBe` False

    describe "communicatingClasses details" $ do
        it "records members, periods, and closedness for the seven-state chain" $ do
            let cs = communicatingClasses sevenState
            map (map getFinite . classMembers) cs
                `shouldBe` [[0, 1], [2, 3, 4, 5], [6]]
            map classPeriod cs `shouldBe` [Just 2, Just 1, Just 1]
            map classClosed cs `shouldBe` [True, False, True]

    describe "absorbingStates" $ do
        it "finds the absorbing states" $ do
            map getFinite (absorbingStates sevenState) `shouldBe` [6]
            map getFinite (absorbingStates identityThree) `shouldBe` [0, 1, 2]
            map getFinite (absorbingStates threeCycle) `shouldBe` []

        prop "absorbing states have only a self-loop (random @4)" $
            forAll (genTransitionRows 4) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 4)) of
                    Right p ->
                        conjoin
                            [ counterexample (show i) $
                                [j | j <- finites :: [Finite 4], supportEdge p i j] === [i]
                            | i <- absorbingStates p
                            ]
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

    describe "whole-chain queries agree with the class summaries" $ do
        -- These are not restatements of one definition: the left-hand sides
        -- reach the support graph through G.components and G.componentPeriod,
        -- the right-hand sides through G.periodOf and per-class closedness.
        prop "on random @4 chains" $
            forAll (genTransitionRows 4) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 4)) of
                    Right p ->
                        let cs = communicatingClasses p
                            closed = filter classClosed cs
                            open = filter (not . classClosed) cs
                         in conjoin
                                [ counterexample "irreducible" $
                                    irreducible p === (length cs == 1)
                                , counterexample "aperiodic" $
                                    aperiodic p
                                        === (not (null cs) && all ((== Just 1) . classPeriod) cs)
                                , counterexample "ergodic" $
                                    ergodic p === (irreducible p && aperiodic p)
                                , counterexample "recurrentStates" $
                                    recurrentStates p
                                        === concatMap classMembers closed
                                , counterexample "transientStates" $
                                    transientStates p
                                        === concatMap classMembers open
                                , counterexample "chainPeriod" $
                                    chainPeriod p
                                        === case cs of
                                            [singleClass] -> classPeriod singleClass
                                            _ -> Nothing
                                ]
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

        prop "chainPeriod is the shared period of an irreducible chain (random @4)" $
            forAll (genTransitionRows 4) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 4)) of
                    Right p
                        | irreducible p ->
                            conjoin
                                [ counterexample (show i) (chainPeriod p === period p i)
                                | i <- finites :: [Finite 4]
                                ]
                        | otherwise -> property True
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

    describe "recurrence and transience" $ do
        it "matches the closed classes of the seven-state chain" $ do
            map getFinite (recurrentStates sevenState) `shouldBe` [0, 1, 6]
            map getFinite (transientStates sevenState) `shouldBe` [2, 3, 4, 5]

        it "marks every state of the irreducible three-cycle recurrent" $ do
            map getFinite (recurrentStates threeCycle) `shouldBe` [0, 1, 2]
            transientStates threeCycle `shouldBe` []

        it "marks every state of the identity chain recurrent" $ do
            map getFinite (recurrentStates identityThree) `shouldBe` [0, 1, 2]
            transientStates identityThree `shouldBe` []

        prop "recurrent and transient states partition the state space (random @4)" $
            forAll (genTransitionRows 4) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 4)) of
                    Right p ->
                        sort
                            ( map getFinite (recurrentStates p)
                                <> map getFinite (transientStates p)
                            )
                            === [0 .. 3]
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

        prop "every finite chain has a recurrent state (random @4)" $
            forAll (genTransitionRows 4) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 4)) of
                    Right p ->
                        property (not (null (recurrentStates p)))
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

        prop "transient iff some reachable state cannot reach back (random @4)" $
            forAll (genTransitionRows 4) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 4)) of
                    Right p ->
                        let states = finites :: [Finite 4]
                         in conjoin
                                [ transientState p i
                                    === or
                                        [ accessible p i j && not (accessible p j i)
                                        | j <- states
                                        ]
                                | i <- states
                                ]
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

        prop "predicates agree with the state lists (random @4)" $
            forAll (genTransitionRows 4) $ \matrix ->
                case fromRows matrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 4)) of
                    Right p ->
                        let states = finites :: [Finite 4]
                         in conjoin
                                [ recurrentState p i === (i `elem` recurrentStates p)
                                | i <- states
                                ]
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

    describe "named finite states" $ do
        it "reports communication and periods with named constructors" $ do
            map classMembers (communicatingClasses namedThreeCycle)
                `shouldBe` [[ClassA, ClassB, ClassC]]
            map (period namedThreeCycle) [ClassA, ClassB, ClassC]
                `shouldBe` replicate 3 (Just 3)

        it "returns named recurrent states in canonical order" $
            recurrentStates namedThreeCycle
                `shouldBe` [ClassA, ClassB, ClassC]

        it "answers whole-chain queries with named constructors" $ do
            map classMembers (communicatingClasses namedThreeCycle)
                `shouldBe` [[ClassA, ClassB, ClassC]]
            recurrentStates namedThreeCycle `shouldBe` [ClassA, ClassB, ClassC]
            absorbingStates namedThreeCycle `shouldBe` []
            chainPeriod namedThreeCycle `shouldBe` Just 3
            ergodic namedThreeCycle `shouldBe` False
