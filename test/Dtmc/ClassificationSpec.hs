{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.ClassificationSpec (
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
import Dtmc.Classification (
    CommClass (..),
    absorbingStates,
    accessible,
    aperiodic,
    chainPeriod,
    classesOf,
    classify,
    communicates,
    communicatingClasses,
    cyclicClasses,
    irreducible,
    isAperiodic,
    isErgodic,
    isIrreducible,
    period,
    recurrentState,
    recurrentStates,
    recurrentStatesOf,
    supportEdge,
    transientState,
    transientStates,
    transientStatesOf,
    unIrreducible,
    witnessIrreducible,
 )
import Dtmc.TestSupport (
    genTransitionMatrix,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    mkTransitionMatrix,
    unTransitionMatrix,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
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

fromRows :: (Show e) => Either e (TransitionMatrix n) -> TransitionMatrix n
fromRows = either (error . show) id

threeCycle :: TransitionMatrix 3
threeCycle =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 0, 1, 0
                , 0, 0, 1
                , 1, 0, 0
                ]
            )

selfLoopTwo :: TransitionMatrix 2
selfLoopTwo =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 0.5, 0.5
                , 1.0, 0.0
                ]
            )

bipartiteTwo :: TransitionMatrix 2
bipartiteTwo =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 0, 1
                , 1, 0
                ]
            )

sevenState :: TransitionMatrix 7
sevenState =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 0, 1, 0, 0, 0, 0, 0
                , 1, 0, 0, 0, 0, 0, 0
                , 0, 0.4, 0, 0.6, 0, 0, 0
                , 0, 0, 0.3, 0, 0.7, 0, 0
                , 0, 0, 0.3, 0.4, 0, 0.3, 0
                , 0, 0, 0, 0, 0.2, 0, 0.8
                , 0, 0, 0, 0, 0, 0, 1
                ]
            )

identityThree :: TransitionMatrix 3
identityThree =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 1, 0, 0
                , 0, 1, 0
                , 0, 0, 1
                ]
            )

-- Exercise 3.2.2: irreducible, period 2, cyclic classes {A,B} and {C,D}.
fourStateCyclic :: TransitionMatrix 4
fourStateCyclic =
    fromRows $
        mkTransitionMatrix
            ( S.matrix
                [ 0, 0, 1, 0
                , 0, 0, 0, 1
                , 0.5, 0.5, 0, 0
                , 1, 0, 0, 0
                ]
            )

matrixSupport :: (KnownNat n) => TransitionMatrix n -> [[Bool]]
matrixSupport =
    map (map (> 0)) . LA.toLists . S.extract . unTransitionMatrix

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

classesAsInts :: (KnownNat n) => TransitionMatrix n -> [[Integer]]
classesAsInts = map (map getFinite) . communicatingClasses

cyclicClassesAsInts :: (KnownNat n) => TransitionMatrix n -> Maybe [[Integer]]
cyclicClassesAsInts = fmap (map (map getFinite)) . cyclicClasses

sortUnique :: (Ord a) => [a] -> [a]
sortUnique = foldr insert []
  where
    insert x [] = [x]
    insert x (y : ys)
        | x < y = x : y : ys
        | x == y = y : ys
        | otherwise = y : insert x ys

periodMatchesReference :: (KnownNat n) => TransitionMatrix n -> [Finite n] -> Property
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
            forAll (genTransitionMatrix @4) $ \matrix ->
                case mkTransitionMatrix matrix of
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
            forAll (genTransitionMatrix @4) $ \matrix ->
                case mkTransitionMatrix matrix of
                    Right p ->
                        conjoin
                            [ property (accessible p i i)
                            | i <- finites :: [Finite 4]
                            ]
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

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
            forAll (genTransitionMatrix @4) $ \matrix ->
                case mkTransitionMatrix matrix of
                    Right p -> periodMatchesReference p (finites :: [Finite 4])
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

        prop "agrees with the gcd of return-time lengths (random @3)" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                case mkTransitionMatrix matrix of
                    Right p -> periodMatchesReference p (finites :: [Finite 3])
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

    describe "communicatingClasses" $ do
        it "splits the seven-state chain into {A,B}, {C,D,E,F}, {G}" $
            classesAsInts sevenState `shouldBe` [[0, 1], [2, 3, 4, 5], [6]]

        it "returns a single class for the irreducible three-cycle" $
            classesAsInts threeCycle `shouldBe` [[0, 1, 2]]

        prop "the classes partition the state space (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                case mkTransitionMatrix matrix of
                    Right p ->
                        property (sortUnique (concat (classesAsInts p)) == [0 .. 3])
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

        prop "communication agrees with the class partition (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                case mkTransitionMatrix matrix of
                    Right p ->
                        let states = finites :: [Finite 4]
                            classIx = communicatingClasses p
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
            forAll (genTransitionMatrix @4) $ \matrix ->
                case mkTransitionMatrix matrix of
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

    describe "classify" $ do
        it "records members, periods, and closedness for the seven-state chain" $ do
            let cs = classesOf (classify sevenState)
            map (map getFinite . classMembers) cs
                `shouldBe` [[0, 1], [2, 3, 4, 5], [6]]
            map classPeriod cs `shouldBe` [Just 2, Just 1, Just 1]
            map classClosed cs `shouldBe` [True, False, True]

    describe "classify report" $ do
        it "finds the absorbing states" $ do
            map getFinite (absorbingStates (classify sevenState)) `shouldBe` [6]
            map getFinite (absorbingStates (classify identityThree)) `shouldBe` [0, 1, 2]
            map getFinite (absorbingStates (classify threeCycle)) `shouldBe` []

        prop "absorbing states have only a self-loop (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                case mkTransitionMatrix matrix of
                    Right p ->
                        conjoin
                            [ counterexample (show i) $
                                [j | j <- finites :: [Finite 4], supportEdge p i j] === [i]
                            | i <- absorbingStates (classify p)
                            ]
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

        prop "report fields agree with the standalone queries (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                case mkTransitionMatrix matrix of
                    Right p ->
                        let c = classify p
                            firstState = minBound :: Finite 4
                         in conjoin
                                [ counterexample "isIrreducible" $
                                    isIrreducible c === irreducible p
                                , counterexample "isAperiodic" $
                                    isAperiodic c === aperiodic p
                                , counterexample "isErgodic" $
                                    isErgodic c === (irreducible p && aperiodic p)
                                , counterexample "recurrentStatesOf" $
                                    recurrentStatesOf c === recurrentStates p
                                , counterexample "transientStatesOf" $
                                    transientStatesOf c === transientStates p
                                , counterexample "chainPeriod" $
                                    chainPeriod c
                                        === (if irreducible p then period p firstState else Nothing)
                                ]
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
            forAll (genTransitionMatrix @4) $ \matrix ->
                case mkTransitionMatrix matrix of
                    Right p ->
                        sort
                            ( map getFinite (recurrentStates p)
                                <> map getFinite (transientStates p)
                            )
                            === [0 .. 3]
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

        prop "every finite chain has a recurrent state (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                case mkTransitionMatrix matrix of
                    Right p ->
                        property (not (null (recurrentStates p)))
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

        prop "transient iff some reachable state cannot reach back (random @4)" $
            forAll (genTransitionMatrix @4) $ \matrix ->
                case mkTransitionMatrix matrix of
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
            forAll (genTransitionMatrix @4) $ \matrix ->
                case mkTransitionMatrix matrix of
                    Right p ->
                        let states = finites :: [Finite 4]
                         in conjoin
                                [ recurrentState p i === (i `elem` recurrentStates p)
                                | i <- states
                                ]
                    Left err ->
                        counterexample ("generated matrix was rejected: " <> show err) False

    describe "witnessIrreducible" $ do
        it "produces a witness exactly for irreducible chains" $ do
            fmap (sameMatrix threeCycle . unIrreducible) (witnessIrreducible threeCycle)
                `shouldBe` Just True
            (witnessIrreducible sevenState >> Just ()) `shouldBe` Nothing

sameMatrix :: (KnownNat n) => TransitionMatrix n -> TransitionMatrix n -> Bool
sameMatrix a b =
    flat a == flat b
  where
    flat = LA.toList . LA.flatten . S.extract . unTransitionMatrix
