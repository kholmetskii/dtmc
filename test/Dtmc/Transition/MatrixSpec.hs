{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Transition.MatrixSpec (
    spec,
) where

import Data.Finite (
    Finite,
    finites,
    getFinite,
 )
import Dtmc.Distribution (
    distributionWeights,
    probabilityAt,
 )
import Dtmc.Distribution.Vector (
    mkDistributionVector,
    unDistributionVector,
 )
import Dtmc.Probability (
    transitionProbability,
    transitionProbabilityN,
 )
import Dtmc.Simplex (SimplexError (..))
import Dtmc.State (
    FiniteState,
 )
import Dtmc.TestSupport (
    approxEq,
    approxTransitionMatrixEq,
    bumpSmallestInFirstRow,
    genTransitionMatrix,
    modifyMatrixRows,
    setFirstEntry,
    testTolerance,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    TransitionMatrixError (..),
    identityMatrix,
    matrixPower,
    mkTransitionMatrix,
    mulTransitionMatrix,
    rowAt,
    unTransitionMatrix,
 )
import GHC.Generics (
    Generic,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )
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
    choose,
    conjoin,
    counterexample,
    forAll,
    property,
    (===),
 )

data NamedPhase = PhaseA | PhaseB | PhaseC
    deriving (Eq, Ord, Show, Generic)

instance FiniteState NamedPhase

cyclicThree :: TransitionMatrix (Finite 3)
cyclicThree =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [0, 1, 0, 0, 0, 1, 1, 0, 0]
            )

namedCycle :: TransitionMatrix NamedPhase
namedCycle =
    either (error . show) id $
        mkTransitionMatrix @NamedPhase
            (S.matrix [0, 1, 0, 0, 0, 1, 1, 0, 0] :: S.Sq 3)

twoState :: TransitionMatrix (Finite 2)
twoState =
    either (error . show) id $
        mkTransitionMatrix
            (S.matrix [0.9, 0.1, 0.4, 0.6] :: S.Sq 2)

twoStateSquared :: TransitionMatrix (Finite 2)
twoStateSquared =
    either (error . show) id $
        mkTransitionMatrix
            (S.matrix [0.85, 0.15, 0.6, 0.4] :: S.Sq 2)

-- Assignment 1: states ordered [A, B, C, D, E].
assignment1 :: TransitionMatrix (Finite 5)
assignment1 =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 0
                , 0
                , 0
                , 1
                , 0
                , 1 / 3
                , 0
                , 0
                , 0
                , 2 / 3
                , 0
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 1 / 3
                , 2 / 3
                , 0
                , 1 / 4
                , 1 / 4
                , 0
                , 0
                , 1 / 2
                ] ::
                S.Sq 5
            )

-- Assignment 2: an unnamed three-state chain.
assignment2 :: TransitionMatrix (Finite 3)
assignment2 =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 0.1
                , 0.5
                , 0.4
                , 0.1
                , 0.8
                , 0.1
                , 0.0
                , 0.5
                , 0.5
                ] ::
                S.Sq 3
            )

-- Assignment 3: the seven-state fruit chain, ordered
-- [apple, pear, banana, mango, kiwi, watermelon, grapefruit]. It has period 3
-- with cyclic classes {apple, pear}, {banana, mango},
-- {kiwi, watermelon, grapefruit}.
assignment3 :: TransitionMatrix (Finite 7)
assignment3 =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 0
                , 0
                , 1 / 2
                , 1 / 2
                , 0
                , 0
                , 0 -- apple
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0 -- pear
                , 0
                , 0
                , 0
                , 0
                , 1 / 3
                , 1 / 3
                , 1 / 3 -- banana
                , 0
                , 0
                , 0
                , 0
                , 0
                , 2 / 3
                , 1 / 3 -- mango
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0 -- kiwi
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0 -- watermelon
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0 -- grapefruit
                ] ::
                S.Sq 7
            )

negSixth :: Double
negSixth = -(1 / 6)

-- | Assignment 2 closed form for @P^n(2, 0)@.
formula2 :: Int -> Double
formula2 n =
    5 / 63 + 5 / 18 * (0.1 ^ n) - 5 / 14 * (0.3 ^ n)

-- | Assignment 3 closed form for @P^(3n+1)(apple, mango)@.
formulaA :: Int -> Double
formulaA n =
    5 / 7 - 3 / 14 * (negSixth ^ n)

-- | Assignment 3 closed form for @P^(3n+2)(mango, pear)@.
formulaB :: Int -> Double
formulaB n =
    3 / 7 - 2 / 21 * (negSixth ^ n)

spec :: Spec
spec = do
    describe "mkTransitionMatrix" $ do
        prop "preserves the validated matrix" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                case mkTransitionMatrix @(Finite 3) matrix of
                    Right transitionMatrix ->
                        S.extract (unTransitionMatrix transitionMatrix)
                            === S.extract matrix
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

        prop "identifies a row whose sum is invalid" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                let invalid =
                        modifyMatrixRows
                            (bumpSmallestInFirstRow 1e-6)
                            matrix
                 in case mkTransitionMatrix @(Finite 3) invalid of
                        Left (InRow 0 (SumOffBy _)) ->
                            property True
                        result ->
                            counterexample
                                ("expected InRow 0 SumOffBy, got " <> show result)
                                False

        prop "identifies a negative entry by row and column" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                let invalid =
                        modifyMatrixRows
                            (setFirstEntry (-1e-6))
                            matrix
                 in case mkTransitionMatrix @(Finite 3) invalid of
                        Left (InRow 0 (NegativeEntry 0 _)) ->
                            property True
                        result ->
                            counterexample
                                ("expected InRow 0 NegativeEntry 0, got " <> show result)
                                False

    describe "mulTransitionMatrix" $ do
        prop "is closed under multiplication"
            $ forAll
                ((,) <$> genTransitionMatrix @3 <*> genTransitionMatrix @3)
            $ \(left, right) ->
                case (mkTransitionMatrix @(Finite 3) left, mkTransitionMatrix @(Finite 3) right) of
                    (Right leftMatrix, Right rightMatrix) ->
                        case mkTransitionMatrix @(Finite 3)
                            ( unTransitionMatrix
                                (mulTransitionMatrix leftMatrix rightMatrix)
                            ) of
                            Right _ ->
                                property True
                            Left err ->
                                counterexample
                                    ("matrix product was rejected: " <> show err)
                                    False
                    result ->
                        counterexample
                            ("generated matrix was rejected: " <> show result)
                            False

        prop "approximately equals itself at zero tolerance" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                case mkTransitionMatrix @(Finite 3) matrix of
                    Right transitionMatrix ->
                        property
                            ( approxTransitionMatrixEq
                                0
                                transitionMatrix
                                transitionMatrix
                            )
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

    describe "TransitionMatrix Semigroup" $ do
        prop "composition is approximately associative"
            $ forAll
                ( (,,)
                    <$> genTransitionMatrix @3
                    <*> genTransitionMatrix @3
                    <*> genTransitionMatrix @3
                )
            $ \(matrixA, matrixB, matrixC) ->
                case ( mkTransitionMatrix @(Finite 3) matrixA
                     , mkTransitionMatrix @(Finite 3) matrixB
                     , mkTransitionMatrix @(Finite 3) matrixC
                     ) of
                    (Right a, Right b, Right c) ->
                        property $
                            approxTransitionMatrixEq
                                1e-9
                                ((a <> b) <> c)
                                (a <> (b <> c))
                    result ->
                        counterexample
                            ("generated matrices were rejected: " <> show result)
                            False

    describe "TransitionMatrix Monoid" $ do
        prop "has a left identity" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                case mkTransitionMatrix @(Finite 3) matrix of
                    Right p ->
                        property $
                            approxTransitionMatrixEq
                                1e-12
                                (mempty <> p)
                                p
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

        prop "has a right identity" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                case mkTransitionMatrix @(Finite 3) matrix of
                    Right p ->
                        property $
                            approxTransitionMatrixEq
                                1e-12
                                (p <> mempty)
                                p
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

        it "uses the identity transition matrix as mempty" $
            approxTransitionMatrixEq
                1e-12
                (mempty :: TransitionMatrix (Finite 2))
                identityMatrix
                `shouldBe` True

    describe "matrixPower" $ do
        it "returns the identity at exponent zero" $
            approxTransitionMatrixEq
                1e-12
                (matrixPower 0 twoState)
                identityMatrix
                `shouldBe` True

        it "returns the matrix itself at exponent one" $
            approxTransitionMatrixEq
                1e-12
                (matrixPower 1 twoState)
                twoState
                `shouldBe` True

        it "matches a hand-computed square" $
            approxTransitionMatrixEq
                1e-9
                (matrixPower 2 twoState)
                twoStateSquared
                `shouldBe` True

        prop "stays stochastic for small exponents" $
            forAll ((,) <$> choose (0, 6 :: Int) <*> genTransitionMatrix @3) $
                \(k, matrix) ->
                    case mkTransitionMatrix @(Finite 3) matrix of
                        Right p ->
                            case mkTransitionMatrix @(Finite 3)
                                ( unTransitionMatrix
                                    (matrixPower (fromIntegral k) p)
                                ) of
                                Right _ ->
                                    property True
                                Left err ->
                                    counterexample
                                        ("power left the stochastic set: " <> show err)
                                        False
                        Left err ->
                            counterexample
                                ("generated matrix was rejected: " <> show err)
                                False

        prop "satisfies the power addition law"
            $ forAll
                ( (,,)
                    <$> choose (0, 6 :: Int)
                    <*> choose (0, 6 :: Int)
                    <*> genTransitionMatrix @3
                )
            $ \(m, n, matrix) ->
                case mkTransitionMatrix @(Finite 3) matrix of
                    Right p ->
                        property $
                            approxTransitionMatrixEq
                                1e-9
                                (matrixPower (fromIntegral (m + n)) p)
                                ( matrixPower (fromIntegral m) p
                                    <> matrixPower (fromIntegral n) p
                                )
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

    describe "rowAt" $ do
        it "reads rows rather than columns" $
            LA.toList
                (S.extract (unDistributionVector (rowAt cyclicThree 0)))
                `shouldBe` [0, 1, 0]

        it "returns each row of the three-cycle" $ do
            let row index =
                    LA.toList
                        (S.extract (unDistributionVector (rowAt cyclicThree index)))

            row 0 `shouldBe` [0, 1, 0]
            row 1 `shouldBe` [0, 0, 1]
            row 2 `shouldBe` [1, 0, 0]

        prop "always returns a valid distribution" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                case mkTransitionMatrix @(Finite 3) matrix of
                    Right transitionMatrix ->
                        conjoin
                            [ case mkDistributionVector @(Finite 3)
                                (unDistributionVector (rowAt transitionMatrix index)) of
                                Right _ ->
                                    property True
                                Left err ->
                                    counterexample
                                        ("row was rejected: " <> show err)
                                        False
                            | index <- [0 .. 2]
                            ]
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

    describe "transitionProbability" $ do
        prop "agrees with rowAt then probabilityAt" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                case mkTransitionMatrix @(Finite 3) matrix of
                    Right p ->
                        conjoin
                            [ transitionProbability p i j
                                === probabilityAt (rowAt p i) j
                            | i <- finites
                            , j <- finites
                            ]
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

    describe "transitionProbabilityN" $ do
        it "is the Kronecker delta at exponent zero" $
            let ijs =
                    [(i, j) | i <- finites, j <- finites] ::
                        [(Finite 2, Finite 2)]
             in map (uncurry (transitionProbabilityN 0 twoState)) ijs
                    `shouldBe` map (\(i, j) -> if i == j then 1 else 0) ijs

        prop "agrees with transitionProbability at exponent one" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                case mkTransitionMatrix @(Finite 3) matrix of
                    Right p ->
                        conjoin
                            [ property $
                                approxEq
                                    testTolerance
                                    (transitionProbabilityN 1 p i j)
                                    (transitionProbability p i j)
                            | i <- finites
                            , j <- finites
                            ]
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

        it "matches a hand-computed square at exponent two" $ do
            approxEq
                testTolerance
                (transitionProbabilityN 2 twoState 0 0)
                (probabilityAt (rowAt twoStateSquared 0) 0)
                `shouldBe` True
            approxEq
                testTolerance
                (transitionProbabilityN 2 twoState 0 1)
                (probabilityAt (rowAt twoStateSquared 0) 1)
                `shouldBe` True
            approxEq
                testTolerance
                (transitionProbabilityN 2 twoState 1 0)
                (probabilityAt (rowAt twoStateSquared 1) 0)
                `shouldBe` True
            approxEq
                testTolerance
                (transitionProbabilityN 2 twoState 1 1)
                (probabilityAt (rowAt twoStateSquared 1) 1)
                `shouldBe` True

        prop "agrees with the corresponding matrixPower entry" $
            forAll ((,) <$> choose (0, 6 :: Int) <*> genTransitionMatrix @3) $
                \(k, matrix) ->
                    case mkTransitionMatrix @(Finite 3) matrix of
                        Right p ->
                            let power = matrixPower (fromIntegral k) p
                                raw = S.extract (unTransitionMatrix power)
                             in conjoin
                                    [ property $
                                        approxEq
                                            testTolerance
                                            ( transitionProbabilityN
                                                (fromIntegral k)
                                                p
                                                i
                                                j
                                            )
                                            ( raw
                                                `LA.atIndex` ( fromIntegral (getFinite i)
                                                             , fromIntegral (getFinite j)
                                                             )
                                            )
                                    | i <- finites
                                    , j <- finites
                                    ]
                        Left err ->
                            counterexample
                                ("generated matrix was rejected: " <> show err)
                                False

    describe "transitionProbabilityN assignment regressions" $ do
        -- Assignment 1: P^3(E, D) = 3/8 (E = 4, D = 3).
        it "assignment 1: P^3(E, D) = 3/8" $
            approxEq
                testTolerance
                (transitionProbabilityN 3 assignment1 4 3)
                (3 / 8)
                `shouldBe` True

        -- Assignment 2: P^n(2, 0) = 5/63 + 5/18 (1/10)^n - 5/14 (3/10)^n.
        it "assignment 2: P^n(2, 0) closed form" $
            mapM_
                ( \n ->
                    approxEq
                        testTolerance
                        (transitionProbabilityN n assignment2 2 0)
                        (formula2 (fromIntegral n))
                        `shouldBe` True
                )
                ([0, 1, 2, 3, 5, 10, 20] :: [Natural])

        -- Assignment 3: P^(3n+1)(apple, mango) = 5/7 - 3/14 (-1/6)^n, apple = 0,
        -- mango = 3. Includes the large index n = 675 (matrix power 2026).
        it "assignment 3: apple -> mango closed form" $
            mapM_
                ( \n ->
                    approxEq
                        testTolerance
                        (transitionProbabilityN (3 * n + 1) assignment3 0 3)
                        (formulaA (fromIntegral n))
                        `shouldBe` True
                )
                ([0, 1, 2, 3, 675] :: [Natural])

        -- Assignment 3: P^(3n+2)(mango, pear) = 3/7 - 2/21 (-1/6)^n, mango = 3,
        -- pear = 1.
        it "assignment 3: mango -> pear closed form" $
            mapM_
                ( \n ->
                    approxEq
                        testTolerance
                        (transitionProbabilityN (3 * n + 2) assignment3 3 1)
                        (formulaB (fromIntegral n))
                        `shouldBe` True
                )
                ([0, 1, 2, 3, 4] :: [Natural])

    describe "named finite states" $ do
        it "returns a row labelled by named constructors" $
            distributionWeights (rowAt namedCycle PhaseA)
                `shouldBe` [(PhaseB, 1)]

        it "uses named constructors for one-step probabilities" $
            transitionProbability namedCycle PhaseB PhaseC
                `shouldBe` 1

        it "preserves the named state type through powers" $
            transitionProbability (matrixPower 2 namedCycle) PhaseA PhaseC
                `shouldBe` 1

        it "provides a named identity matrix" $
            transitionProbability (identityMatrix @NamedPhase) PhaseB PhaseB
                `shouldBe` 1

        it "composes matrices without changing their named state type" $
            approxTransitionMatrixEq
                0
                (namedCycle <> identityMatrix)
                namedCycle
                `shouldBe` True
