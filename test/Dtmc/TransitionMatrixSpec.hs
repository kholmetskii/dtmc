{-# LANGUAGE TypeApplications #-}

module Dtmc.TransitionMatrixSpec (
    spec,
) where

import Dtmc (
    SimplexError (..),
    TransitionMatrixError (..),
    TransitionMatrix,
    identityMatrix,
    matrixPower,
    mkDistribution,
    mkTransitionMatrix,
    mulTransitionMatrix,
    rowAt,
    unDistribution,
    unTransitionMatrix,
 )
import Dtmc.TestSupport (
    approxTransitionMatrixEq,
    bumpSmallestInFirstRow,
    genTransitionMatrix,
    modifyMatrixRows,
    setFirstEntry,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
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

cyclicThree :: TransitionMatrix 3
cyclicThree =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [0, 1, 0, 0, 0, 1, 1, 0, 0]
            )

twoState :: TransitionMatrix 2
twoState =
    either (error . show) id $
        mkTransitionMatrix
            (S.matrix [0.9, 0.1, 0.4, 0.6] :: S.Sq 2)

twoStateSquared :: TransitionMatrix 2
twoStateSquared =
    either (error . show) id $
        mkTransitionMatrix
            (S.matrix [0.85, 0.15, 0.6, 0.4] :: S.Sq 2)

spec :: Spec
spec = do
    describe "mkTransitionMatrix" $ do
        prop "preserves the validated matrix" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                case mkTransitionMatrix matrix of
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
                 in case mkTransitionMatrix invalid of
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
                 in case mkTransitionMatrix invalid of
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
                case (mkTransitionMatrix left, mkTransitionMatrix right) of
                    (Right leftMatrix, Right rightMatrix) ->
                        case mkTransitionMatrix
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
                case mkTransitionMatrix matrix of
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
                case ( mkTransitionMatrix matrixA
                     , mkTransitionMatrix matrixB
                     , mkTransitionMatrix matrixC
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
                case mkTransitionMatrix matrix of
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
                case mkTransitionMatrix matrix of
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
                (mempty :: TransitionMatrix 2)
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
                    case mkTransitionMatrix matrix of
                        Right p ->
                            case mkTransitionMatrix
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
                case mkTransitionMatrix matrix of
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
                (S.extract (unDistribution (rowAt cyclicThree 0)))
                `shouldBe` [0, 1, 0]

        it "returns each row of the three-cycle" $ do
            let row index =
                    LA.toList
                        (S.extract (unDistribution (rowAt cyclicThree index)))

            row 0 `shouldBe` [0, 1, 0]
            row 1 `shouldBe` [0, 0, 1]
            row 2 `shouldBe` [1, 0, 0]

        prop "always returns a valid distribution" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                case mkTransitionMatrix matrix of
                    Right transitionMatrix ->
                        conjoin
                            [ case mkDistribution
                                (unDistribution (rowAt transitionMatrix index)) of
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
