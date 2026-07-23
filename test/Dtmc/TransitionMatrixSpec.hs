{-# LANGUAGE TypeApplications #-}

module Dtmc.TransitionMatrixSpec (
    spec,
) where

import Dtmc (
    SimplexError (..),
    TransitionMatrixError (..),
    TransitionMatrix,
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
