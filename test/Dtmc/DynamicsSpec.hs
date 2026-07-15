{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.DynamicsSpec (
    spec,
) where

import Data.Proxy (
    Proxy (..),
 )
import Dtmc.Distribution (
    mkDistribution,
    unDistribution,
 )
import Dtmc.Dynamics (
    evolve,
    evolveN,
    identityMatrix,
    matrixPower,
 )
import Dtmc.TestSupport (
    approxDistributionEq,
    approxTransitionMatrixEq,
    genSimplexPoint,
    genTransitionMatrix,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    mkTransitionMatrix,
    unTransitionMatrix,
 )
import GHC.TypeNats (
    KnownNat,
    natVal,
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
    Gen,
    choose,
    counterexample,
    forAll,
    property,
 )

genDistribution :: forall n. (KnownNat n) => Gen (S.R n)
genDistribution = do
    entries <- genSimplexPoint (fromIntegral (natVal (Proxy @n)))
    pure (S.vector entries)

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

    describe "evolve" $ do
        prop "keeps the distribution on the simplex" $
            forAll ((,) <$> genDistribution @3 <*> genTransitionMatrix @3) $
                \(vector, matrix) ->
                    case (mkDistribution vector, mkTransitionMatrix matrix) of
                        (Right mu, Right p) ->
                            case mkDistribution (unDistribution (evolve mu p)) of
                                Right _ ->
                                    property True
                                Left err ->
                                    counterexample
                                        ("evolved distribution left the simplex: " <> show err)
                                        False
                        result ->
                            counterexample
                                ("generated input was rejected: " <> show result)
                                False

        it "matches a hand-computed two-state step" $ do
            let mu =
                    either (error . show) id $
                        mkDistribution
                            (S.vector [1, 0] :: S.R 2)

            LA.toList (S.extract (unDistribution (evolve mu twoState)))
                `shouldBe` [0.9, 0.1]

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

    describe "evolveN" $ do
        it "leaves a distribution unchanged after zero steps" $ do
            let mu =
                    either (error . show) id $
                        mkDistribution
                            (S.vector [0.25, 0.75] :: S.R 2)

            approxDistributionEq
                1e-12
                (evolveN 0 mu twoState)
                mu
                `shouldBe` True

        prop "agrees with iterating evolve"
            $ forAll
                ( (,,)
                    <$> choose (0, 6 :: Int)
                    <*> genDistribution @3
                    <*> genTransitionMatrix @3
                )
            $ \(k, vector, matrix) ->
                case (mkDistribution vector, mkTransitionMatrix matrix) of
                    (Right mu, Right p) ->
                        let iterated =
                                iterate (`evolve` p) mu !! k
                         in property $
                                approxDistributionEq
                                    1e-9
                                    (evolveN (fromIntegral k) mu p)
                                    iterated
                    result ->
                        counterexample
                            ("generated input was rejected: " <> show result)
                            False

        prop "composes m steps then n steps"
            $ forAll
                ( (,,,)
                    <$> choose (0, 4 :: Int)
                    <*> choose (0, 4 :: Int)
                    <*> genDistribution @3
                    <*> genTransitionMatrix @3
                )
            $ \(m, n, vector, matrix) ->
                case (mkDistribution vector, mkTransitionMatrix matrix) of
                    (Right mu, Right p) ->
                        property $
                            approxDistributionEq
                                1e-9
                                (evolveN (fromIntegral (m + n)) mu p)
                                ( evolveN
                                    (fromIntegral n)
                                    (evolveN (fromIntegral m) mu p)
                                    p
                                )
                    result ->
                        counterexample
                            ("generated input was rejected: " <> show result)
                            False
