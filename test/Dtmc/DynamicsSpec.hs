{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.DynamicsSpec (
    spec,
) where

import Data.Finite (
    finites,
 )
import Data.Proxy (
    Proxy (..),
 )
import Dtmc.Distribution (
    Distribution,
    mkDistribution,
    probabilityAt,
    unDistribution,
 )
import Dtmc.Dynamics (
    evolve,
    evolveN,
 )
import Dtmc.Probability (
    probabilityAtTime,
 )
import Dtmc.TestSupport (
    approxDistributionEq,
    approxEq,
    genSimplexPoint,
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    mkTransitionMatrix,
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
    conjoin,
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

-- Assignment 1: states ordered [A, B, C, D, E].
assignment1 :: TransitionMatrix 5
assignment1 =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 0, 0, 0, 1, 0
                , 1 / 3, 0, 0, 0, 2 / 3
                , 0, 0, 0, 0, 1
                , 0, 0, 1 / 3, 2 / 3, 0
                , 1 / 4, 1 / 4, 0, 0, 1 / 2
                ] ::
                S.Sq 5
            )

-- Assignment 1 initial law lambda = [1/4, 1/2, 0, 1/4, 0].
assignment1Initial :: Distribution 5
assignment1Initial =
    either (error . show) id $
        mkDistribution (S.vector [1 / 4, 1 / 2, 0, 1 / 4, 0] :: S.R 5)

spec :: Spec
spec = do
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

    describe "probabilityAtTime" $ do
        it "returns the initial probability at time zero" $ do
            let mu =
                    either (error . show) id $
                        mkDistribution (S.vector [0.25, 0.75] :: S.R 2)

            conjoin
                [ property $
                    approxEq
                        testTolerance
                        (probabilityAtTime 0 mu twoState j)
                        (probabilityAt mu j)
                | j <- finites
                ]

        prop "agrees with probabilityAt of evolveN"
            $ forAll
                ( (,,)
                    <$> choose (0, 6 :: Int)
                    <*> genDistribution @3
                    <*> genTransitionMatrix @3
                )
            $ \(k, vector, matrix) ->
                case (mkDistribution vector, mkTransitionMatrix matrix) of
                    (Right mu, Right p) ->
                        conjoin
                            [ property $
                                approxEq
                                    testTolerance
                                    (probabilityAtTime (fromIntegral k) mu p j)
                                    (probabilityAt (evolveN (fromIntegral k) mu p) j)
                            | j <- finites
                            ]
                    result ->
                        counterexample
                            ("generated input was rejected: " <> show result)
                            False

        prop "agrees with repeated evolve for small exponents"
            $ forAll
                ( (,,)
                    <$> choose (0, 6 :: Int)
                    <*> genDistribution @3
                    <*> genTransitionMatrix @3
                )
            $ \(k, vector, matrix) ->
                case (mkDistribution vector, mkTransitionMatrix matrix) of
                    (Right mu, Right p) ->
                        let iterated = iterate (`evolve` p) mu !! k
                         in conjoin
                                [ property $
                                    approxEq
                                        testTolerance
                                        (probabilityAtTime (fromIntegral k) mu p j)
                                        (probabilityAt iterated j)
                                | j <- finites
                                ]
                    result ->
                        counterexample
                            ("generated input was rejected: " <> show result)
                            False

    describe "probabilityAtTime assignment regressions" $ do
        -- Assignment 1: P(X_2 = C) = 5/36 (C = 2).
        it "assignment 1: P(X_2 = C) = 5/36" $
            approxEq
                testTolerance
                (probabilityAtTime 2 assignment1Initial assignment1 2)
                (5 / 36)
                `shouldBe` True

        -- Assignment 1: P(X_3 = E) = 23/72 (E = 4).
        it "assignment 1: P(X_3 = E) = 23/72" $
            approxEq
                testTolerance
                (probabilityAtTime 3 assignment1Initial assignment1 4)
                (23 / 72)
                `shouldBe` True
