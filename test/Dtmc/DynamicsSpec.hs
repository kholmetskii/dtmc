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
 )
import Dtmc.TestSupport (
    approxDistributionEq,
    genSimplexPoint,
    genTransitionMatrix,
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
