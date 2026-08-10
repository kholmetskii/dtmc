{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Distribution.VectorSpec (
    spec,
) where

import Data.Finite (
    Finite,
 )
import Dtmc.Distribution (
    Distribution (..),
    DistributionError (..),
 )
import Dtmc.Distribution.Vector (
    mkDistributionVector,
    unDistributionVector,
 )
import Dtmc.Simplex (
    SimplexError (..),
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.TestSupport (
    approxEq,
    bumpSmallest,
    genSimplexPoint,
    testTolerance,
 )
import GHC.Generics (
    Generic,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    expectationFailure,
    it,
    shouldBe,
 )
import Test.Hspec.QuickCheck (
    prop,
 )
import Test.QuickCheck (
    counterexample,
    forAll,
    property,
    (===),
 )

data NamedState = NamedA | NamedB | NamedC
    deriving (Eq, Ord, Show, Generic)

instance FiniteState NamedState

spec :: Spec
spec = do
    describe "mkDistributionVector" $ do
        it "rejects an empty vector" $
            case mkDistributionVector @(Finite 0) (S.vector [] :: S.R 0) of
                Left err ->
                    err `shouldBe` DistributionError (SumOffBy 0)
                Right _ ->
                    expectationFailure "expected rejection"

        it "accepts a tiny negative rounding error" $
            case mkDistributionVector @(Finite 2) (S.vector [-1e-17, 1] :: S.R 2) of
                Right _ ->
                    pure ()
                Left err ->
                    expectationFailure
                        ("expected acceptance, got " <> show err)

        it "reports an entry above one" $
            case mkDistributionVector @(Finite 2) (S.vector [1.5, -0.5] :: S.R 2) of
                Left err ->
                    err
                        `shouldBe` DistributionError (EntryAboveOne 0 1.5)
                Right _ ->
                    expectationFailure "expected rejection"

        prop "accepts normalised vectors" $
            forAll (genSimplexPoint 3) $ \entries ->
                case mkDistributionVector @(Finite 3) (S.vector entries :: S.R 3) of
                    Right _ ->
                        property True
                    Left err ->
                        counterexample
                            ("generated vector was rejected: " <> show err)
                            False

        prop "rejects vectors whose sum is too large" $
            forAll (genSimplexPoint 3) $ \entries ->
                case mkDistributionVector @(Finite 3)
                    (S.vector (bumpSmallest 1e-6 entries) :: S.R 3) of
                    Left (DistributionError (SumOffBy _)) ->
                        property True
                    result ->
                        counterexample
                            ("expected DistributionError SumOffBy, got " <> show result)
                            False

        prop "rejects genuinely negative entries" $
            forAll (genSimplexPoint 3) $ \entries ->
                let invalid =
                        case entries of
                            _ : rest -> (-1e-6) : rest
                            [] -> []
                 in case mkDistributionVector @(Finite 3) (S.vector invalid :: S.R 3) of
                        Left (DistributionError (NegativeEntry 0 _)) ->
                            property True
                        result ->
                            counterexample
                                ("expected DistributionError NegativeEntry 0, got " <> show result)
                                False

        prop "preserves the validated vector" $
            forAll (genSimplexPoint 3) $ \entries ->
                let simplexVector = S.vector entries :: S.R 3
                 in case mkDistributionVector @(Finite 3) simplexVector of
                        Right distribution ->
                            S.extract (unDistributionVector distribution)
                                === S.extract simplexVector
                        Left err ->
                            counterexample
                                ("generated vector was rejected: " <> show err)
                                False

    describe "probabilityAt" $ do
        let known =
                either (error . show) id $
                    mkDistributionVector @(Finite 3) (S.vector [0.2, 0.5, 0.3] :: S.R 3)

        it "returns each coordinate of a known distribution" $ do
            approxEq testTolerance (probabilityAt known 0) 0.2 `shouldBe` True
            approxEq testTolerance (probabilityAt known 1) 0.5 `shouldBe` True
            approxEq testTolerance (probabilityAt known 2) 0.3 `shouldBe` True

        it "reads the first and last valid states" $ do
            approxEq testTolerance (probabilityAt known minBound) 0.2
                `shouldBe` True
            approxEq testTolerance (probabilityAt known maxBound) 0.3
                `shouldBe` True

        it "returns tolerated stored values without clamping" $ do
            let tolerated =
                    either (error . show) id $
                        mkDistributionVector @(Finite 2) (S.vector [-1e-17, 1] :: S.R 2)

            probabilityAt tolerated 0 `shouldBe` (-1e-17)
            probabilityAt tolerated 1 `shouldBe` 1

    describe "named finite states" $ do
        let namedDistribution =
                either (error . show) id $
                    mkDistributionVector @NamedState
                        (S.vector [0.2, 0, 0.8] :: S.R 3)
            indexedDistribution =
                either (error . show) id $
                    mkDistributionVector @(Finite 3)
                        (S.vector [0.2, 0, 0.8] :: S.R 3)

        it "indexes coordinates by state constructors" $ do
            probabilityAt namedDistribution NamedA `shouldBe` 0.2
            probabilityAt namedDistribution NamedB `shouldBe` 0
            probabilityAt namedDistribution NamedC `shouldBe` 0.8

        it "reports weights and support in constructor order" $ do
            distributionWeights namedDistribution
                `shouldBe` [(NamedA, 0.2), (NamedC, 0.8)]
            support namedDistribution `shouldBe` [NamedA, NamedC]

        it "matches the low-level indexed representation coordinate for coordinate" $
            map (probabilityAt namedDistribution) [NamedA, NamedB, NamedC]
                `shouldBe` map (probabilityAt indexedDistribution) [0, 1, 2]
