module Dtmc.DistributionSpec (
    spec,
) where

import Data.Finite (
    Finite,
 )
import Data.Map.Strict qualified as Map
import Dtmc.Distribution (
    Distribution (..),
    DistributionError (..),
 )
import Dtmc.Distribution.Map (
    DistributionMap,
    mkDistributionMap,
    pointMass,
    toDistributionMap,
    unDistributionMap,
 )
import Dtmc.Distribution.Vector (
    mkDistributionVector,
    unDistributionVector,
 )
import Dtmc.Simplex (
    SimplexError (..),
 )
import Dtmc.TestSupport (
    approxDistributionEq,
    approxEq,
    bumpSmallest,
    genSimplexPoint,
    testTolerance,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    expectationFailure,
    it,
    shouldBe,
    shouldSatisfy,
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

spec :: Spec
spec = do
    describe "mkDistributionVector" $ do
        it "rejects an empty vector" $
            case mkDistributionVector (S.vector [] :: S.R 0) of
                Left err ->
                    err `shouldBe` DistributionError (SumOffBy 0)
                Right _ ->
                    expectationFailure "expected rejection"

        it "accepts a tiny negative rounding error" $
            case mkDistributionVector (S.vector [-1e-17, 1] :: S.R 2) of
                Right _ ->
                    pure ()
                Left err ->
                    expectationFailure
                        ("expected acceptance, got " <> show err)

        it "reports an entry above one" $
            case mkDistributionVector (S.vector [1.5, -0.5] :: S.R 2) of
                Left err ->
                    err
                        `shouldBe` DistributionError (EntryAboveOne 0 1.5)
                Right _ ->
                    expectationFailure "expected rejection"

        prop "accepts normalised vectors" $
            forAll (genSimplexPoint 3) $ \entries ->
                case mkDistributionVector (S.vector entries :: S.R 3) of
                    Right _ ->
                        property True
                    Left err ->
                        counterexample
                            ("generated vector was rejected: " <> show err)
                            False

        prop "rejects vectors whose sum is too large" $
            forAll (genSimplexPoint 3) $ \entries ->
                case mkDistributionVector
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
                 in case mkDistributionVector (S.vector invalid :: S.R 3) of
                        Left (DistributionError (NegativeEntry 0 _)) ->
                            property True
                        result ->
                            counterexample
                                ("expected DistributionError NegativeEntry 0, got " <> show result)
                                False

        prop "preserves the validated vector" $
            forAll (genSimplexPoint 3) $ \entries ->
                let simplexVector = S.vector entries :: S.R 3
                 in case mkDistributionVector simplexVector of
                        Right distribution ->
                            S.extract (unDistributionVector distribution)
                                === S.extract simplexVector
                        Left err ->
                            counterexample
                                ("generated vector was rejected: " <> show err)
                                False

    describe "DistributionMap" $ do
        it "combines duplicates and stores canonical ascending entries" $ do
            let distribution =
                    either (error . show) id $
                        mkDistributionMap
                            [('b', 0.2), ('a', 0.5), ('b', 0.3), ('c', 0)]
            distributionWeights distribution `shouldBe` [('a', 0.5), ('b', 0.5)]
            support distribution `shouldBe` ['a', 'b']
            Map.toAscList (unDistributionMap distribution)
                `shouldBe` [('a', 0.5), ('b', 0.5)]

        it "returns zero for an absent state" $
            probabilityAt (pointMass "present") "absent"
                `shouldBe` 0

        it "rejects an empty law" $
            (mkDistributionMap [] :: Either DistributionError (DistributionMap Int))
                `shouldSatisfy` either (const True) (const False)

        it "uses the shared error type" $
            mkDistributionMap ([] :: [(Int, Double)])
                `shouldBe` Left (DistributionError (SumOffBy 0))

    describe "Distribution abstraction" $ do
        let vector =
                either (error . show) id $
                    mkDistributionVector (S.vector [0.2, 0, 0.8] :: S.R 3)
            mapDistribution =
                either (error . show) id $
                    ( mkDistributionMap [(0, 0.2), (2, 0.8)] ::
                        Either DistributionError (DistributionMap (Finite 3))
                    )

        it "exposes the same weights and support for both representations" $ do
            distributionWeights vector `shouldBe` distributionWeights mapDistribution
            support vector `shouldBe` support mapDistribution

        it "converts both representations to the same canonical map" $ do
            toDistributionMap vector `shouldBe` mapDistribution
            toDistributionMap mapDistribution `shouldBe` mapDistribution

    describe "probabilityAt" $ do
        let known =
                either (error . show) id $
                    mkDistributionVector (S.vector [0.2, 0.5, 0.3] :: S.R 3)

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
                        mkDistributionVector (S.vector [-1e-17, 1] :: S.R 2)

            probabilityAt tolerated 0 `shouldBe` (-1e-17)
            probabilityAt tolerated 1 `shouldBe` 1

    describe "approxDistributionEq" $
        prop "is reflexive at zero tolerance" $
            forAll (genSimplexPoint 3) $ \entries ->
                case mkDistributionVector (S.vector entries :: S.R 3) of
                    Right distribution ->
                        property
                            (approxDistributionEq 0 distribution distribution)
                    Left err ->
                        counterexample
                            ("generated vector was rejected: " <> show err)
                            False
