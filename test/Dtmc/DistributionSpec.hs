module Dtmc.DistributionSpec (
    spec,
) where

import Dtmc.Distribution (
    DistributionError (..),
    SparseDistribution,
    SparseDistributionError,
    mkDistribution,
    mkSparseDistribution,
    pointMass,
    probabilityAt,
    sparseEntries,
    sparseProbabilityAt,
    sparseSupport,
    unDistribution,
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
    describe "mkDistribution" $ do
        it "rejects an empty vector" $
            case mkDistribution (S.vector [] :: S.R 0) of
                Left err ->
                    err `shouldBe` DistributionError (SumOffBy 0)
                Right _ ->
                    expectationFailure "expected rejection"

        it "accepts a tiny negative rounding error" $
            case mkDistribution (S.vector [-1e-17, 1] :: S.R 2) of
                Right _ ->
                    pure ()
                Left err ->
                    expectationFailure
                        ("expected acceptance, got " <> show err)

        it "reports an entry above one" $
            case mkDistribution (S.vector [1.5, -0.5] :: S.R 2) of
                Left err ->
                    err
                        `shouldBe` DistributionError (EntryAboveOne 0 1.5)
                Right _ ->
                    expectationFailure "expected rejection"

        prop "accepts normalised vectors" $
            forAll (genSimplexPoint 3) $ \entries ->
                case mkDistribution (S.vector entries :: S.R 3) of
                    Right _ ->
                        property True
                    Left err ->
                        counterexample
                            ("generated vector was rejected: " <> show err)
                            False

        prop "rejects vectors whose sum is too large" $
            forAll (genSimplexPoint 3) $ \entries ->
                case mkDistribution
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
                 in case mkDistribution (S.vector invalid :: S.R 3) of
                        Left (DistributionError (NegativeEntry 0 _)) ->
                            property True
                        result ->
                            counterexample
                                ("expected DistributionError NegativeEntry 0, got " <> show result)
                                False

        prop "preserves the validated vector" $
            forAll (genSimplexPoint 3) $ \entries ->
                let simplexVector = S.vector entries :: S.R 3
                 in case mkDistribution simplexVector of
                        Right distribution ->
                            S.extract (unDistribution distribution)
                                === S.extract simplexVector
                        Left err ->
                            counterexample
                                ("generated vector was rejected: " <> show err)
                                False

    describe "SparseDistribution" $ do
        it "combines duplicates and stores canonical ascending entries" $ do
            let distribution =
                    either (error . show) id $
                        mkSparseDistribution
                            [('b', 0.2), ('a', 0.5), ('b', 0.3), ('c', 0)]
            sparseEntries distribution `shouldBe` [('a', 0.5), ('b', 0.5)]
            sparseSupport distribution `shouldBe` ['a', 'b']

        it "returns zero for an absent state" $
            sparseProbabilityAt (pointMass "present") "absent"
                `shouldBe` 0

        it "rejects an empty law" $
            (mkSparseDistribution [] :: Either SparseDistributionError (SparseDistribution Int))
                `shouldSatisfy` either (const True) (const False)

    describe "probabilityAt" $ do
        let known =
                either (error . show) id $
                    mkDistribution (S.vector [0.2, 0.5, 0.3] :: S.R 3)

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
                        mkDistribution (S.vector [-1e-17, 1] :: S.R 2)

            probabilityAt tolerated 0 `shouldBe` (-1e-17)
            probabilityAt tolerated 1 `shouldBe` 1

    describe "approxDistributionEq" $
        prop "is reflexive at zero tolerance" $
            forAll (genSimplexPoint 3) $ \entries ->
                case mkDistribution (S.vector entries :: S.R 3) of
                    Right distribution ->
                        property
                            (approxDistributionEq 0 distribution distribution)
                    Left err ->
                        counterexample
                            ("generated vector was rejected: " <> show err)
                            False
