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
 )
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.Distribution.Vector (
    DistributionVectorError (..),
    fromList,
    toList,
 )
import Dtmc.Simplex (
    SimplexError (..),
 )
import Dtmc.State (
    FiniteState,
    finiteStates,
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
 )

data NamedState = NamedA | NamedB | NamedC
    deriving (Eq, Ord, Show, Generic)

instance FiniteState NamedState

spec :: Spec
spec = do
    describe "fromList" $ do
        it "reports too few weights against the state cardinality" $
            case fromList @NamedState [0.5, 0.5] of
                Left err ->
                    err `shouldBe` WrongLength 3 2
                Right _ ->
                    expectationFailure "expected rejection"

        it "checks the length before the simplex invariant" $
            case fromList @NamedState [0.5, 0.5, 0.5, 0.5] of
                Left err ->
                    err `shouldBe` WrongLength 3 4
                Right _ ->
                    expectationFailure "expected rejection"

        prop "rejects any length other than the state cardinality" $
            forAll (genSimplexPoint 3) $ \entries ->
                case fromList @(Finite 3) (take 2 entries) of
                    Left err ->
                        counterexample (show err) (err == WrongLength 3 2)
                    Right _ ->
                        counterexample "expected rejection" False

        it "reports a total outside tolerance" $
            case fromList @NamedState [0.8, 0, 0] of
                Left (InWeights (SumOffBy total)) ->
                    total `shouldBe` 0.8
                result ->
                    expectationFailure
                        ("expected InWeights SumOffBy, got " <> show result)

        it "rejects an empty vector" $
            case fromList @(Finite 0) [] of
                Left err ->
                    err `shouldBe` InWeights (SumOffBy 0)
                Right _ ->
                    expectationFailure "expected rejection"

        it "clamps a tiny negative rounding error" $
            case fromList @(Finite 2) [-1e-17, 1] of
                Right distribution ->
                    toList distribution
                        `shouldBe` [0, 1]
                Left err ->
                    expectationFailure
                        ("expected acceptance, got " <> show err)

        it "normalises an accepted total near one" $
            case fromList @(Finite 2) [0.5, 0.5 - 5e-10] of
                Right distribution -> do
                    let stored =
                            toList distribution
                    approxEq 1e-12 (sum stored) 1 `shouldBe` True
                    stored == [0.5, 0.5 - 5e-10] `shouldBe` False
                Left err ->
                    expectationFailure
                        ("expected acceptance, got " <> show err)

        it "reports NaN at its coordinate" $
            case fromList @(Finite 2) [0 / 0, 1] of
                Left err ->
                    err `shouldBe` InWeights (NonFiniteEntry 0)
                Right _ ->
                    expectationFailure "expected rejection"

        it "reports infinity at its coordinate" $
            case fromList @(Finite 2) [1, 1 / 0] of
                Left err ->
                    err `shouldBe` InWeights (NonFiniteEntry 1)
                Right _ ->
                    expectationFailure "expected rejection"

        it "reports an entry above one" $
            case fromList @(Finite 2) [1.5, -0.5] of
                Left err ->
                    err
                        `shouldBe` InWeights (EntryAboveOne 0 1.5)
                Right _ ->
                    expectationFailure "expected rejection"

        prop "accepts normalised vectors" $
            forAll (genSimplexPoint 3) $ \entries ->
                case fromList @(Finite 3) entries of
                    Right _ ->
                        property True
                    Left err ->
                        counterexample
                            ("generated vector was rejected: " <> show err)
                            False

        prop "rejects vectors whose sum is too large" $
            forAll (genSimplexPoint 3) $ \entries ->
                case fromList @(Finite 3) (bumpSmallest 1e-6 entries) of
                    Left (InWeights (SumOffBy _)) ->
                        property True
                    result ->
                        counterexample
                            ("expected InWeights SumOffBy, got " <> show result)
                            False

        prop "rejects genuinely negative entries" $
            forAll (genSimplexPoint 3) $ \entries ->
                let invalid =
                        case entries of
                            _ : rest -> (-1e-6) : rest
                            [] -> []
                 in case fromList @(Finite 3) invalid of
                        Left (InWeights (NegativeEntry 0 _)) ->
                            property True
                        result ->
                            counterexample
                                ("expected InWeights NegativeEntry 0, got " <> show result)
                                False

        prop "stores a canonical vector close to the accepted input" $
            forAll (genSimplexPoint 3) $ \entries ->
                case fromList @(Finite 3) entries of
                    Right distribution ->
                        let stored = toList distribution
                         in counterexample ("stored vector: " <> show stored) $
                                property
                                    ( all (\entry -> entry >= 0 && entry <= 1) stored
                                        && approxEq 1e-12 (sum stored) 1
                                        && and
                                            ( zipWith
                                                (approxEq testTolerance)
                                                stored
                                                entries
                                            )
                                    )
                    Left err ->
                        counterexample
                            ("generated vector was rejected: " <> show err)
                            False

    describe "fromList and toList are a positional pair" $ do
        prop "fromList accepts whatever toList produced (random @3)" $
            forAll (genSimplexPoint 3) $ \entries ->
                case fromList @(Finite 3) entries of
                    Right distribution ->
                        case fromList @(Finite 3) (toList distribution) of
                            Right again ->
                                counterexample (show (toList again)) $
                                    property
                                        ( and
                                            ( zipWith
                                                (approxEq testTolerance)
                                                (toList again)
                                                (toList distribution)
                                            )
                                        )
                            Left err ->
                                counterexample
                                    ("round trip was rejected: " <> show err)
                                    False
                    Left err ->
                        counterexample
                            ("generated vector was rejected: " <> show err)
                            False

    describe "labelled construction through the sparse representation" $ do
        it "combines duplicates and fills missing states with zero" $
            case DistributionMap.fromList
                [(NamedC, 0.5), (NamedA, 0.25), (NamedA, 0.25)] of
                Left err ->
                    expectationFailure
                        ("expected acceptance, got " <> show err)
                Right sparse ->
                    case fromList
                        [probabilityAt sparse state | state <- finiteStates] of
                        Right distribution -> do
                            toList distribution `shouldBe` [0.5, 0, 0.5]
                            distributionWeights distribution
                                `shouldBe` [(NamedA, 0.5), (NamedC, 0.5)]
                        Left err ->
                            expectationFailure
                                ("expected acceptance, got " <> show err)

        prop "agrees with the sparse representation coordinate for coordinate" $
            forAll (genSimplexPoint 3) $ \entries ->
                case DistributionMap.fromList
                    (zip (finiteStates @NamedState) entries) of
                    Left err ->
                        counterexample ("sparse rejected: " <> show err) False
                    Right sparse ->
                        case fromList @NamedState entries of
                            Right dense ->
                                counterexample (show (toList dense)) $
                                    property
                                        ( and
                                            [ approxEq
                                                testTolerance
                                                (probabilityAt sparse state)
                                                (probabilityAt dense state)
                                            | state <- finiteStates
                                            ]
                                        )
                            Left err ->
                                counterexample
                                    ("dense rejected: " <> show err)
                                    False

    describe "probabilityAt" $ do
        let known =
                either (error . show) id $
                    fromList @(Finite 3) [0.2, 0.5, 0.3]

        it "returns each coordinate of a known distribution" $ do
            approxEq testTolerance (probabilityAt known 0) 0.2 `shouldBe` True
            approxEq testTolerance (probabilityAt known 1) 0.5 `shouldBe` True
            approxEq testTolerance (probabilityAt known 2) 0.3 `shouldBe` True

        it "reads the first and last valid states" $ do
            approxEq testTolerance (probabilityAt known minBound) 0.2
                `shouldBe` True
            approxEq testTolerance (probabilityAt known maxBound) 0.3
                `shouldBe` True

        it "returns canonical stored values after tolerated repair" $ do
            let tolerated =
                    either (error . show) id $
                        fromList @(Finite 2) [-1e-17, 1]

            probabilityAt tolerated 0 `shouldBe` 0
            probabilityAt tolerated 1 `shouldBe` 1

    describe "named finite states" $ do
        let namedDistribution =
                either (error . show) id $
                    fromList @NamedState [0.2, 0, 0.8]
            indexedDistribution =
                either (error . show) id $
                    fromList @(Finite 3) [0.2, 0, 0.8]

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
