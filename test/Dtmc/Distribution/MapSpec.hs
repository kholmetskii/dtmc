module Dtmc.Distribution.MapSpec (
    spec,
) where

import Data.Map.Strict qualified as Map
import Dtmc.Distribution (
    Distribution (..),
    DistributionError (..),
 )
import Dtmc.Distribution.Map (
    DistributionMap,
    fromList,
    mapStates,
    pointMass,
    toMap,
 )
import Dtmc.Simplex (
    SimplexError (..),
 )
import Dtmc.TestSupport (
    approxEq,
 )
import Test.Hspec (
    Spec,
    describe,
    expectationFailure,
    it,
    shouldBe,
    shouldSatisfy,
 )

spec :: Spec
spec =
    describe "DistributionMap" $ do
        it "combines duplicates and stores canonical ascending entries" $ do
            let distribution =
                    either (error . show) id $
                        fromList
                            [('b', 0.2), ('a', 0.5), ('b', 0.3), ('c', 0)]
            distributionWeights distribution `shouldBe` [('a', 0.5), ('b', 0.5)]
            support distribution `shouldBe` ['a', 'b']
            Map.toAscList (toMap distribution)
                `shouldBe` [('a', 0.5), ('b', 0.5)]

        it "returns zero for an absent state" $
            probabilityAt (pointMass "present") "absent"
                `shouldBe` 0

        it "pushes weights through a state mapping" $
            case
                ( fromList [(-1, 0.5), (1, 0.5)] ::
                    Either DistributionError (DistributionMap Int)
                )
            of
                Right steps ->
                    distributionWeights (mapStates (+ 10) steps)
                        `shouldBe` [(9, 0.5), (11, 0.5)]
                Left err ->
                    expectationFailure
                        ("expected acceptance, got " <> show err)

        it "combines weights whose states map to the same target" $
            case
                ( fromList [(0, 0.25), (1, 0.25), (2, 0.5)] ::
                    Either DistributionError (DistributionMap Int)
                )
            of
                Right distribution ->
                    distributionWeights (mapStates (`mod` 2) distribution)
                        `shouldBe` [(0, 0.75), (1, 0.25)]
                Left err ->
                    expectationFailure
                        ("expected acceptance, got " <> show err)

        it "rejects an empty law" $
            (fromList [] :: Either DistributionError (DistributionMap Int))
                `shouldSatisfy` either (const True) (const False)

        it "uses the shared error type" $
            fromList ([] :: [(Int, Double)])
                `shouldBe` Left (DistributionError (SumOffBy 0))

        it "removes weights repaired to zero" $
            case fromList [('a', -1e-17), ('b', 1)] of
                Right distribution ->
                    Map.toAscList (toMap distribution)
                        `shouldBe` [('b', 1)]
                Left err ->
                    expectationFailure
                        ("expected acceptance, got " <> show err)

        it "normalises an accepted combined total near one" $
            case fromList [('a', 0.5), ('b', 0.5 - 5e-10)] of
                Right distribution ->
                    approxEq
                        1e-12
                        (sum (Map.elems (toMap distribution)))
                        1
                        `shouldBe` True
                Left err ->
                    expectationFailure
                        ("expected acceptance, got " <> show err)

        it "reports a non-finite combined weight by ascending state index" $
            case fromList [('b', 1), ('a', 0 / 0), ('a', 0)] of
                Left err ->
                    err `shouldBe` DistributionError (NonFiniteEntry 0)
                Right _ ->
                    expectationFailure "expected rejection"
