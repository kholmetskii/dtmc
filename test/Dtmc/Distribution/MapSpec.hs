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
    mkDistributionMap,
    pointMass,
    unDistributionMap,
 )
import Dtmc.Simplex (
    SimplexError (..),
 )
import Test.Hspec (
    Spec,
    describe,
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
