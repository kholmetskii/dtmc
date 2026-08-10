{-# LANGUAGE TypeApplications #-}

module Dtmc.Distribution.InterfaceSpec (
    spec,
) where

import Data.Finite (
    Finite,
 )
import Dtmc.Distribution (
    Distribution (..),
    DistributionError,
 )
import Dtmc.Distribution.Map (
    DistributionMap,
    mkDistributionMap,
    toDistributionMap,
 )
import Dtmc.Distribution.Vector (
    mkDistributionVector,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
 )

spec :: Spec
spec =
    describe "Distribution interface" $ do
        let vector =
                either (error . show) id $
                    mkDistributionVector @(Finite 3) (S.vector [0.2, 0, 0.8] :: S.R 3)
            mapDistribution =
                either (error . show) id
                    ( mkDistributionMap [(0, 0.2), (2, 0.8)] ::
                        Either DistributionError (DistributionMap (Finite 3))
                    )

        it "exposes the same weights and support for both representations" $ do
            distributionWeights vector `shouldBe` distributionWeights mapDistribution
            support vector `shouldBe` support mapDistribution

        it "converts both representations to the same canonical map" $ do
            toDistributionMap vector `shouldBe` mapDistribution
            toDistributionMap mapDistribution `shouldBe` mapDistribution

        it "converts a dense law without changing its weights" $
            distributionWeights (toDistributionMap vector)
                `shouldBe` [(0, 0.2), (2, 0.8)]
