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
    fromList,
    fromDistribution,
 )
import Dtmc.Distribution.Vector.HMatrix (
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
                either
                    (error . show)
                    id
                    ( fromList [(0, 0.2), (2, 0.8)] ::
                        Either DistributionError (DistributionMap (Finite 3))
                    )

        it "exposes the same weights and support for both representations" $ do
            distributionWeights vector `shouldBe` distributionWeights mapDistribution
            support vector `shouldBe` support mapDistribution

        it "converts both representations to the same canonical map" $ do
            fromDistribution vector `shouldBe` mapDistribution
            fromDistribution mapDistribution `shouldBe` mapDistribution

        it "converts a dense law without changing its weights" $
            distributionWeights (fromDistribution vector)
                `shouldBe` [(0, 0.2), (2, 0.8)]
