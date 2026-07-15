{-# LANGUAGE TypeApplications #-}

module Dtmc.ApproxSpec (
    spec,
) where

import Dtmc.Approx (
    approxEq,
    approxEqDist,
    approxEqMatrix,
    approxEqR,
    tolerance,
 )
import Dtmc.Distribution (
    Distribution,
    mkDistribution,
 )
import Dtmc.TestSupport (
    approxDistributionEq,
    approxTransitionMatrixEq,
    genSimplexPoint,
    genTransitionMatrix,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    mkTransitionMatrix,
 )
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
    choose,
    counterexample,
    forAll,
    (===),
 )

fromValid :: (Show e) => Either e a -> a
fromValid = either (error . show) id

uniformTwo :: Distribution 2
uniformTwo = fromValid (mkDistribution (S.vector [0.5, 0.5]))

shiftedTwo :: Double -> Distribution 2
shiftedTwo eps =
    fromValid (mkDistribution (S.vector [0.5 + eps, 0.5 - eps]))

twoState :: TransitionMatrix 2
twoState =
    fromValid (mkTransitionMatrix (S.matrix [0.9, 0.1, 0.4, 0.6]))

twoStateShifted :: Double -> TransitionMatrix 2
twoStateShifted eps =
    fromValid
        (mkTransitionMatrix (S.matrix [0.9 + eps, 0.1 - eps, 0.4, 0.6]))

spec :: Spec
spec = do
    describe "approxEq" $ do
        it "accepts a difference of exactly the tolerance" $
            approxEq 0 tolerance `shouldBe` True

        it "rejects a difference beyond the tolerance" $
            approxEq 0 (2 * tolerance) `shouldBe` False

        prop "is reflexive on probabilities" $
            forAll (choose (0, 1 :: Double)) $ \x ->
                approxEq x x

    describe "approxEqR" $ do
        it "accepts vectors differing within tolerance" $
            approxEqR
                (S.vector [0.5, 0.5] :: S.R 2)
                (S.vector [0.5 + 1e-10, 0.5 - 1e-10])
                `shouldBe` True

        it "rejects vectors differing beyond tolerance" $
            approxEqR
                (S.vector [0.5, 0.5] :: S.R 2)
                (S.vector [0.5 + 1e-6, 0.5 - 1e-6])
                `shouldBe` False

    describe "approxEqDist" $ do
        it "accepts distributions differing within tolerance" $
            approxEqDist uniformTwo (shiftedTwo 1e-10) `shouldBe` True

        it "rejects distributions differing beyond tolerance" $
            approxEqDist uniformTwo (shiftedTwo 1e-6) `shouldBe` False

        prop "agrees with the independent test comparator at tolerance" $
            forAll ((,) <$> genSimplexPoint 3 <*> genSimplexPoint 3) $
                \(xs, ys) ->
                    case ( mkDistribution (S.vector xs :: S.R 3)
                         , mkDistribution (S.vector ys :: S.R 3)
                         ) of
                        (Right a, Right b) ->
                            approxEqDist a b
                                === approxDistributionEq tolerance a b
                        result ->
                            counterexample
                                ("generated vectors were rejected: " <> show result)
                                False

    describe "approxEqMatrix" $ do
        it "accepts matrices differing within tolerance" $
            approxEqMatrix twoState (twoStateShifted 1e-10) `shouldBe` True

        it "rejects matrices differing beyond tolerance" $
            approxEqMatrix twoState (twoStateShifted 1e-6) `shouldBe` False

        prop "agrees with the independent test comparator at tolerance" $
            forAll ((,) <$> genTransitionMatrix @3 <*> genTransitionMatrix @3) $
                \(m1, m2) ->
                    case (mkTransitionMatrix m1, mkTransitionMatrix m2) of
                        (Right a, Right b) ->
                            approxEqMatrix a b
                                === approxTransitionMatrixEq tolerance a b
                        result ->
                            counterexample
                                ("generated matrices were rejected: " <> show result)
                                False
