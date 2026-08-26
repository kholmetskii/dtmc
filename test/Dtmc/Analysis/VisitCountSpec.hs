{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.VisitCountSpec (
    spec,
) where

import Data.Finite (
    Finite,
    finites,
 )
import Dtmc.Analysis.FixedTime (
    probabilityAtTime,
    transitionProbability,
 )
import Dtmc.Analysis.VisitCount (
    visitCountDistributionBefore,
    visitCountExpectationBefore,
    visitCountProbabilityBefore,
 )
import Dtmc.Distribution (
    distributionWeights,
 )
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.TestSupport (
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.Transition.Kernel qualified as Kernel
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    mkTransitionMatrix,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
    shouldSatisfy,
 )
import Test.Hspec.QuickCheck (
    prop,
 )
import Test.QuickCheck (
    choose,
    conjoin,
    counterexample,
    forAll,
    property,
 )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

twoCycle :: TransitionMatrix (Finite 2)
twoCycle =
    checked $
        mkTransitionMatrix
            ( S.matrix
                [ 0
                , 1
                , 1
                , 0
                ] ::
                S.Sq 2
            )

mixedInitial :: DistributionMap.DistributionMap (Finite 2)
mixedInitial =
    checked (DistributionMap.mkDistributionMap [(0, 0.25), (1, 0.75)])

asKernel ::
    TransitionMatrix (Finite 2) ->
    Kernel.TransitionKernel (Finite 2)
asKernel matrix =
    Kernel.transitionKernel $ \source ->
        checked $
            DistributionMap.mkDistributionMap
                [ (destination, transitionProbability matrix source destination)
                | destination <- finites
                ]

simpleRandomWalk :: Kernel.TransitionKernel Integer
simpleRandomWalk =
    Kernel.transitionKernel $ \state ->
        checked $
            DistributionMap.mkDistributionMap
                [(state - 1, 0.5), (state + 1, 0.5)]

closeTo :: Double -> Double -> Bool
closeTo expected actual = abs (actual - expected) <= testTolerance

spec :: Spec
spec = do
    describe "visitCountDistributionBefore" $ do
        it "is a point mass at zero for bound zero" $
            distributionWeights
                ( visitCountDistributionBefore
                    0
                    (DistributionMap.pointMass (0 :: Finite 2))
                    twoCycle
                    (== 0)
                )
                `shouldBe` [(0, 1)]

        it "counts the initial state at a positive bound" $
            distributionWeights
                (visitCountDistributionBefore 1 mixedInitial twoCycle (== 0))
                `shouldBe` [(0, 0.75), (1, 0.25)]

        it "counts deterministic visits at times zero through bound minus one" $ do
            let initial = DistributionMap.pointMass (0 :: Finite 2)
            distributionWeights (visitCountDistributionBefore 1 initial twoCycle (== 0))
                `shouldBe` [(1, 1)]
            distributionWeights (visitCountDistributionBefore 2 initial twoCycle (== 0))
                `shouldBe` [(1, 1)]
            distributionWeights (visitCountDistributionBefore 3 initial twoCycle (== 0))
                `shouldBe` [(2, 1)]

        it "computes an exact law on an infinite random walk" $
            distributionWeights
                ( visitCountDistributionBefore
                    3
                    (DistributionMap.pointMass (0 :: Integer))
                    simpleRandomWalk
                    (== 0)
                )
                `shouldBe` [(1, 0.5), (2, 0.5)]

        prop "has total mass one and no count above the bound (random @3)" $
            forAll (choose (0, 5 :: Int)) $ \rawBound ->
                forAll (genTransitionMatrix @3) $ \rawMatrix ->
                    case mkTransitionMatrix rawMatrix of
                        Left err -> counterexample (show err) False
                        Right matrix ->
                            let bound = fromIntegral rawBound
                                law =
                                    visitCountDistributionBefore
                                        bound
                                        (DistributionMap.pointMass (0 :: Finite 3))
                                        matrix
                                        (== 0)
                                weights = distributionWeights law
                             in conjoin
                                    [ property (closeTo 1 (sum (map snd weights)))
                                    , property (all ((<= bound) . fst) weights)
                                    ]

    describe "visitCountProbabilityBefore" $ do
        it "looks up one coordinate of the count distribution" $ do
            let initial = DistributionMap.pointMass (0 :: Integer)
            visitCountProbabilityBefore 3 initial simpleRandomWalk (== 0) 1
                `shouldSatisfy` closeTo 0.5
            visitCountProbabilityBefore 3 initial simpleRandomWalk (== 0) 2
                `shouldSatisfy` closeTo 0.5
            visitCountProbabilityBefore 3 initial simpleRandomWalk (== 0) 3
                `shouldBe` 0

        it "agrees for a matrix and its equivalent kernel" $ do
            let initial = DistributionMap.pointMass (0 :: Finite 2)
                kernel = asKernel twoCycle
            sequence_
                [ visitCountProbabilityBefore bound initial twoCycle (== 0) count
                    `shouldBe` visitCountProbabilityBefore bound initial kernel (== 0) count
                | bound <- [0 .. 5]
                , count <- [0 .. bound]
                ]

    describe "visitCountExpectationBefore" $ do
        it "is the expectation of the random-walk count law" $
            visitCountExpectationBefore
                3
                (DistributionMap.pointMass (0 :: Integer))
                simpleRandomWalk
                (== 0)
                `shouldSatisfy` closeTo 1.5

        prop "equals the sum of finite-time visit probabilities (random @3)" $
            forAll (choose (0, 5 :: Int)) $ \rawBound ->
                forAll (genTransitionMatrix @3) $ \rawMatrix ->
                    case mkTransitionMatrix rawMatrix of
                        Left err -> counterexample (show err) False
                        Right matrix ->
                            let bound = fromIntegral rawBound
                                initial = DistributionMap.pointMass (0 :: Finite 3)
                                expectation =
                                    visitCountExpectationBefore bound initial matrix (== 0)
                                marginalSum =
                                    sum
                                        [ probabilityAtTime (fromIntegral time) initial matrix 0
                                        | time <- [0 .. rawBound - 1]
                                        ]
                             in property (closeTo marginalSum expectation)
