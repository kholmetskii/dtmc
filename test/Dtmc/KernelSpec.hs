{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.KernelSpec (
    spec,
) where

import Control.Monad.ST (
    runST,
 )
import Data.Finite (
    Finite,
    finites,
 )
import Data.List.NonEmpty (
    NonEmpty ((:|)),
 )
import Dtmc.Distribution (
    DistributionVector,
    mkDistributionVector,
 )
import Dtmc.Distribution qualified as Distribution
import Dtmc.Dynamics qualified as Dynamics
import Dtmc.Hitting qualified as Hitting
import Dtmc.Kernel qualified as Kernel
import Dtmc.Probability qualified as Probability
import Dtmc.Simulation qualified as Simulation
import Dtmc.TestSupport (
    approxEq,
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    mkTransitionMatrix,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra.Static qualified as S
import System.Random.MWC qualified as MWC
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
    conjoin,
    counterexample,
    forAll,
    property,
 )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

finiteChain :: TransitionMatrix 3
finiteChain =
    checked $
        mkTransitionMatrix
            ( S.matrix
                [ 0.5
                , 0.5
                , 0
                , 0
                , 0.2
                , 0.8
                , 1
                , 0
                , 0
                ] ::
                S.Sq 3
            )

finiteInitial :: DistributionVector 3
finiteInitial =
    checked (mkDistributionVector (S.vector [0.6, 0.3, 0.1] :: S.R 3))

finiteCycle :: TransitionMatrix 3
finiteCycle =
    checked $
        mkTransitionMatrix
            ( S.matrix
                [ 0
                , 1
                , 0
                , 0
                , 0
                , 1
                , 1
                , 0
                , 0
                ] ::
                S.Sq 3
            )

asSparseKernel ::
    forall n.
    (KnownNat n) =>
    TransitionMatrix n ->
    Kernel.TransitionKernel (Finite n)
asSparseKernel matrix =
    Kernel.transitionKernel $ \source ->
        checked $
            Distribution.mkSparseDistribution
                [ (destination, Probability.transitionProbability matrix source destination)
                | destination <- finites
                ]

asSparseDistribution ::
    forall n.
    (KnownNat n) =>
    DistributionVector n ->
    Distribution.SparseDistribution (Finite n)
asSparseDistribution distribution =
    checked $
        Distribution.mkSparseDistribution
            [ (state, Distribution.probabilityAt distribution state)
            | state <- finites
            ]

sparseChain :: Kernel.TransitionKernel (Finite 3)
sparseChain = asSparseKernel finiteChain

sparseInitial :: Distribution.SparseDistribution (Finite 3)
sparseInitial = asSparseDistribution finiteInitial

simpleRandomWalk :: Kernel.TransitionKernel Integer
simpleRandomWalk =
    Kernel.transitionKernel $ \state ->
        checked (Distribution.mkSparseDistribution [(state - 1, 0.5), (state + 1, 0.5)])

closeTo :: Double -> Double -> Bool
closeTo = approxEq testTolerance

rightCloseTo :: Either error Double -> Either error Double -> Bool
rightCloseTo (Right left) (Right right) = closeTo left right
rightCloseTo (Left _) (Left _) = True
rightCloseTo _ _ = False

spec :: Spec
spec = do
    describe "shared MarkovKernel interface" $ do
        it "exposes a finite matrix row as a validated finite-support law" $
            Distribution.distributionWeights (Kernel.transitionLaw finiteChain 1)
                `shouldBe` [(1, 0.2), (2, 0.8)]

        it "converts a dense finite initial law without changing its weights" $
            Distribution.distributionWeights (Distribution.toSparseDistribution finiteInitial)
                `shouldBe` [(0, 0.6), (1, 0.3), (2, 0.1)]

        it "runs one probability function over both kernel representations" $ do
            let finiteResult =
                    Probability.probabilityAtTime
                        4
                        finiteInitial
                        finiteChain
                        2
                sparseResult =
                    Probability.probabilityAtTime
                        4
                        sparseInitial
                        sparseChain
                        2
            sparseResult `shouldSatisfy` closeTo finiteResult

        it "runs shared bounded stopping queries directly on a finite matrix" $ do
            Hitting.hittingTimeProbabilityBefore finiteChain (== 2) 0 4
                `shouldSatisfy` closeTo 0.68
            Hitting.returnTimeProbabilityAt finiteChain 1 1
                `shouldSatisfy` closeTo 0.2

        it "simulates a finite matrix through the same interface" $
            let trajectory = runST $ do
                    generator <- MWC.create
                    Simulation.simulateN 4 finiteCycle 0 generator
             in trajectory `shouldBe` [0, 1, 2, 0, 1]

        prop "gives the same transition powers for both representations" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left problem -> counterexample (show problem) False
                    Right matrix ->
                        let kernel = asSparseKernel matrix
                         in conjoin
                                [ counterexample (show (source, destination, time)) $
                                    property $
                                        closeTo
                                            (Probability.transitionProbabilityN time matrix source destination)
                                            (Probability.transitionProbabilityN time kernel source destination)
                                | source <- finites :: [Finite 3]
                                , destination <- finites :: [Finite 3]
                                , time <- [0 .. 4]
                                ]

    describe "locally finite dynamics" $ do
        it "evolves an infinite-state random walk without enumerating its state space" $
            Distribution.distributionWeights
                (Dynamics.evolveSparseN 2 (Distribution.pointMass 0) simpleRandomWalk)
                `shouldBe` [(-2, 0.25), (0, 0.5), (2, 0.25)]

        it "agrees with finite matrix dynamics at several horizons" $
            sequence_
                [ Distribution.probabilityAt
                    (Dynamics.evolveSparseN time sparseInitial sparseChain)
                    state
                    `shouldSatisfy` closeTo
                        ( Distribution.probabilityAt
                            (Dynamics.evolveN time finiteInitial finiteChain)
                            state
                        )
                | time <- [0 .. 4]
                , state <- finites :: [Finite 3]
                ]

    describe "finite-horizon probabilities" $ do
        it "computes transition probabilities on an infinite state type" $ do
            Probability.transitionProbabilityN 2 simpleRandomWalk 0 0
                `shouldSatisfy` closeTo 0.5
            Probability.transitionProbabilityN 3 simpleRandomWalk 0 0
                `shouldBe` 0

        it "matches finite trajectory and observation queries" $ do
            Probability.pathProbability sparseInitial sparseChain (0 :| [1, 2])
                `shouldSatisfy` closeTo
                    (Probability.pathProbability finiteInitial finiteChain (0 :| [1, 2]))
            Probability.probability
                sparseInitial
                sparseChain
                [Probability.At 3 2, Probability.At 0 0, Probability.At 1 1]
                `shouldSatisfy` closeTo
                    ( Probability.probability
                        finiteInitial
                        finiteChain
                        [Probability.At 3 2, Probability.At 0 0, Probability.At 1 1]
                    )

        it "matches finite conditional probability queries" $
            rightCloseTo
                ( Probability.conditionalProbability
                    sparseInitial
                    sparseChain
                    [Probability.At 2 2]
                    [Probability.At 0 0]
                )
                ( Probability.conditionalProbability
                    finiteInitial
                    finiteChain
                    [Probability.At 2 2]
                    [Probability.At 0 0]
                )
                `shouldBe` True

    describe "bounded hitting and first return" $ do
        it "uses strict hitting bounds on the infinite random walk" $ do
            Hitting.hittingTimeProbabilityAt simpleRandomWalk (== 2) 0 2
                `shouldSatisfy` closeTo 0.25
            Hitting.hittingTimeProbabilityBefore simpleRandomWalk (== 2) 0 2
                `shouldBe` 0
            Hitting.hittingTimeProbabilityBefore simpleRandomWalk (== 2) 0 3
                `shouldSatisfy` closeTo 0.25

        it "distinguishes return time from time-zero hitting" $ do
            Hitting.hittingTimeProbabilityAt simpleRandomWalk (== 0) 0 0
                `shouldBe` 1
            Hitting.returnTimeProbabilityAt simpleRandomWalk 0 0
                `shouldBe` 0
            Hitting.returnTimeProbabilityAt simpleRandomWalk 0 2
                `shouldSatisfy` closeTo 0.5
            Hitting.returnTimeProbabilityBefore simpleRandomWalk 0 2
                `shouldBe` 0
            Hitting.returnTimeProbabilityBefore simpleRandomWalk 0 3
                `shouldSatisfy` closeTo 0.5

        prop "matches finite exact and bounded hitting/return queries" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left problem -> counterexample (show problem) False
                    Right matrix ->
                        let kernel = asSparseKernel matrix
                            target state = state == (2 :: Finite 3)
                         in conjoin
                                [ counterexample (show (state, time)) $
                                    property $
                                        and
                                            [ closeTo
                                                (Hitting.hittingTimeProbabilityAt matrix target state time)
                                                (Hitting.hittingTimeProbabilityAt kernel target state time)
                                            , closeTo
                                                (Hitting.hittingTimeProbabilityBefore matrix target state time)
                                                (Hitting.hittingTimeProbabilityBefore kernel target state time)
                                            , closeTo
                                                (Hitting.returnTimeProbabilityAt matrix state time)
                                                (Hitting.returnTimeProbabilityAt kernel state time)
                                            , closeTo
                                                (Hitting.returnTimeProbabilityBefore matrix state time)
                                                (Hitting.returnTimeProbabilityBefore kernel state time)
                                            ]
                                | state <- finites :: [Finite 3]
                                , time <- [0 .. 4]
                                ]

    describe "simulation" $
        it "returns the initial state plus exactly the requested transitions" $
            let trajectory = runST $ do
                    generator <- MWC.create
                    Simulation.simulateN
                        4
                        (Kernel.deterministicKernel (\state -> (state + 1) `mod` (3 :: Int)))
                        0
                        generator
             in trajectory `shouldBe` [0, 1, 2, 0, 1]
