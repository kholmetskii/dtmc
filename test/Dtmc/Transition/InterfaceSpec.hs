{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Transition.InterfaceSpec (
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
import Dtmc.Distribution qualified as Distribution
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.Distribution.Vector (
    DistributionVector,
    mkDistributionVector,
 )
import Dtmc.Dynamics qualified as Dynamics
import Dtmc.Hitting qualified as Hitting
import Dtmc.Probability qualified as Probability
import Dtmc.Simulation qualified as Simulation
import Dtmc.State qualified as State
import Dtmc.TestSupport (
    approxEq,
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.Transition qualified as Transition
import Dtmc.Transition.Kernel qualified as Kernel
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    mkTransitionMatrix,
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

finiteChain :: TransitionMatrix (Finite 3)
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

finiteInitial :: DistributionVector (Finite 3)
finiteInitial =
    checked (mkDistributionVector (S.vector [0.6, 0.3, 0.1] :: S.R 3))

finiteCycle :: TransitionMatrix (Finite 3)
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

asTransitionKernel ::
    forall state.
    (State.FiniteState state) =>
    TransitionMatrix state ->
    Kernel.TransitionKernel state
asTransitionKernel matrix =
    Kernel.transitionKernel $ \source ->
        checked $
            DistributionMap.mkDistributionMap
                [ (destination, Probability.transitionProbability matrix source destination)
                | destination <- State.finiteStates
                ]

asDistributionMap ::
    forall state.
    (State.FiniteState state) =>
    DistributionVector state ->
    DistributionMap.DistributionMap state
asDistributionMap distribution =
    checked $
        DistributionMap.mkDistributionMap
            [ (state, Distribution.probabilityAt distribution state)
            | state <- State.finiteStates
            ]

kernelChain :: Kernel.TransitionKernel (Finite 3)
kernelChain = asTransitionKernel finiteChain

mapInitial :: DistributionMap.DistributionMap (Finite 3)
mapInitial = asDistributionMap finiteInitial

closeTo :: Double -> Double -> Bool
closeTo = approxEq testTolerance

rightCloseTo :: Either error Double -> Either error Double -> Bool
rightCloseTo (Right left) (Right right) = closeTo left right
rightCloseTo (Left _) (Left _) = True
rightCloseTo _ _ = False

spec :: Spec
spec = do
    describe "shared Transition interface" $ do
        it "exposes a finite matrix row as a validated finite-support law" $
            Distribution.distributionWeights (Transition.transitionLaw finiteChain 1)
                `shouldBe` [(1, 0.2), (2, 0.8)]

        it "converts a dense finite initial law without changing its weights" $
            Distribution.distributionWeights (DistributionMap.toDistributionMap finiteInitial)
                `shouldBe` [(0, 0.6), (1, 0.3), (2, 0.1)]

        it "runs one probability function over both kernel representations" $ do
            let finiteResult =
                    Probability.probabilityAtTime
                        4
                        finiteInitial
                        finiteChain
                        2
                kernelResult =
                    Probability.probabilityAtTime
                        4
                        mapInitial
                        kernelChain
                        2
            kernelResult `shouldSatisfy` closeTo finiteResult

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
                        let kernel = asTransitionKernel matrix
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

    describe "cross-representation dynamics" $ do
        it "agrees with finite matrix dynamics at several horizons" $
            sequence_
                [ Distribution.probabilityAt
                    (Dynamics.evolveN time mapInitial kernelChain)
                    state
                    `shouldSatisfy` closeTo
                        ( Distribution.probabilityAt
                            (Dynamics.evolveVectorN time finiteInitial finiteChain)
                            state
                        )
                | time <- [0 .. 4]
                , state <- finites :: [Finite 3]
                ]

    describe "cross-representation finite-horizon probabilities" $ do
        it "matches finite trajectory and observation queries" $ do
            Probability.pathProbability mapInitial kernelChain (0 :| [1, 2])
                `shouldSatisfy` closeTo
                    (Probability.pathProbability finiteInitial finiteChain (0 :| [1, 2]))
            Probability.probability
                mapInitial
                kernelChain
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
                    mapInitial
                    kernelChain
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

    describe "cross-representation bounded hitting and first return" $ do
        prop "matches finite exact and bounded hitting/return queries" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left problem -> counterexample (show problem) False
                    Right matrix ->
                        let kernel = asTransitionKernel matrix
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
