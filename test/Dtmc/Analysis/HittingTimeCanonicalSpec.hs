{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.HittingTimeCanonicalSpec (
    spec,
) where

import Data.Finite (
    Finite,
    finites,
 )
import Data.Maybe (
    fromMaybe,
 )
import Dtmc.Analysis.Event (
    DiscreteEvent (..),
 )
import Dtmc.Analysis.Expectation (
    Expectation (..),
 )
import Dtmc.Analysis.HittingTime qualified as Hit
import Dtmc.Analysis.ProbabilityOracle qualified as Oracle
import Dtmc.Distribution.Map (
    mkDistributionMap,
 )
import Dtmc.TestSupport (
    approxEq,
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.Transition.Kernel (
    TransitionKernel,
    transitionKernel,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    mkTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
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
    counterexample,
    forAll,
    property,
 )

terminalChain :: TransitionMatrix (Finite 3)
terminalChain =
    checked
        ( mkTransitionMatrix
            ( S.matrix
                [ 0, 0.5, 0.5
                , 0, 0, 1
                , 0, 0, 1
                ] ::
                S.Sq 3
            )
        )

simpleRandomWalk :: TransitionKernel Integer
simpleRandomWalk =
    transitionKernel $ \state ->
        checked
            (mkDistributionMap [(state - 1, 0.5), (state + 1, 0.5)])

tinySurvival :: Double
tinySurvival = 1e-12

tinySurvivalKernel :: TransitionKernel Int
tinySurvivalKernel =
    transitionKernel $ \state ->
        case state of
            0 ->
                checked
                    ( mkDistributionMap
                        [(1, 1 - tinySurvival), (2, tinySurvival)]
                    )
            _ -> deterministicLaw state
  where
    deterministicLaw state =
        checked (mkDistributionMap [(state, 1)])

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

entries :: S.R 3 -> [Double]
entries = LA.toList . S.extract

close :: Double -> Double -> Bool
close = approxEq testTolerance

known :: Maybe Double -> Double
known = fromMaybe (error "oracle horizon does not determine this event")

eventsThrough :: Integer -> [DiscreteEvent]
eventsThrough rawHorizon =
    [EqualTo time | time <- [0 .. horizon]]
        <> [LessThan time | time <- [0 .. horizon + 1]]
        <> [AtMost time | time <- [0 .. horizon]]
        <> [GreaterThan time | time <- [0 .. horizon]]
        <> [AtLeast time | time <- [0 .. horizon + 1]]
  where
    horizon = fromInteger rawHorizon

generatedChecks :: TransitionMatrix (Finite 3) -> Bool
generatedChecks matrix =
    and
        [ let law = Oracle.hittingLaw 4 matrix isTarget initial
              oracle = known (Oracle.lawProbability event law)
              scalar = Hit.probability event matrix isTarget initial
              dense = entries (Hit.probabilityByState event matrix [2])
           in close scalar oracle
                && close (dense !! fromIntegral initial) oracle
        | initial <- finites
        , event <- eventsThrough 4
        ]
  where
    isTarget state = state == (2 :: Finite 3)

spec :: Spec
spec = do
    describe "canonical hitting probability" $ do
        it "implements every relation and carries the infinity atom in upper tails" $ do
            let target state = state == (1 :: Finite 3)
            Hit.probability (EqualTo 0) terminalChain target 0 `shouldBe` 0
            Hit.probability (EqualTo 1) terminalChain target 0 `shouldBe` 0.5
            Hit.probability (LessThan 1) terminalChain target 0 `shouldBe` 0
            Hit.probability (AtMost 1) terminalChain target 0 `shouldBe` 0.5
            Hit.probability (GreaterThan 0) terminalChain target 0 `shouldBe` 1
            Hit.probability (GreaterThan 1) terminalChain target 0 `shouldBe` 0.5
            Hit.probability (AtLeast 0) terminalChain target 0 `shouldBe` 1
            Hit.probability (AtLeast 1) terminalChain target 0 `shouldBe` 1
            Hit.probability (AtLeast 2) terminalChain target 0 `shouldBe` 0.5
            entries (Hit.probabilityByState (GreaterThan 1) terminalChain [1])
                `shouldBe` [0.5, 0, 1]
            entries (Hit.probabilityByState (AtMost 1) terminalChain [1])
                `shouldBe` [0.5, 1, 0]

        it "keeps empty-target and time-zero boundaries structural" $ do
            entries (Hit.probabilityByState (EqualTo 3) terminalChain [])
                `shouldBe` [0, 0, 0]
            entries (Hit.probabilityByState (AtMost 3) terminalChain [])
                `shouldBe` [0, 0, 0]
            entries (Hit.probabilityByState (GreaterThan 3) terminalChain [])
                `shouldBe` [1, 1, 1]
            entries (Hit.probabilityByState (AtLeast 0) terminalChain [1])
                `shouldBe` [1, 1, 1]
            Hit.probability (EqualTo 0) terminalChain (== 1) 1 `shouldBe` 1
            Hit.probability (GreaterThan 0) terminalChain (== 1) 1 `shouldBe` 0

        it "preserves locally finite kernels and tiny survivor mass directly" $ do
            Hit.probability (EqualTo 2) simpleRandomWalk (== 2) 0
                `shouldBe` 0.25
            Hit.probability (AtMost 2) simpleRandomWalk (== 2) 0
                `shouldBe` 0.25
            Hit.probability (GreaterThan 2) simpleRandomWalk (== 2) 0
                `shouldBe` 0.75
            Hit.probability (AtLeast 3) simpleRandomWalk (== 2) 0
                `shouldBe` 0.75
            Hit.probability (GreaterThan 1) tinySurvivalKernel (== 1) 0
                `shouldBe` tinySurvival

        prop "matches the path oracle for every relation (random @3)" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left problem -> counterexample (show problem) False
                    Right matrix -> property (generatedChecks matrix)

    describe "canonical eventual, race, and expectation names" $ do
        it "match the completed defective hitting law" $ do
            let states = finites :: [Finite 3]
            case Hit.eventualProbabilityByState terminalChain [1] of
                Left problem -> error (show problem)
                Right values -> entries values `shouldBe` [0.5, 1, 0]
            mapM_
                ( \(state, expected) ->
                    Hit.eventualProbability terminalChain [1] state
                        `shouldBe` Right expected
                )
                (zip states [0.5, 1, 0])
            case Hit.raceProbabilityByState terminalChain [1] [2] of
                Left problem -> error (show problem)
                Right values -> entries values `shouldBe` [0.5, 1, 0]
            mapM_
                ( \(state, expected) ->
                    Hit.raceProbability terminalChain [1] [2] state
                        `shouldBe` Right expected
                )
                (zip states [0.5, 1, 0])
            Hit.expectationByState terminalChain [1]
                `shouldBe` Right
                    [ InfiniteExpectation
                    , FiniteExpectation 0
                    , InfiniteExpectation
                    ]
            mapM_
                ( \(state, expected) ->
                    Hit.expectation terminalChain [1] state
                        `shouldBe` Right expected
                )
                ( zip
                    states
                    [ InfiniteExpectation
                    , FiniteExpectation 0
                    , InfiniteExpectation
                    ]
                )
