{-# LANGUAGE DataKinds #-}

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
    fromList,
 )
import Dtmc.TestSupport
import Dtmc.Transition.Kernel (
    TransitionKernel,
    fromLaws,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    TransitionMatrixError,
    fromRows,
 )
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
        ( fromRows
            ( chunksOf
                3
                [ 0
                , 0.5
                , 0.5
                , 0
                , 0
                , 1
                , 0
                , 0
                , 1
                ]
            )
        )

simpleRandomWalk :: TransitionKernel Integer
simpleRandomWalk =
    fromLaws $ \state ->
        checked
            (fromList [(state - 1, 0.5), (state + 1, 0.5)])

tinySurvival :: Double
tinySurvival = 1e-12

tinySurvivalKernel :: TransitionKernel Int
tinySurvivalKernel =
    fromLaws $ \state ->
        case state of
            0 ->
                checked
                    ( fromList
                        [(1, 1 - tinySurvival), (2, tinySurvival)]
                    )
            _ -> deterministicLaw state
  where
    deterministicLaw state =
        checked (fromList [(state, 1)])

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

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
              scalar = Hit.probabilityGivenInitialState event matrix isTarget initial
              dense = hitProbabilityByState event matrix [2]
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
            Hit.probabilityGivenInitialState (EqualTo 0) terminalChain target 0 `shouldBe` 0
            Hit.probabilityGivenInitialState (EqualTo 1) terminalChain target 0 `shouldBe` 0.5
            Hit.probabilityGivenInitialState (LessThan 1) terminalChain target 0 `shouldBe` 0
            Hit.probabilityGivenInitialState (AtMost 1) terminalChain target 0 `shouldBe` 0.5
            Hit.probabilityGivenInitialState (GreaterThan 0) terminalChain target 0 `shouldBe` 1
            Hit.probabilityGivenInitialState (GreaterThan 1) terminalChain target 0 `shouldBe` 0.5
            Hit.probabilityGivenInitialState (AtLeast 0) terminalChain target 0 `shouldBe` 1
            Hit.probabilityGivenInitialState (AtLeast 1) terminalChain target 0 `shouldBe` 1
            Hit.probabilityGivenInitialState (AtLeast 2) terminalChain target 0 `shouldBe` 0.5
            hitProbabilityByState (GreaterThan 1) terminalChain [1]
                `shouldBe` [0.5, 0, 1]
            hitProbabilityByState (AtMost 1) terminalChain [1]
                `shouldBe` [0.5, 1, 0]

        it "keeps empty-target and time-zero boundaries structural" $ do
            hitProbabilityByState (EqualTo 3) terminalChain []
                `shouldBe` [0, 0, 0]
            hitProbabilityByState (AtMost 3) terminalChain []
                `shouldBe` [0, 0, 0]
            hitProbabilityByState (GreaterThan 3) terminalChain []
                `shouldBe` [1, 1, 1]
            hitProbabilityByState (AtLeast 0) terminalChain [1]
                `shouldBe` [1, 1, 1]
            Hit.probabilityGivenInitialState (EqualTo 0) terminalChain (== 1) 1 `shouldBe` 1
            Hit.probabilityGivenInitialState (GreaterThan 0) terminalChain (== 1) 1 `shouldBe` 0

        it "preserves locally finite kernels and tiny survivor mass directly" $ do
            Hit.probabilityGivenInitialState (EqualTo 2) simpleRandomWalk (== 2) 0
                `shouldBe` 0.25
            Hit.probabilityGivenInitialState (AtMost 2) simpleRandomWalk (== 2) 0
                `shouldBe` 0.25
            Hit.probabilityGivenInitialState (GreaterThan 2) simpleRandomWalk (== 2) 0
                `shouldBe` 0.75
            Hit.probabilityGivenInitialState (AtLeast 3) simpleRandomWalk (== 2) 0
                `shouldBe` 0.75
            Hit.probabilityGivenInitialState (GreaterThan 1) tinySurvivalKernel (== 1) 0
                `shouldBe` tinySurvival

        prop "matches the path oracle for every relation (random @3)" $
            forAll (genTransitionRows 3) $ \rawMatrix ->
                case fromRows rawMatrix ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 3)) of
                    Left problem -> counterexample (show problem) False
                    Right matrix -> property (generatedChecks matrix)

    describe "canonical eventual, race, and expectation names" $ do
        it "match the completed defective hitting law" $ do
            let states = finites :: [Finite 3]
            case hitEventualProbabilityByState terminalChain [1] of
                Left problem -> error (show problem)
                Right values -> values `shouldBe` [0.5, 1, 0]
            mapM_
                ( \(state, expected) ->
                    Hit.eventualProbabilityGivenInitialState terminalChain [1] state
                        `shouldBe` Right expected
                )
                (zip states [0.5, 1, 0])
            case hitRaceProbabilityByState terminalChain [1] [2] of
                Left problem -> error (show problem)
                Right values -> values `shouldBe` [0.5, 1, 0]
            mapM_
                ( \(state, expected) ->
                    Hit.raceProbabilityGivenInitialState terminalChain [1] [2] state
                        `shouldBe` Right expected
                )
                (zip states [0.5, 1, 0])
            hitExpectationByState terminalChain [1]
                `shouldBe` Right
                    [ InfiniteExpectation
                    , FiniteExpectation 0
                    , InfiniteExpectation
                    ]
            mapM_
                ( \(state, expected) ->
                    Hit.expectationGivenInitialState terminalChain [1] state
                        `shouldBe` Right expected
                )
                ( zip
                    states
                    [ InfiniteExpectation
                    , FiniteExpectation 0
                    , InfiniteExpectation
                    ]
                )
