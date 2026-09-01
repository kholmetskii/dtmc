{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.ReturnTimeCanonicalSpec (
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
import Dtmc.Analysis.ProbabilityOracle qualified as Oracle
import Dtmc.Analysis.ReturnTime qualified as Return
import Dtmc.Distribution.Map (
    mkDistributionMap,
 )
import Dtmc.TestSupport
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
                [ 0
                , 0.5
                , 0.5
                , 0
                , 0
                , 1
                , 0
                , 0
                , 1
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

tinyReturnKernel :: TransitionKernel Int
tinyReturnKernel =
    transitionKernel $ \state ->
        case state of
            0 ->
                checked
                    ( mkDistributionMap
                        [(0, 1 - tinySurvival), (1, tinySurvival)]
                    )
            _ -> checked (mkDistributionMap [(state, 1)])

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
        [ let law = Oracle.returnLaw 4 matrix initial
              oracle = known (Oracle.lawProbability event law)
              scalar = Return.probabilityGivenInitialState event matrix initial
              dense = entries (returnProbabilityByState event matrix)
           in close scalar oracle
                && close (dense !! fromIntegral initial) oracle
        | initial <- finites
        , event <- eventsThrough 4
        ]

spec :: Spec
spec = do
    describe "canonical return probability" $ do
        it "enforces the time-zero exclusion exactly" $ do
            entries (returnProbabilityByState (EqualTo 0) terminalChain)
                `shouldBe` [0, 0, 0]
            entries (returnProbabilityByState (LessThan 1) terminalChain)
                `shouldBe` [0, 0, 0]
            entries (returnProbabilityByState (AtMost 0) terminalChain)
                `shouldBe` [0, 0, 0]
            entries (returnProbabilityByState (GreaterThan 0) terminalChain)
                `shouldBe` [1, 1, 1]
            entries (returnProbabilityByState (AtLeast 0) terminalChain)
                `shouldBe` [1, 1, 1]
            entries (returnProbabilityByState (AtLeast 1) terminalChain)
                `shouldBe` [1, 1, 1]

        it "implements every relation and carries non-return mass in upper tails" $ do
            Return.probabilityGivenInitialState (EqualTo 1) terminalChain 2 `shouldBe` 1
            Return.probabilityGivenInitialState (AtMost 1) terminalChain 2 `shouldBe` 1
            Return.probabilityGivenInitialState (GreaterThan 1) terminalChain 2 `shouldBe` 0
            entries (returnProbabilityByState (AtMost 1) terminalChain)
                `shouldBe` [0, 0, 1]
            entries (returnProbabilityByState (GreaterThan 1) terminalChain)
                `shouldBe` [1, 1, 0]
            entries (returnProbabilityByState (AtLeast 2) terminalChain)
                `shouldBe` [1, 1, 0]

        it "preserves locally finite kernels and tiny survivor mass directly" $ do
            Return.probabilityGivenInitialState (EqualTo 2) simpleRandomWalk 0 `shouldBe` 0.5
            Return.probabilityGivenInitialState (AtMost 2) simpleRandomWalk 0 `shouldBe` 0.5
            Return.probabilityGivenInitialState (GreaterThan 2) simpleRandomWalk 0 `shouldBe` 0.5
            Return.probabilityGivenInitialState (AtLeast 3) simpleRandomWalk 0 `shouldBe` 0.5
            Return.probabilityGivenInitialState (GreaterThan 1) tinyReturnKernel 0
                `shouldBe` tinySurvival

        prop "matches the path oracle for every relation (random @3)" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left problem -> counterexample (show problem) False
                    Right matrix -> property (generatedChecks matrix)

    describe "canonical eventual and expectation names" $ do
        it "match the completed defective return laws" $ do
            let states = finites :: [Finite 3]
            case returnEventualProbabilityByState terminalChain of
                Left problem -> error (show problem)
                Right values -> entries values `shouldBe` [0, 0, 1]
            mapM_
                ( \(state, expected) ->
                    Return.eventualProbabilityGivenInitialState terminalChain state
                        `shouldBe` Right expected
                )
                (zip states [0, 0, 1])
            returnExpectationByState terminalChain
                `shouldBe` Right
                    [ InfiniteExpectation
                    , InfiniteExpectation
                    , FiniteExpectation 1
                    ]
            mapM_
                ( \(state, expected) ->
                    Return.expectationGivenInitialState terminalChain state
                        `shouldBe` Right expected
                )
                ( zip
                    states
                    [ InfiniteExpectation
                    , InfiniteExpectation
                    , FiniteExpectation 1
                    ]
                )
