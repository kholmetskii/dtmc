{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.VisitCountCanonicalSpec (
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
import Dtmc.Analysis.VisitCount qualified as Visit
import Dtmc.Distribution (
    distributionWeights,
 )
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.TestSupport
import Dtmc.Transition.Kernel (
    TransitionKernel,
    transitionKernel,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
 )
import Dtmc.Transition.Matrix.HMatrix (
    mkTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )
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
    counterexample,
    forAll,
    property,
 )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

transientVisitChain :: TransitionMatrix (Finite 3)
transientVisitChain =
    checked
        ( mkTransitionMatrix
            ( S.matrix
                [ 1 / 4
                , 0
                , 3 / 4
                , 1 / 2
                , 0
                , 1 / 2
                , 0
                , 0
                , 1
                ] ::
                S.Sq 3
            )
        )

recurrentVisitChain :: TransitionMatrix (Finite 4)
recurrentVisitChain =
    checked
        ( mkTransitionMatrix
            ( S.matrix
                [ 0
                , 1 / 2
                , 1 / 2
                , 0
                , 1 / 2
                , 0
                , 0
                , 1 / 2
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 1
                ] ::
                S.Sq 4
            )
        )

tinyReturn :: Double
tinyReturn = 1e-12

tinyVisitChain :: TransitionMatrix (Finite 2)
tinyVisitChain =
    checked
        ( mkTransitionMatrix
            ( S.matrix
                [ tinyReturn
                , 1 - tinyReturn
                , 0
                , 1
                ] ::
                S.Sq 2
            )
        )

simpleRandomWalk :: TransitionKernel Integer
simpleRandomWalk =
    transitionKernel $ \state ->
        checked
            ( DistributionMap.fromList
                [(state - 1, 0.5), (state + 1, 0.5)]
            )

close :: Double -> Double -> Bool
close = approxEq testTolerance

entries :: S.R 3 -> [Double]
entries = LA.toList . S.extract

known :: Maybe Double -> Double
known = fromMaybe (error "oracle horizon does not determine this event")

eventsThrough :: Natural -> [DiscreteEvent]
eventsThrough horizon =
    [EqualTo count | count <- [0 .. horizon]]
        <> [LessThan count | count <- [0 .. horizon + 1]]
        <> [AtMost count | count <- [0 .. horizon]]
        <> [GreaterThan count | count <- [0 .. horizon]]
        <> [AtLeast count | count <- [0 .. horizon + 1]]

generatedTotalChecks :: TransitionMatrix (Finite 3) -> Bool
generatedTotalChecks matrix =
    and
        [ let scalar = checked (Visit.totalProbabilityGivenInitialState event matrix 0 initial)
              dense = entries (checked (visitTotalProbabilityByState event matrix 0))
           in close (dense !! fromIntegral initial) scalar
                && scalar >= negate testTolerance
                && scalar <= 1 + testTolerance
        | initial <- finites
        , event <- eventsThrough 4
        ]

generatedBoundedChecks :: TransitionMatrix (Finite 3) -> Bool
generatedBoundedChecks matrix =
    and
        [ let initial = DistributionMap.pointMass (0 :: Finite 3)
              oracleLaw = Oracle.visitLawBefore bound [(0, 1)] matrix (== 0)
              expected = known (Oracle.lawProbability event oracleLaw)
              actual = Visit.boundedProbability bound event initial matrix (== 0)
           in close actual expected
        | bound <- [0 .. 4]
        , event <- eventsThrough bound
        ]

spec :: Spec
spec = do
    describe "canonical total visit count" $ do
        it "implements every relation for a transient geometric law" $ do
            let probability event =
                    checked (Visit.totalProbabilityGivenInitialState event transientVisitChain 0 1)
            probability (EqualTo 0) `shouldSatisfy` close (1 / 2)
            probability (EqualTo 1) `shouldSatisfy` close (3 / 8)
            probability (LessThan 2) `shouldSatisfy` close (7 / 8)
            probability (AtMost 1) `shouldSatisfy` close (7 / 8)
            probability (GreaterThan 1) `shouldSatisfy` close (1 / 8)
            probability (AtLeast 2) `shouldSatisfy` close (1 / 8)
            probability (AtLeast 0) `shouldBe` 1

        it "places recurrent positive-count mass structurally at infinity" $ do
            let probabilities event =
                    LA.toList
                        (S.extract (checked (visitTotalProbabilityByState event recurrentVisitChain 2)))
                expectedHit = [2 / 3, 1 / 3, 1, 0]
                expectedMiss = [1 / 3, 2 / 3, 0, 1]
            sequence_
                [ actual `shouldSatisfy` close expected
                | (actual, expected) <- zip (probabilities (GreaterThan 3)) expectedHit
                ]
            sequence_
                [ actual `shouldSatisfy` close expected
                | (actual, expected) <- zip (probabilities (AtMost 3)) expectedMiss
                ]
            probabilities (EqualTo 2) `shouldBe` [0, 0, 0, 0]
            probabilities (AtLeast 0) `shouldBe` [1, 1, 1, 1]

        it "evaluates a tiny upper tail without complement subtraction" $ do
            let actual =
                    checked
                        ( Visit.totalProbabilityGivenInitialState
                            (GreaterThan 1)
                            tinyVisitChain
                            0
                            0
                        )
            actual `shouldSatisfy` (> 0)
            actual `shouldSatisfy` (\value -> abs (value - tinyReturn) < 1e-15)

        prop "keeps scalar and all-state event queries consistent (random @3)" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left problem -> counterexample (show problem) False
                    Right matrix -> property (generatedTotalChecks matrix)

    describe "canonical bounded visit count" $ do
        it "supports every event relation on a locally finite kernel" $ do
            let initial = DistributionMap.pointMass (0 :: Integer)
                probability event =
                    Visit.boundedProbability 3 event initial simpleRandomWalk (== 0)
            distributionWeights (Visit.boundedLaw 3 initial simpleRandomWalk (== 0))
                `shouldBe` [(1, 0.5), (2, 0.5)]
            probability (EqualTo 1) `shouldBe` 0.5
            probability (LessThan 2) `shouldBe` 0.5
            probability (AtMost 1) `shouldBe` 0.5
            probability (GreaterThan 1) `shouldBe` 0.5
            probability (AtLeast 2) `shouldBe` 0.5
            Visit.boundedExpectation 3 initial simpleRandomWalk (== 0)
                `shouldBe` 1.5

        prop "matches independent path enumeration for every relation (random @3)" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left problem -> counterexample (show problem) False
                    Right matrix -> property (generatedBoundedChecks matrix)

    describe "canonical infinite and expectation names" $ do
        it "match the completed total-visit law" $ do
            let infiniteValues =
                    LA.toList
                        ( S.extract
                            ( checked
                                (visitInfiniteProbabilityByState recurrentVisitChain 2)
                            )
                        )
            sequence_
                [ actual `shouldSatisfy` close expected
                | (actual, expected) <- zip infiniteValues [2 / 3, 1 / 3, 1, 0]
                ]
            checked (Visit.infiniteProbabilityGivenInitialState recurrentVisitChain 2 0)
                `shouldSatisfy` close (2 / 3)
            visitTotalExpectationByState recurrentVisitChain 2
                `shouldBe` Right
                    [ InfiniteExpectation
                    , InfiniteExpectation
                    , InfiniteExpectation
                    , FiniteExpectation 0
                    ]
            Visit.totalExpectationGivenInitialState recurrentVisitChain 2 3
                `shouldBe` Right (FiniteExpectation 0)
