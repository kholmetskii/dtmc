{-# LANGUAGE DeriveGeneric #-}

module Dtmc.IntegrationSpec (
    spec,
) where

import Dtmc.Analysis.Classification (
    absorbingStates,
    classify,
    reachesAny,
 )
import Dtmc.Analysis.Event (
    DiscreteEvent (..),
 )
import Dtmc.Analysis.Expectation (
    Expectation (..),
 )
import Dtmc.Analysis.FiniteTime qualified as FT
import Dtmc.Analysis.HittingTime qualified as Hit
import Dtmc.Analysis.Stationary (
    stationaryDistributions,
 )
import Dtmc.Analysis.VisitCount qualified as Visit
import Dtmc.Distribution (
    Distribution (..),
 )
import Dtmc.Distribution.Vector (
    DistributionVector,
 )
import Dtmc.Distribution.Vector.HMatrix (
    mkDistributionVector,
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
 )
import Dtmc.Transition.Matrix.HMatrix (
    mkTransitionMatrix,
 )
import GHC.Generics (
    Generic,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
 )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

data CafeState
    = Thinking
    | Menu
    | Drink
    | Food
    | PlainWaffle
    | ChocolateWaffle
    | Leave
    deriving (Eq, Ord, Show, Generic)

instance FiniteState CafeState

data FruitState
    = Apple
    | Pear
    | Banana
    | Mango
    | Kiwi
    | Watermelon
    | Grapefruit
    deriving (Eq, Ord, Show, Generic)

instance FiniteState FruitState

data Weather = Dry | Wet
    deriving (Eq, Ord, Show, Generic)

instance FiniteState Weather

weatherTransition :: TransitionMatrix Weather
weatherTransition =
    checked
        ( mkTransitionMatrix
            (S.matrix [0.9, 0.1, 0.4, 0.6] :: S.Sq 2)
        )

weatherStationary :: DistributionVector Weather
weatherStationary =
    case checked (stationaryDistributions weatherTransition) of
        [(_, distribution)] -> distribution
        _ -> error "weather transition does not have a unique stationary distribution"

fruitTransition :: TransitionMatrix FruitState
fruitTransition =
    checked
        ( mkTransitionMatrix
            ( S.matrix
                [ 0
                , 0
                , 1 / 2
                , 1 / 2
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1 / 3
                , 1 / 3
                , 1 / 3
                , 0
                , 0
                , 0
                , 0
                , 0
                , 2 / 3
                , 1 / 3
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0
                ] ::
                S.Sq 7
            )
        )

appleToMangoProbability :: Int -> Double
appleToMangoProbability n =
    5 / 7 - 3 / 14 * ((-(1 / 6)) ^ n)

mangoToPearProbability :: Int -> Double
mangoToPearProbability n =
    3 / 7 - 2 / 21 * ((-(1 / 6)) ^ n)

cafeInitial :: DistributionVector CafeState
cafeInitial =
    checked
        ( mkDistributionVector
            (S.vector [1, 0, 0, 0, 0, 0, 0] :: S.R 7)
        )

cafeTransition :: TransitionMatrix CafeState
cafeTransition =
    checked
        ( mkTransitionMatrix
            ( S.matrix
                [ 0
                , 1 / 5
                , 0
                , 1 / 5
                , 1 / 5
                , 1 / 5
                , 1 / 5
                , 1 / 5
                , 0
                , 2 / 5
                , 0
                , 2 / 5
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1 / 2
                , 0
                , 0
                , 1 / 2
                , 1 / 2
                , 0
                , 0
                , 0
                , 0
                , 1 / 2
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 0
                , 0
                , 0
                , 0
                , 1
                ] ::
                S.Sq 7
            )
        )

spec :: Spec
spec =
    describe "public module integration" $ do
        it "computes a stationary distribution" $ do
            abs (probabilityAt weatherStationary Dry - 0.8) < 1e-12
                `shouldBe` True
            abs (probabilityAt weatherStationary Wet - 0.2) < 1e-12
                `shouldBe` True

        it "matches the apple-to-mango transition closed form" $
            mapM_
                ( \n ->
                    abs
                        ( FT.nStepProbability
                            (3 * n + 1)
                            fruitTransition
                            Apple
                            Mango
                            - appleToMangoProbability (fromIntegral n)
                        )
                        < 1e-12
                        `shouldBe` True
                )
                ([0, 1, 2, 3, 675] :: [Natural])

        it "matches the mango-to-pear transition closed form" $
            mapM_
                ( \n ->
                    abs
                        ( FT.nStepProbability
                            (3 * n + 2)
                            fruitTransition
                            Mango
                            Pear
                            - mangoToPearProbability (fromIntegral n)
                        )
                        < 1e-12
                        `shouldBe` True
                )
                ([0, 1, 2, 3, 4] :: [Natural])

        it "runs the seven-state cafe analysis entirely with named states" $ do
            probabilityAt cafeInitial Thinking `shouldBe` 1
            reachesAny cafeTransition Thinking [Leave] `shouldBe` True
            absorbingStates (classify cafeTransition) `shouldBe` [Leave]
            abs
                (checked (Hit.eventualProbabilityGivenInitialState cafeTransition [Leave] Thinking) - 1)
                < 1e-12
                `shouldBe` True
            abs
                ( checked
                    (Hit.eventualProbabilityGivenInitialState cafeTransition [Drink] Thinking)
                    - 4 / 43
                )
                < 1e-12
                `shouldBe` True
            abs
                ( checked
                    ( Hit.raceProbabilityGivenInitialState
                        cafeTransition
                        [PlainWaffle, ChocolateWaffle]
                        [Drink, Leave]
                        Thinking
                    )
                    - 29 / 43
                )
                < 1e-12
                `shouldBe` True

        it "uses qualified finite-horizon visit-count analysis" $
            Visit.boundedExpectation 1 cafeInitial cafeTransition (== Thinking)
                `shouldBe` 1

        it "uses qualified infinite-horizon total visit-count analysis" $ do
            checked (Visit.infiniteProbabilityGivenInitialState weatherTransition Dry Wet)
                `shouldBe` 1
            checked (Visit.totalProbabilityGivenInitialState (EqualTo 1) weatherTransition Dry Wet)
                `shouldBe` 0
            checked (Visit.totalExpectationGivenInitialState weatherTransition Dry Wet)
                `shouldBe` InfiniteExpectation

        it "uses qualified conditional-probability errors" $
            FT.probabilityGiven
                cafeInitial
                cafeTransition
                []
                [FT.At 0 Leave]
                `shouldBe` Left FT.ZeroProbabilityCondition

        it "uses qualified timed-observation probabilities" $
            FT.probability cafeInitial cafeTransition [FT.At 0 Thinking]
                `shouldBe` 1
