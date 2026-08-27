{-# LANGUAGE DeriveGeneric #-}

module Dtmc.FacadeSpec (
    spec,
) where

import Dtmc
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
    checked
        ( stationaryDistribution
            ( case witnessIrreducible weatherTransition of
                Nothing -> error "weather transition is not irreducible"
                Just witness -> witness
            )
        )

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
    describe "Dtmc facade" $ do
        it "computes a stationary distribution" $ do
            abs (probabilityAt weatherStationary Dry - 0.8) < 1e-12
                `shouldBe` True
            abs (probabilityAt weatherStationary Wet - 0.2) < 1e-12
                `shouldBe` True

        it "matches the apple-to-mango transition closed form" $
            mapM_
                ( \n ->
                    abs
                        ( transitionProbabilityN
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
                        ( transitionProbabilityN
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
                (checked (hittingProbability cafeTransition [Leave] Thinking) - 1)
                < 1e-12
                `shouldBe` True
            abs
                ( checked
                    (hittingProbability cafeTransition [Drink] Thinking)
                    - 4 / 43
                )
                < 1e-12
                `shouldBe` True
            abs
                ( checked
                    ( hittingBeforeProbability
                        cafeTransition
                        [PlainWaffle, ChocolateWaffle]
                        [Drink, Leave]
                        Thinking
                    )
                    - 29 / 43
                )
                < 1e-12
                `shouldBe` True

        it "exposes finite-horizon visit-count analysis" $
            visitCountExpectationBefore 1 cafeInitial cafeTransition (== Thinking)
                `shouldBe` 1

        it "exposes infinite-horizon total visit-count analysis" $ do
            checked (visitCountProbability weatherTransition Dry InfiniteVisits Wet)
                `shouldBe` 1
            checked (visitCountProbability weatherTransition Dry (FiniteVisits 1) Wet)
                `shouldBe` 0
            checked (visitCountExpectation weatherTransition Dry Wet)
                `shouldBe` InfiniteMeanCount

        it "exposes conditional-probability errors" $
            conditionalProbability
                cafeInitial
                cafeTransition
                []
                [At 0 Leave]
                `shouldBe` Left ZeroProbabilityCondition

        it "exposes joint probabilities of timed observations" $
            jointProbability cafeInitial cafeTransition [At 0 Thinking]
                `shouldBe` 1
