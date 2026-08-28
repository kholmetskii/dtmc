module Dtmc.Analysis.EventSpec (
    spec,
) where

import Dtmc.Analysis.Event (
    DiscreteEvent (..),
    includesInfiniteOutcome,
    matchesDiscreteEvent,
 )
import Numeric.Natural (
    Natural,
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
    NonNegative (..),
    property,
 )

asNatural :: NonNegative Integer -> Natural
asNatural (NonNegative value) = fromInteger value

spec :: Spec
spec = do
    describe "matchesDiscreteEvent" $ do
        it "implements every comparison at and around its threshold" $ do
            let values = [2, 3, 4]
            map (matchesDiscreteEvent (EqualTo 3)) values
                `shouldBe` [False, True, False]
            map (matchesDiscreteEvent (LessThan 3)) values
                `shouldBe` [True, False, False]
            map (matchesDiscreteEvent (AtMost 3)) values
                `shouldBe` [True, True, False]
            map (matchesDiscreteEvent (GreaterThan 3)) values
                `shouldBe` [False, False, True]
            map (matchesDiscreteEvent (AtLeast 3)) values
                `shouldBe` [False, True, True]

        it "has the structural zero-threshold boundaries" $ do
            matchesDiscreteEvent (LessThan 0) 0 `shouldBe` False
            matchesDiscreteEvent (AtMost 0) 0 `shouldBe` True
            matchesDiscreteEvent (GreaterThan 0) 0 `shouldBe` False
            matchesDiscreteEvent (AtLeast 0) 0 `shouldBe` True

        prop "AtMost n is LessThan (n + 1)" $ \rawThreshold rawValue ->
            let threshold = asNatural rawThreshold
                value = asNatural rawValue
             in matchesDiscreteEvent (AtMost threshold) value
                    == matchesDiscreteEvent (LessThan (threshold + 1)) value

        prop "AtLeast (n + 1) is GreaterThan n" $ \rawThreshold rawValue ->
            let threshold = asNatural rawThreshold
                value = asNatural rawValue
             in matchesDiscreteEvent (AtLeast (threshold + 1)) value
                    == matchesDiscreteEvent (GreaterThan threshold) value

        prop "upper and lower complements partition every finite value" $
            \rawThreshold rawValue ->
                let threshold = asNatural rawThreshold
                    value = asNatural rawValue
                 in property $
                        and
                            [ matchesDiscreteEvent (GreaterThan threshold) value
                                /= matchesDiscreteEvent (AtMost threshold) value
                            , matchesDiscreteEvent (AtLeast threshold) value
                                /= matchesDiscreteEvent (LessThan threshold) value
                            ]

    describe "includesInfiniteOutcome" $ do
        it "includes infinity exactly in upper-tail events" $ do
            map
                includesInfiniteOutcome
                [ EqualTo 3
                , LessThan 3
                , AtMost 3
                , GreaterThan 3
                , AtLeast 3
                ]
                `shouldBe` [False, False, False, True, True]

        prop "is independent of the finite threshold" $ \rawThreshold ->
            let threshold = asNatural rawThreshold
             in property $
                    and
                        [ not (includesInfiniteOutcome (EqualTo threshold))
                        , not (includesInfiniteOutcome (LessThan threshold))
                        , not (includesInfiniteOutcome (AtMost threshold))
                        , includesInfiniteOutcome (GreaterThan threshold)
                        , includesInfiniteOutcome (AtLeast threshold)
                        ]
