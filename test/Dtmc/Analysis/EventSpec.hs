module Dtmc.Analysis.EventSpec (
    spec,
) where

import Dtmc.Analysis.Event (
    DiscreteEvent (..),
    includesInfiniteOutcome,
    matches,
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
    describe "matches" $ do
        it "implements every comparison at and around its threshold" $ do
            let values = [2, 3, 4]
            map (matches (EqualTo 3)) values
                `shouldBe` [False, True, False]
            map (matches (LessThan 3)) values
                `shouldBe` [True, False, False]
            map (matches (AtMost 3)) values
                `shouldBe` [True, True, False]
            map (matches (GreaterThan 3)) values
                `shouldBe` [False, False, True]
            map (matches (AtLeast 3)) values
                `shouldBe` [False, True, True]

        it "has the structural zero-threshold boundaries" $ do
            matches (LessThan 0) 0 `shouldBe` False
            matches (AtMost 0) 0 `shouldBe` True
            matches (GreaterThan 0) 0 `shouldBe` False
            matches (AtLeast 0) 0 `shouldBe` True

        prop "AtMost n is LessThan (n + 1)" $ \rawThreshold rawValue ->
            let threshold = asNatural rawThreshold
                value = asNatural rawValue
             in matches (AtMost threshold) value
                    == matches (LessThan (threshold + 1)) value

        prop "AtLeast (n + 1) is GreaterThan n" $ \rawThreshold rawValue ->
            let threshold = asNatural rawThreshold
                value = asNatural rawValue
             in matches (AtLeast (threshold + 1)) value
                    == matches (GreaterThan threshold) value

        prop "upper and lower complements partition every finite value" $
            \rawThreshold rawValue ->
                let threshold = asNatural rawThreshold
                    value = asNatural rawValue
                 in property $
                        and
                            [ matches (GreaterThan threshold) value
                                /= matches (AtMost threshold) value
                            , matches (AtLeast threshold) value
                                /= matches (LessThan threshold) value
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
