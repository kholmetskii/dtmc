{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE EmptyDataDecls #-}
{-# LANGUAGE EmptyDataDeriving #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.StateSpec (
    spec,
) where

import Data.Finite (
    Finite,
    finites,
    getFinite,
 )
import Data.List (
    sort,
 )
import Data.Proxy (
    Proxy (Proxy),
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
    finiteStates,
    stateAt,
    stateIndex,
 )
import GHC.Generics (
    Generic,
 )
import GHC.TypeNats (
    natVal,
 )
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
 )

data Empty
    deriving (Eq, Ord, Show, Generic)

data One = One
    deriving (Eq, Ord, Show, Generic)

data Three = A | B | C
    deriving (Eq, Ord, Show, Generic, FiniteState)

instance FiniteState Empty

instance FiniteState One

spec :: Spec
spec = do
    describe "generic FiniteState" $ do
        it "can be derived directly in the state declaration" $
            finiteStates @Three `shouldBe` [A, B, C]

        it "derives cardinalities for empty, singleton, and sum types" $ do
            natVal (Proxy @(Cardinality Empty)) `shouldBe` 0
            natVal (Proxy @(Cardinality One)) `shouldBe` 1
            natVal (Proxy @(Cardinality Three)) `shouldBe` 3

        it "enumerates states in constructor order" $
            finiteStates @Three `shouldBe` [A, B, C]

        it "uses the same order as a stock Ord instance" $
            finiteStates @Three `shouldBe` sort (finiteStates @Three)

        it "round-trips every state through its finite index" $
            map (stateAt . stateIndex) (finiteStates @Three)
                `shouldBe` finiteStates @Three

        it "round-trips every finite index through its state" $
            map (stateIndex . stateAt @Three) finites `shouldBe` finites

        it "assigns consecutive zero-based indices" $
            map (getFinite . stateIndex) (finiteStates @Three)
                `shouldBe` [0, 1, 2]

        it "supports empty and singleton state types" $ do
            finiteStates @Empty `shouldBe` []
            finiteStates @One `shouldBe` [One]
            stateAt (stateIndex One) `shouldBe` One

        describe "Empty laws" (finiteStateLaws @Empty)
        describe "One laws" (finiteStateLaws @One)
        describe "Three laws" (finiteStateLaws @Three)

    describe "Finite identity instance" $ do
        it "preserves enumeration and both conversions" $ do
            finiteStates @(Finite 3) `shouldBe` finites
            map stateIndex (finites @3) `shouldBe` finites
            map stateAt (finites @3) `shouldBe` finites

        it "supports Finite 0" $
            finiteStates @(Finite 0) `shouldBe` []

        describe "Finite 0 laws" (finiteStateLaws @(Finite 0))
        describe "Finite 3 laws" (finiteStateLaws @(Finite 3))

    describe "base instances" $ do
        it "use their standard constructor order" $ do
            finiteStates @() `shouldBe` [()]
            finiteStates @Bool `shouldBe` [False, True]
            finiteStates @Ordering `shouldBe` [LT, EQ, GT]

        describe "() laws" (finiteStateLaws @())
        describe "Bool laws" (finiteStateLaws @Bool)
        describe "Ordering laws" (finiteStateLaws @Ordering)

finiteStateLaws ::
    forall state.
    (FiniteState state, Show state) =>
    Spec
finiteStateLaws = do
    it "enumerates every finite index in canonical order" $
        finiteStates @state
            `shouldBe` map (stateAt @state) (finites @(Cardinality state))

    it "round-trips every state through its index" $
        map (stateAt . stateIndex) (finiteStates @state)
            `shouldBe` finiteStates @state

    it "round-trips every index through its state" $
        map (stateIndex . stateAt @state) (finites @(Cardinality state))
            `shouldBe` finites @(Cardinality state)

    it "enumerates states in strictly ascending order" $
        and (zipWith (<) states (drop 1 states)) `shouldBe` True
  where
    states = finiteStates @state
