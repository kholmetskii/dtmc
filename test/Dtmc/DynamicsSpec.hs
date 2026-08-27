{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.DynamicsSpec (
    spec,
) where

import Data.Finite (
    Finite,
    finites,
 )
import Data.Proxy (
    Proxy (..),
 )
import Dtmc.Analysis.FiniteTime (
    transitionProbability,
 )
import Dtmc.Distribution (
    distributionWeights,
    probabilityAt,
 )
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.Distribution.Vector (
    DistributionVector,
    mkDistributionVector,
    unDistributionVector,
 )
import Dtmc.Dynamics (
    evolveN,
    evolveVector,
    evolveVectorN,
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.TestSupport (
    approxDistributionEq,
    approxEq,
    genSimplexPoint,
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.Transition.Kernel qualified as Kernel
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    mkTransitionMatrix,
 )
import GHC.Generics (
    Generic,
 )
import GHC.TypeNats (
    KnownNat,
    natVal,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
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
    Gen,
    choose,
    counterexample,
    forAll,
    property,
 )

data NamedPosition = LowerPosition | UpperPosition
    deriving (Eq, Ord, Show, Generic)

instance FiniteState NamedPosition

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

kernelChain :: Kernel.TransitionKernel (Finite 3)
kernelChain =
    Kernel.transitionKernel $ \source ->
        checked $
            DistributionMap.mkDistributionMap
                [ (destination, transitionProbability finiteChain source destination)
                | destination <- finites
                ]

mapInitial :: DistributionMap.DistributionMap (Finite 3)
mapInitial =
    checked $
        DistributionMap.mkDistributionMap
            [(state, probabilityAt finiteInitial state) | state <- finites]

simpleRandomWalk :: Kernel.TransitionKernel Integer
simpleRandomWalk =
    Kernel.transitionKernel $ \state ->
        checked
            (DistributionMap.mkDistributionMap [(state - 1, 0.5), (state + 1, 0.5)])

closeTo :: Double -> Double -> Bool
closeTo = approxEq testTolerance

genDistribution :: forall n. (KnownNat n) => Gen (S.R n)
genDistribution = do
    entries <- genSimplexPoint (fromIntegral (natVal (Proxy @n)))
    pure (S.vector entries)

twoState :: TransitionMatrix (Finite 2)
twoState =
    either (error . show) id $
        mkTransitionMatrix
            (S.matrix [0.9, 0.1, 0.4, 0.6] :: S.Sq 2)

namedInitial :: DistributionVector NamedPosition
namedInitial =
    either (error . show) id $
        mkDistributionVector @NamedPosition
            (S.vector [1, 0] :: S.R 2)

namedTwoState :: TransitionMatrix NamedPosition
namedTwoState =
    either (error . show) id $
        mkTransitionMatrix @NamedPosition
            (S.matrix [0.9, 0.1, 0.4, 0.6] :: S.Sq 2)

spec :: Spec
spec = do
    describe "evolveVector" $ do
        prop "keeps the distribution on the simplex" $
            forAll ((,) <$> genDistribution @3 <*> genTransitionMatrix @3) $
                \(vector, matrix) ->
                    case (mkDistributionVector @(Finite 3) vector, mkTransitionMatrix matrix) of
                        (Right mu, Right p) ->
                            case mkDistributionVector @(Finite 3) (unDistributionVector (evolveVector mu p)) of
                                Right _ ->
                                    property True
                                Left err ->
                                    counterexample
                                        ("evolved distribution left the simplex: " <> show err)
                                        False
                        result ->
                            counterexample
                                ("generated input was rejected: " <> show result)
                                False

        it "matches a hand-computed two-state step" $ do
            let mu =
                    either (error . show) id $
                        mkDistributionVector @(Finite 2)
                            (S.vector [1, 0] :: S.R 2)

            LA.toList (S.extract (unDistributionVector (evolveVector mu twoState)))
                `shouldBe` [0.9, 0.1]

        it "preserves named states while evolving the dense vector" $ do
            probabilityAt (evolveVector namedInitial namedTwoState) LowerPosition
                `shouldBe` 0.9
            probabilityAt (evolveVector namedInitial namedTwoState) UpperPosition
                `shouldBe` 0.1

    describe "evolveVectorN" $ do
        it "leaves a distribution unchanged after zero steps" $ do
            let mu =
                    either (error . show) id $
                        mkDistributionVector @(Finite 2)
                            (S.vector [0.25, 0.75] :: S.R 2)

            approxDistributionEq
                1e-12
                (evolveVectorN 0 mu twoState)
                mu
                `shouldBe` True

        prop "agrees with iterating evolveVector"
            $ forAll
                ( (,,)
                    <$> choose (0, 6 :: Int)
                    <*> genDistribution @3
                    <*> genTransitionMatrix @3
                )
            $ \(k, vector, matrix) ->
                case (mkDistributionVector @(Finite 3) vector, mkTransitionMatrix matrix) of
                    (Right mu, Right p) ->
                        let iterated =
                                iterate (`evolveVector` p) mu !! k
                         in property $
                                approxDistributionEq
                                    1e-9
                                    (evolveVectorN (fromIntegral k) mu p)
                                    iterated
                    result ->
                        counterexample
                            ("generated input was rejected: " <> show result)
                            False

        prop "composes m steps then n steps"
            $ forAll
                ( (,,,)
                    <$> choose (0, 4 :: Int)
                    <*> choose (0, 4 :: Int)
                    <*> genDistribution @3
                    <*> genTransitionMatrix @3
                )
            $ \(m, n, vector, matrix) ->
                case (mkDistributionVector @(Finite 3) vector, mkTransitionMatrix matrix) of
                    (Right mu, Right p) ->
                        property $
                            approxDistributionEq
                                1e-9
                                (evolveVectorN (fromIntegral (m + n)) mu p)
                                ( evolveVectorN
                                    (fromIntegral n)
                                    (evolveVectorN (fromIntegral m) mu p)
                                    p
                                )
                    result ->
                        counterexample
                            ("generated input was rejected: " <> show result)
                            False

    describe "evolve/evolveN" $ do
        it "evolves an infinite-state random walk without enumerating its state space" $
            distributionWeights
                (evolveN 2 (DistributionMap.pointMass 0) simpleRandomWalk)
                `shouldBe` [(-2, 0.25), (0, 0.5), (2, 0.25)]

        it "agrees across equivalent matrix and kernel representations" $
            sequence_
                [ probabilityAt (evolveN time mapInitial kernelChain) state
                    `shouldSatisfy` closeTo
                        (probabilityAt (evolveVectorN time finiteInitial finiteChain) state)
                | time <- [0 .. 4]
                , state <- finites :: [Finite 3]
                ]
