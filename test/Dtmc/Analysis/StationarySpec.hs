{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.StationarySpec (
    spec,
) where

import Data.Finite (
    Finite,
 )
import Dtmc.Analysis.Classification (
    Irreducible,
    witnessIrreducible,
 )
import Dtmc.Analysis.Stationary (
    LinearSystemError (IllConditionedSystem),
    stationaryDistribution,
 )
import Dtmc.Distribution (
    probabilityAt,
 )
import Dtmc.Distribution.Vector (
    DistributionVector,
    unDistributionVector,
 )
import Dtmc.Dynamics (
    evolveVector,
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.TestSupport (
    approxDistributionEq,
    approxEq,
    testTolerance,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    mkTransitionMatrix,
 )
import GHC.Generics (
    Generic,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    expectationFailure,
    it,
    shouldBe,
 )
import Test.Hspec.QuickCheck (
    prop,
 )
import Test.QuickCheck (
    Gen,
    Property,
    choose,
    conjoin,
    counterexample,
    forAll,
    property,
    vectorOf,
 )

data Weather = Dry | Wet
    deriving (Eq, Ord, Show, Generic)

instance FiniteState Weather

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

certified :: TransitionMatrix state -> Irreducible state
certified matrix =
    case witnessIrreducible matrix of
        Nothing -> error "test matrix is not irreducible"
        Just witness -> witness

twoState :: TransitionMatrix (Finite 2)
twoState =
    checked
        ( mkTransitionMatrix
            (S.matrix [0.9, 0.1, 0.4, 0.6] :: S.Sq 2)
        )

singleton :: TransitionMatrix (Finite 1)
singleton =
    checked
        ( mkTransitionMatrix
            (S.matrix [1] :: S.Sq 1)
        )

threeCycle :: TransitionMatrix (Finite 3)
threeCycle =
    checked
        ( mkTransitionMatrix
            (S.matrix [0, 1, 0, 0, 0, 1, 1, 0, 0] :: S.Sq 3)
        )

namedTwoState :: TransitionMatrix Weather
namedTwoState =
    checked
        ( mkTransitionMatrix
            (S.matrix [0.9, 0.1, 0.4, 0.6] :: S.Sq 2)
        )

entries :: (FiniteState state) => DistributionVector state -> [Double]
entries = LA.toList . S.extract . unDistributionVector

genPositiveTransitionMatrix :: Gen (S.Sq 3)
genPositiveTransitionMatrix = do
    rows <- vectorOf 3 positiveSimplex
    pure (S.matrix (concat rows))
  where
    positiveSimplex = do
        weights <- vectorOf 3 (choose (1, 1000 :: Double))
        let total = sum weights
        pure (map (/ total) weights)

stationaryLawsHold :: S.Sq 3 -> Property
stationaryLawsHold raw =
    case mkTransitionMatrix @(Finite 3) raw of
        Left err -> counterexample (show err) (property False)
        Right matrix ->
            case witnessIrreducible matrix of
                Nothing -> counterexample "positive matrix was reducible" (property False)
                Just witness ->
                    case stationaryDistribution witness of
                        Left err -> counterexample (show err) (property False)
                        Right distribution ->
                            conjoin
                                [ counterexample "pi P /= pi" $
                                    property
                                        ( approxDistributionEq
                                            testTolerance
                                            (evolveVector distribution matrix)
                                            distribution
                                        )
                                , counterexample "sum pi /= 1" $
                                    property
                                        (approxEq testTolerance (sum (entries distribution)) 1)
                                ]

spec :: Spec
spec =
    describe "stationaryDistribution" $ do
        it "returns the point mass for a singleton chain" $
            entries (checked (stationaryDistribution (certified singleton)))
                `shouldBe` [1]

        it "matches the closed form for a two-state chain" $
            and
                ( zipWith
                    (approxEq testTolerance)
                    (entries (checked (stationaryDistribution (certified twoState))))
                    [0.8, 0.2]
                )
                `shouldBe` True

        it "is uniform for a periodic three-cycle" $
            and
                [ approxEq testTolerance actual (1 / 3)
                | actual <- entries (checked (stationaryDistribution (certified threeCycle)))
                ]
                `shouldBe` True

        it "preserves named-state coordinates" $ do
            let distribution =
                    checked (stationaryDistribution (certified namedTwoState))
            approxEq testTolerance (probabilityAt distribution Dry) 0.8
                `shouldBe` True
            approxEq testTolerance (probabilityAt distribution Wet) 0.2
                `shouldBe` True

        prop "satisfies the balance and normalization equations" $
            forAll genPositiveTransitionMatrix stationaryLawsHold

        it "reports an ill-conditioned balance system explicitly" $ do
            let epsilon = 1e-14
                matrix =
                    checked
                        ( mkTransitionMatrix @(Finite 2)
                            ( S.matrix
                                [ 1 - epsilon
                                , epsilon
                                , epsilon
                                , 1 - epsilon
                                ] ::
                                S.Sq 2
                            )
                        )
            case stationaryDistribution (certified matrix) of
                Left IllConditionedSystem{} -> pure ()
                result ->
                    expectationFailure
                        ("expected IllConditionedSystem, got " ++ show result)
