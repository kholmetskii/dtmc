{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.StationarySpec (
    spec,
) where

import Data.Finite (
    Finite,
 )
import Dtmc.Analysis.Expectation (
    Expectation (..),
 )
import Dtmc.Analysis.ReturnTime qualified as Return
import Dtmc.Analysis.Stationary (
    stationaryDistributions,
 )
import Dtmc.Distribution (
    probabilityAt,
 )
import Dtmc.Distribution.Vector qualified as Vector
import Dtmc.Dynamics (
    evolveVector,
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.TestSupport (
    approxDistributionEq,
    approxEq,
    chunksOf,
    genTransitionRows,
    testTolerance,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    fromRows,
 )
import GHC.Generics (
    Generic,
 )
import Test.Hspec (
    Spec,
    describe,
    expectationFailure,
    it,
    shouldBe,
    shouldSatisfy,
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

onlyStationary ::
    (FiniteState state) =>
    TransitionMatrix state ->
    Vector.DistributionVector state
onlyStationary matrix =
    case checked (stationaryDistributions matrix) of
        [(_, distribution)] -> distribution
        _ -> error "test matrix does not have a unique stationary distribution"

twoState :: TransitionMatrix (Finite 2)
twoState =
    checked
        ( fromRows
            (chunksOf 2 [0.9, 0.1, 0.4, 0.6])
        )

singleton :: TransitionMatrix (Finite 1)
singleton =
    checked
        ( fromRows
            (chunksOf 1 [1])
        )

threeCycle :: TransitionMatrix (Finite 3)
threeCycle =
    checked
        ( fromRows
            (chunksOf 3 [0, 1, 0, 0, 0, 1, 1, 0, 0])
        )

namedTwoState :: TransitionMatrix Weather
namedTwoState =
    checked
        ( fromRows
            (chunksOf 2 [0.9, 0.1, 0.4, 0.6])
        )

genPositiveTransitionMatrix :: Gen [[Double]]
genPositiveTransitionMatrix = vectorOf 3 positiveSimplex
  where
    positiveSimplex = do
        weights <- vectorOf 3 (choose (1, 1000 :: Double))
        let total = sum weights
        pure (map (/ total) weights)

stationaryLawsHold :: [[Double]] -> Property
stationaryLawsHold raw =
    case fromRows @(Finite 3) raw of
        Left err -> counterexample (show err) (property False)
        Right matrix ->
            case stationaryDistributions matrix of
                Left err -> counterexample (show err) (property False)
                Right [(_, distribution)] ->
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
                                (approxEq testTolerance (sum (Vector.toList distribution)) 1)
                        ]
                Right _ -> counterexample "positive matrix was not uniquely stationary" (property False)

spec :: Spec
spec = do
    describe "stationaryDistributions" $ do
        it "returns the point mass for a singleton chain" $
            Vector.toList (onlyStationary singleton)
                `shouldBe` [1]

        it "matches the closed form for a two-state chain" $
            and
                ( zipWith
                    (approxEq testTolerance)
                    (Vector.toList (onlyStationary twoState))
                    [0.8, 0.2]
                )
                `shouldBe` True

        it "is uniform for a periodic three-cycle" $
            and
                [ approxEq testTolerance actual (1 / 3)
                | actual <- Vector.toList (onlyStationary threeCycle)
                ]
                `shouldBe` True

        it "preserves named-state coordinates" $ do
            let distribution =
                    onlyStationary namedTwoState
            approxEq testTolerance (probabilityAt distribution Dry) 0.8
                `shouldBe` True
            approxEq testTolerance (probabilityAt distribution Wet) 0.2
                `shouldBe` True

        prop "satisfies the balance and normalization equations" $
            forAll genPositiveTransitionMatrix stationaryLawsHold

        it "solves a symmetric nearly uncoupled chain exactly" $ do
            -- The balance system is hopelessly ill conditioned here, but GTH
            -- never forms it: the exit mass is accumulated rather than taken
            -- as 1 - P(k,k), so the answer comes out bit-exact.
            let epsilon = 1e-14
                matrix =
                    checked
                        ( fromRows @(Finite 2)
                            ( chunksOf
                                2
                                [ 1 - epsilon
                                , epsilon
                                , epsilon
                                , 1 - epsilon
                                ]
                            )
                        )
            Vector.toList (onlyStationary matrix)
                `shouldBe` [0.5, 0.5]

        it "solves an asymmetric nearly uncoupled chain" $ do
            -- For [[1-a, a], [b, 1-b]] the stationary law is
            -- (b, a) / (a + b), here (3/4, 1/4) at a scale where forming
            -- transpose(P) - I would destroy every significant digit.
            let leaving = 1e-14
                returning = 3e-14
                matrix =
                    checked
                        ( fromRows @(Finite 2)
                            ( chunksOf
                                2
                                [ 1 - leaving
                                , leaving
                                , returning
                                , 1 - returning
                                ]
                            )
                        )
            Vector.toList (onlyStationary matrix)
                `shouldSatisfy` allCloseTo [0.75, 0.25]

        it "normalises extreme finite GTH weights without overflow" $ do
            let epsilon = 5e-309
                matrix =
                    checked
                        ( fromRows @(Finite 3)
                            ( chunksOf
                                3
                                [ 0
                                , 0.5
                                , 0.5
                                , epsilon
                                , 0
                                , 1
                                , epsilon
                                , 1
                                , 0
                                ]
                            )
                        )
                weights = Vector.toList (onlyStationary matrix)
            weights `shouldSatisfy` all isFinite
            sum weights `shouldSatisfy` approxEq testTolerance 1
            weights `shouldSatisfy` allCloseTo [0, 0.5, 0.5]
            case weights of
                first : _ -> first `shouldSatisfy` (> 0)
                [] -> expectationFailure "expected three stationary weights"

    describe "multiple recurrent classes" $ do
        it "returns one distribution per recurrent class, by least member" $
            fmap (map fst) (stationaryDistributions twoClosedClasses)
                `shouldBe` Right [[0], [1, 2]]

        it "matches the closed form of the notes" $
            case stationaryDistributions twoClosedClasses of
                Right [(_, onFirst), (_, onSecond)] -> do
                    Vector.toList onFirst `shouldSatisfy` allCloseTo [1, 0, 0]
                    Vector.toList onSecond `shouldSatisfy` allCloseTo [0, 5 / 11, 6 / 11]
                other -> expectationFailure ("unexpected result: " ++ show other)

        it "puts exact zero on a transient state" $
            case stationaryDistributions withTransient of
                Right [(members, only)] -> do
                    members `shouldBe` [1, 2]
                    take 1 (Vector.toList only) `shouldBe` [0]
                    Vector.toList only `shouldSatisfy` allCloseTo [0, 5 / 11, 6 / 11]
                other -> expectationFailure ("unexpected result: " ++ show other)

        it "returns one distribution for an irreducible chain" $
            case stationaryDistributions twoState of
                Right [(_, only)] ->
                    Vector.toList only `shouldSatisfy` allCloseTo [0.8, 0.2]
                other -> expectationFailure ("unexpected result: " ++ show other)

        it "inverts the mean return time" $
            -- pi_i m_i = 1 for state 1 of the recurrent class {1, 2}
            case stationaryDistributions twoClosedClasses of
                Right [_, (_, onSecond)] ->
                    Return.expectationGivenInitialState twoClosedClasses 1
                        `shouldSatisfy` inverts (Vector.toList onSecond !! 1)
                other -> expectationFailure ("unexpected result: " ++ show other)

        prop "every returned distribution is stationary and normalised" $
            forAll (genTransitionRows 3) $ \raw ->
                case fromRows @(Finite 3) raw of
                    Left err -> counterexample (show err) (property False)
                    Right matrix ->
                        case stationaryDistributions matrix of
                            -- A refused solve is a documented outcome.
                            Left _ -> property True
                            Right results ->
                                conjoin
                                    [ conjoin
                                        [ counterexample "pi P /= pi" $
                                            property
                                                ( approxDistributionEq
                                                    testTolerance
                                                    (evolveVector d matrix)
                                                    d
                                                )
                                        , counterexample "sum pi /= 1" $
                                            property
                                                (approxEq testTolerance (sum (Vector.toList d)) 1)
                                        ]
                                    | (_, d) <- results
                                    ]

-- Section 4.1: two closed classes, hence infinitely many stationary
-- distributions for the chain as a whole.
twoClosedClasses :: TransitionMatrix (Finite 3)
twoClosedClasses =
    checked
        ( fromRows
            (chunksOf 3 [1, 0, 0, 0, 0.4, 0.6, 0, 0.5, 0.5])
        )

-- State 0 is transient; {1, 2} is the only recurrent class.
withTransient :: TransitionMatrix (Finite 3)
withTransient =
    checked
        ( fromRows
            (chunksOf 3 [0, 0.5, 0.5, 0, 0.4, 0.6, 0, 0.5, 0.5])
        )

allCloseTo :: [Double] -> [Double] -> Bool
allCloseTo expected actual =
    length expected == length actual
        && and (zipWith (approxEq testTolerance) expected actual)

isFinite :: Double -> Bool
isFinite value = not (isNaN value || isInfinite value)

inverts :: Double -> Either error Expectation -> Bool
inverts probability (Right (FiniteExpectation mean)) =
    approxEq testTolerance (probability * mean) 1
inverts _ _ = False
