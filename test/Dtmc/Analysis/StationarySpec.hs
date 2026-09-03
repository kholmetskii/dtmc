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
    genTransitionMatrix,
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

onlyStationary :: (FiniteState state) => TransitionMatrix state -> DistributionVector state
onlyStationary matrix =
    case checked (stationaryDistributions matrix) of
        [(_, distribution)] -> distribution
        _ -> error "test matrix does not have a unique stationary distribution"

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
                                (approxEq testTolerance (sum (entries distribution)) 1)
                        ]
                Right _ -> counterexample "positive matrix was not uniquely stationary" (property False)

spec :: Spec
spec = do
    describe "stationaryDistributions" $ do
        it "returns the point mass for a singleton chain" $
            entries (onlyStationary singleton)
                `shouldBe` [1]

        it "matches the closed form for a two-state chain" $
            and
                ( zipWith
                    (approxEq testTolerance)
                    (entries (onlyStationary twoState))
                    [0.8, 0.2]
                )
                `shouldBe` True

        it "is uniform for a periodic three-cycle" $
            and
                [ approxEq testTolerance actual (1 / 3)
                | actual <- entries (onlyStationary threeCycle)
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
            entries (onlyStationary matrix)
                `shouldBe` [0.5, 0.5]

        it "solves an asymmetric nearly uncoupled chain" $ do
            -- For [[1-a, a], [b, 1-b]] the stationary law is
            -- (b, a) / (a + b), here (3/4, 1/4) at a scale where forming
            -- transpose(P) - I would destroy every significant digit.
            let leaving = 1e-14
                returning = 3e-14
                matrix =
                    checked
                        ( mkTransitionMatrix @(Finite 2)
                            ( S.matrix
                                [ 1 - leaving
                                , leaving
                                , returning
                                , 1 - returning
                                ] ::
                                S.Sq 2
                            )
                        )
            entries (onlyStationary matrix)
                `shouldSatisfy` allCloseTo [0.75, 0.25]

    describe "multiple recurrent classes" $ do
        it "returns one distribution per recurrent class, by least member" $
            fmap (map fst) (stationaryDistributions twoClosedClasses)
                `shouldBe` Right [[0], [1, 2]]

        it "matches the closed form of the notes" $
            case stationaryDistributions twoClosedClasses of
                Right [(_, onFirst), (_, onSecond)] -> do
                    entries onFirst `shouldSatisfy` allCloseTo [1, 0, 0]
                    entries onSecond `shouldSatisfy` allCloseTo [0, 5 / 11, 6 / 11]
                other -> expectationFailure ("unexpected result: " ++ show other)

        it "puts exact zero on a transient state" $
            case stationaryDistributions withTransient of
                Right [(members, only)] -> do
                    members `shouldBe` [1, 2]
                    take 1 (entries only) `shouldBe` [0]
                    entries only `shouldSatisfy` allCloseTo [0, 5 / 11, 6 / 11]
                other -> expectationFailure ("unexpected result: " ++ show other)

        it "returns one distribution for an irreducible chain" $
            case stationaryDistributions twoState of
                Right [(_, only)] ->
                    entries only `shouldSatisfy` allCloseTo [0.8, 0.2]
                other -> expectationFailure ("unexpected result: " ++ show other)

        it "inverts the mean return time" $
            -- pi_i m_i = 1 for state 1 of the recurrent class {1, 2}
            case stationaryDistributions twoClosedClasses of
                Right [_, (_, onSecond)] ->
                    Return.expectationGivenInitialState twoClosedClasses 1
                        `shouldSatisfy` inverts (entries onSecond !! 1)
                other -> expectationFailure ("unexpected result: " ++ show other)

        prop "every returned distribution is stationary and normalised" $
            forAll (genTransitionMatrix @3) $ \raw ->
                case mkTransitionMatrix @(Finite 3) raw of
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
                                                (approxEq testTolerance (sum (entries d)) 1)
                                        ]
                                    | (_, d) <- results
                                    ]

-- Section 4.1: two closed classes, hence infinitely many stationary
-- distributions for the chain as a whole.
twoClosedClasses :: TransitionMatrix (Finite 3)
twoClosedClasses =
    checked
        ( mkTransitionMatrix
            (S.matrix [1, 0, 0, 0, 0.4, 0.6, 0, 0.5, 0.5] :: S.Sq 3)
        )

-- State 0 is transient; {1, 2} is the only recurrent class.
withTransient :: TransitionMatrix (Finite 3)
withTransient =
    checked
        ( mkTransitionMatrix
            (S.matrix [0, 0.5, 0.5, 0, 0.4, 0.6, 0, 0.5, 0.5] :: S.Sq 3)
        )

allCloseTo :: [Double] -> [Double] -> Bool
allCloseTo expected actual =
    length expected == length actual
        && and (zipWith (approxEq testTolerance) expected actual)

inverts :: Double -> Either error Expectation -> Bool
inverts probability (Right (FiniteExpectation mean)) =
    approxEq testTolerance (probability * mean) 1
inverts _ _ = False
