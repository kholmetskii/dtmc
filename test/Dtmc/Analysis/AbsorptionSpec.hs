{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.AbsorptionSpec (
    spec,
) where

import Data.Finite (
    Finite,
 )
import Dtmc.Analysis.Absorption qualified as Absorption
import Dtmc.Analysis.Classification (
    recurrentStates,
    transientStates,
 )
import Dtmc.Analysis.Expectation qualified as E
import Dtmc.Analysis.HittingTime qualified as Hitting
import Dtmc.State (
    FiniteState,
 )
import Dtmc.TestSupport (
    approxEq,
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    TransitionMatrixError,
    mkTransitionMatrix,
 )
import GHC.Generics (
    Generic,
 )
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec
import Test.Hspec.QuickCheck (
    prop,
 )
import Test.QuickCheck

-- The chain of the "transient class {A,B}, recurrent class {C,D}" example in
-- section 3.7 of the notes, which states G, eta and the hit-before
-- probabilities in closed form.
data Four = A | B | C | D
    deriving (Eq, Ord, Show, Generic, FiniteState)

chain :: TransitionMatrix Four
chain =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 0, 1 / 3, 2 / 3, 0
                , 1 / 2, 0, 1 / 8, 3 / 8
                , 0, 0, 1 / 2, 1 / 2
                , 0, 0, 3 / 4, 1 / 4
                ]
            )

-- Two states in one recurrent class: no transient states at all.
twoCycle :: TransitionMatrix Bool
twoCycle =
    either (error . show) id $
        mkTransitionMatrix (S.matrix [0, 1, 1, 0])

closeTo :: Double -> Double -> Bool
closeTo expected actual = approxEq testTolerance expected actual

finiteCloseTo :: Double -> E.Expectation -> Bool
finiteCloseTo expected (E.FiniteExpectation x) = closeTo expected x
finiteCloseTo _ E.InfiniteExpectation = False

rightCloseTo :: Double -> Either err Double -> Bool
rightCloseTo expected (Right actual) = closeTo expected actual
rightCloseTo _ (Left _) = False

spec :: Spec
spec = do
    describe "canonicalOrder" $
        it "splits the example chain into {A,B} and {C,D}" $
            Absorption.canonicalOrder chain `shouldBe` ([A, B], [C, D])

    describe "fundamentalMatrix" $ do
        it "indexes rows and columns by the transient states" $
            fmap fst (Absorption.fundamentalMatrix chain) `shouldBe` Right [A, B]

        it "matches the closed form of the notes" $
            case Absorption.fundamentalMatrix chain of
                Left err -> expectationFailure ("solve failed: " <> show err)
                Right (_, rows) ->
                    concat rows
                        `shouldSatisfy` ( and
                                            . zipWith closeTo [6 / 5, 2 / 5, 3 / 5, 6 / 5]
                                        )

        it "returns an empty block when no state is transient" $
            Absorption.fundamentalMatrix twoCycle `shouldBe` Right ([], [])

    describe "probability" $ do
        it "reproduces the hit-before probabilities of the notes" $ do
            Absorption.probability chain C A `shouldSatisfy` rightCloseTo (17 / 20)
            Absorption.probability chain D A `shouldSatisfy` rightCloseTo (3 / 20)
            Absorption.probability chain C B `shouldSatisfy` rightCloseTo (11 / 20)
            Absorption.probability chain D B `shouldSatisfy` rightCloseTo (9 / 20)

        it "is exact at a recurrent starting state" $ do
            Absorption.probability chain C C `shouldBe` Right 1
            Absorption.probability chain D C `shouldBe` Right 0

        it "is exactly zero for a transient target" $
            Absorption.probability chain A B `shouldBe` Right 0

        it "differs from ever hitting the same state" $ do
            -- The chain reaches C almost surely, but it enters {C,D} at D
            -- with probability 3/20, so B(A,C) is strictly smaller.
            Hitting.eventualProbability chain [C] A `shouldSatisfy` rightCloseTo 1
            Absorption.probability chain C A `shouldSatisfy` rightCloseTo (17 / 20)

        it "agrees with eventual hitting after summing over the class" $
            case (Absorption.probability chain C A, Absorption.probability chain D A) of
                (Right toC, Right toD) ->
                    Hitting.eventualProbability chain [C, D] A
                        `shouldSatisfy` rightCloseTo (toC + toD)
                other -> expectationFailure ("solve failed: " <> show other)

    describe "expectationByState" $ do
        it "matches the closed form of the notes" $
            case Absorption.expectationByState chain of
                Left err -> expectationFailure ("solve failed: " <> show err)
                Right values ->
                    values
                        `shouldSatisfy` ( and
                                            . zipWith finiteCloseTo [8 / 5, 9 / 5, 0, 0]
                                        )

        it "equals the expected hitting time of the recurrent set" $
            Absorption.expectationByState chain
                `shouldBe` Hitting.expectationByState chain [C, D]

    describe "absorption probabilities" $
        prop "sum to one from every transient state" $
            forAll (genTransitionMatrix @3) $ \m ->
                case mkTransitionMatrix m ::
                    Either TransitionMatrixError (TransitionMatrix (Finite 3)) of
                    Left err ->
                        counterexample ("generated matrix rejected: " <> show err) False
                    Right p ->
                        conjoin
                            [ case traverse (\k -> Absorption.probability p k i) (recurrentStates p) of
                                -- A refused solve is a documented outcome, not a
                                -- violated law.
                                Left _ -> property True
                                Right ps ->
                                    counterexample
                                        ("from " <> show i <> ": " <> show ps)
                                        (approxEq 1e-9 1 (sum ps))
                            | i <- transientStates p
                            ]
