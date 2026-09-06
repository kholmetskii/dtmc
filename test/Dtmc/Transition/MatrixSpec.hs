{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Transition.MatrixSpec (
    spec,
) where

import Data.Finite (
    Finite,
 )
import Dtmc.Distribution (
    distributionWeights,
 )
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.Distribution.Vector qualified as Vector
import Dtmc.Simplex (SimplexError (..))
import Dtmc.State (
    FiniteState,
 )
import Dtmc.TestSupport (
    approxEq,
    approxTransitionMatrixEq,
    bumpSmallestInFirstRow,
    chunksOf,
    genTransitionRows,
    setFirstEntry,
    testTolerance,
 )
import Dtmc.Transition.Kernel qualified as Kernel
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    TransitionMatrixError (..),
    compose,
    fromKernel,
    fromRows,
    identity,
    power,
    rowAt,
    toRows,
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
 )
import Test.Hspec.QuickCheck (
    prop,
 )
import Test.QuickCheck (
    choose,
    conjoin,
    counterexample,
    forAll,
    property,
 )

data NamedPhase = PhaseA | PhaseB | PhaseC
    deriving (Eq, Ord, Show, Generic)

instance FiniteState NamedPhase

cyclicThree :: TransitionMatrix (Finite 3)
cyclicThree =
    either (error . show) id $
        fromRows
            (chunksOf 3 [0, 1, 0, 0, 0, 1, 1, 0, 0])

namedCycle :: TransitionMatrix NamedPhase
namedCycle =
    either (error . show) id $
        fromRows @NamedPhase
            (chunksOf 3 [0, 1, 0, 0, 0, 1, 1, 0, 0])

twoState :: TransitionMatrix (Finite 2)
twoState =
    either (error . show) id $
        fromRows
            (chunksOf 2 [0.9, 0.1, 0.4, 0.6])

twoStateSquared :: TransitionMatrix (Finite 2)
twoStateSquared =
    either (error . show) id $
        fromRows
            (chunksOf 2 [0.85, 0.15, 0.6, 0.4])

spec :: Spec
spec = do
    describe "fromKernel" $ do
        let successor state =
                case state of
                    PhaseA -> PhaseB
                    PhaseB -> PhaseC
                    PhaseC -> PhaseA
            materialised =
                fromKernel
                    (Kernel.fromLaws (DistributionMap.pointMass . successor)) ::
                    TransitionMatrix NamedPhase

        it "materialises a finite deterministic kernel" $
            approxTransitionMatrixEq 0 materialised namedCycle `shouldBe` True

        it "exposes rows without hmatrix types" $
            toRows materialised
                `shouldBe` [[0, 1, 0], [0, 0, 1], [1, 0, 0]]

        it "materialises the empty finite chain" $
            toRows
                ( fromKernel (Kernel.fromLaws DistributionMap.pointMass) ::
                    TransitionMatrix (Finite 0)
                )
                `shouldBe` []

    describe "fromRows" $ do
        prop "stores canonical rows close to the accepted input" $
            forAll (genTransitionRows 3) $ \matrix ->
                case fromRows @(Finite 3) matrix of
                    Right transitionMatrix ->
                        let storedRows = toRows transitionMatrix
                            inputRows = matrix
                            closeRow left right =
                                and
                                    ( zipWith
                                        (approxEq testTolerance)
                                        left
                                        right
                                    )
                         in counterexample ("stored rows: " <> show storedRows) $
                                property
                                    ( all
                                        ( \row ->
                                            all
                                                (\entry -> entry >= 0 && entry <= 1)
                                                row
                                                && approxEq 1e-12 (sum row) 1
                                        )
                                        storedRows
                                        && and (zipWith closeRow storedRows inputRows)
                                    )
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

        it "canonicalises tolerated error independently in each row" $
            case fromRows @(Finite 2)
                (chunksOf 2 [-5e-10, 1 + 5e-10, 0.5, 0.5 - 5e-10]) of
                Right transitionMatrix -> do
                    let rows = toRows transitionMatrix
                    case rows of
                        firstRow : secondRow : _ -> do
                            firstRow `shouldBe` [0, 1]
                            approxEq 1e-12 (sum secondRow) 1 `shouldBe` True
                        _ ->
                            expectationFailure "expected two rows"
                Left err ->
                    expectationFailure
                        ("expected acceptance, got " <> show err)

        it "reports a non-finite coordinate with its row and column" $
            case fromRows @(Finite 2)
                (chunksOf 2 [1, 0, 0, 1 / 0]) of
                Left err ->
                    err `shouldBe` InRow 1 (NonFiniteEntry 1)
                Right _ ->
                    expectationFailure "expected rejection"

        prop "identifies a row whose sum is invalid" $
            forAll (genTransitionRows 3) $ \matrix ->
                let invalid = bumpSmallestInFirstRow 1e-6 matrix
                 in case fromRows @(Finite 3) invalid of
                        Left (InRow 0 (SumOffBy _)) ->
                            property True
                        result ->
                            counterexample
                                ("expected InRow 0 SumOffBy, got " <> show result)
                                False

        prop "identifies a negative entry by row and column" $
            forAll (genTransitionRows 3) $ \matrix ->
                let invalid = setFirstEntry (-1e-6) matrix
                 in case fromRows @(Finite 3) invalid of
                        Left (InRow 0 (NegativeEntry 0 _)) ->
                            property True
                        result ->
                            counterexample
                                ("expected InRow 0 NegativeEntry 0, got " <> show result)
                                False

    describe "compose" $ do
        prop "is closed under multiplication"
            $ forAll
                ((,) <$> genTransitionRows 3 <*> genTransitionRows 3)
            $ \(left, right) ->
                case (fromRows @(Finite 3) left, fromRows @(Finite 3) right) of
                    (Right leftMatrix, Right rightMatrix) ->
                        case fromRows @(Finite 3)
                            (toRows (compose leftMatrix rightMatrix)) of
                            Right _ ->
                                property True
                            Left err ->
                                counterexample
                                    ("matrix product was rejected: " <> show err)
                                    False
                    result ->
                        counterexample
                            ("generated matrix was rejected: " <> show result)
                            False

        prop "approximately equals itself at zero tolerance" $
            forAll (genTransitionRows 3) $ \matrix ->
                case fromRows @(Finite 3) matrix of
                    Right transitionMatrix ->
                        property
                            ( approxTransitionMatrixEq
                                0
                                transitionMatrix
                                transitionMatrix
                            )
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

    describe "TransitionMatrix Semigroup" $ do
        prop "composition is approximately associative"
            $ forAll
                ( (,,)
                    <$> genTransitionRows 3
                    <*> genTransitionRows 3
                    <*> genTransitionRows 3
                )
            $ \(matrixA, matrixB, matrixC) ->
                case ( fromRows @(Finite 3) matrixA
                     , fromRows @(Finite 3) matrixB
                     , fromRows @(Finite 3) matrixC
                     ) of
                    (Right a, Right b, Right c) ->
                        property $
                            approxTransitionMatrixEq
                                1e-9
                                ((a <> b) <> c)
                                (a <> (b <> c))
                    result ->
                        counterexample
                            ("generated matrices were rejected: " <> show result)
                            False

    describe "TransitionMatrix Monoid" $ do
        prop "has a left identity" $
            forAll (genTransitionRows 3) $ \matrix ->
                case fromRows @(Finite 3) matrix of
                    Right p ->
                        property $
                            approxTransitionMatrixEq
                                1e-12
                                (mempty <> p)
                                p
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

        prop "has a right identity" $
            forAll (genTransitionRows 3) $ \matrix ->
                case fromRows @(Finite 3) matrix of
                    Right p ->
                        property $
                            approxTransitionMatrixEq
                                1e-12
                                (p <> mempty)
                                p
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

        it "uses the identity transition matrix as mempty" $
            approxTransitionMatrixEq
                1e-12
                (mempty :: TransitionMatrix (Finite 2))
                identity
                `shouldBe` True

    describe "power" $ do
        it "returns the identity at exponent zero" $
            approxTransitionMatrixEq
                1e-12
                (power 0 twoState)
                identity
                `shouldBe` True

        it "returns the matrix itself at exponent one" $
            approxTransitionMatrixEq
                1e-12
                (power 1 twoState)
                twoState
                `shouldBe` True

        it "matches a hand-computed square" $
            approxTransitionMatrixEq
                1e-9
                (power 2 twoState)
                twoStateSquared
                `shouldBe` True

        prop "stays stochastic for small exponents" $
            forAll ((,) <$> choose (0, 6 :: Int) <*> genTransitionRows 3) $
                \(k, matrix) ->
                    case fromRows @(Finite 3) matrix of
                        Right p ->
                            case fromRows @(Finite 3)
                                (toRows (power (fromIntegral k) p)) of
                                Right _ ->
                                    property True
                                Left err ->
                                    counterexample
                                        ("power left the stochastic set: " <> show err)
                                        False
                        Left err ->
                            counterexample
                                ("generated matrix was rejected: " <> show err)
                                False

        prop "satisfies the power addition law"
            $ forAll
                ( (,,)
                    <$> choose (0, 6 :: Int)
                    <*> choose (0, 6 :: Int)
                    <*> genTransitionRows 3
                )
            $ \(m, n, matrix) ->
                case fromRows @(Finite 3) matrix of
                    Right p ->
                        property $
                            approxTransitionMatrixEq
                                1e-9
                                (power (fromIntegral (m + n)) p)
                                ( power (fromIntegral m) p
                                    <> power (fromIntegral n) p
                                )
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

    describe "rowAt" $ do
        it "reads rows rather than columns" $
            Vector.toList (rowAt cyclicThree 0)
                `shouldBe` [0, 1, 0]

        it "returns each row of the three-cycle" $ do
            let row index = Vector.toList (rowAt cyclicThree index)

            row 0 `shouldBe` [0, 1, 0]
            row 1 `shouldBe` [0, 0, 1]
            row 2 `shouldBe` [1, 0, 0]

        prop "always returns a valid distribution" $
            forAll (genTransitionRows 3) $ \matrix ->
                case fromRows @(Finite 3) matrix of
                    Right transitionMatrix ->
                        conjoin
                            [ case Vector.fromList @(Finite 3)
                                (Vector.toList (rowAt transitionMatrix index)) of
                                Right _ ->
                                    property True
                                Left err ->
                                    counterexample
                                        ("row was rejected: " <> show err)
                                        False
                            | index <- [0 .. 2]
                            ]
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

    describe "named finite states" $ do
        it "returns a row labelled by named constructors" $
            distributionWeights (rowAt namedCycle PhaseA)
                `shouldBe` [(PhaseB, 1)]

        it "preserves the named state type through powers" $
            distributionWeights (rowAt (power 2 namedCycle) PhaseA)
                `shouldBe` [(PhaseC, 1)]

        it "provides a named identity matrix" $
            distributionWeights (rowAt (identity @NamedPhase) PhaseB)
                `shouldBe` [(PhaseB, 1)]

        it "composes matrices without changing their named state type" $
            approxTransitionMatrixEq
                0
                (namedCycle <> identity)
                namedCycle
                `shouldBe` True
