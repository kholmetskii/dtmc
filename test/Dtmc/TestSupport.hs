{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.TestSupport (
    testTolerance,
    approxEq,
    approxDistributionEq,
    approxTransitionMatrixEq,
    genSimplexPoint,
    genTransitionMatrix,
    modifyMatrixRows,
    bumpSmallest,
    bumpSmallestInFirstRow,
    setFirstEntry,
    hitProbabilityByState,
    hitEventualProbabilityByState,
    hitRaceProbabilityByState,
    hitExpectationByState,
    returnProbabilityByState,
    returnEventualProbabilityByState,
    returnExpectationByState,
    visitTotalProbabilityByState,
    visitInfiniteProbabilityByState,
    visitTotalExpectationByState,
    absorptionProbabilityByState,
    absorptionExpectationByState,
) where

import Data.Proxy (
    Proxy (..),
 )
import Dtmc.Analysis.Absorption qualified as Absorption
import Dtmc.Analysis.Event (
    DiscreteEvent,
 )
import Dtmc.Analysis.Expectation (
    Expectation,
 )
import Dtmc.Analysis.HittingTime qualified as Hit
import Dtmc.Analysis.LinearSystem (
    LinearSystemError,
 )
import Dtmc.Analysis.ReturnTime qualified as Return
import Dtmc.Analysis.VisitCount qualified as Visit
import Dtmc.Distribution.Vector (
    DistributionVector,
    unDistributionVector,
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
    finiteStates,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    unTransitionMatrix,
 )
import GHC.TypeNats (
    KnownNat,
    natVal,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Test.QuickCheck (
    Gen,
    choose,
    frequency,
    vectorOf,
 )

hitProbabilityByState ::
    forall state.
    (FiniteState state) =>
    DiscreteEvent ->
    TransitionMatrix state ->
    [state] ->
    S.R (Cardinality state)
hitProbabilityByState event matrix targets =
    S.vector
        [ Hit.probabilityGivenInitialState event matrix (`elem` targets) initial
        | initial <- finiteStates
        ]

hitEventualProbabilityByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    [state] ->
    Either LinearSystemError (S.R (Cardinality state))
hitEventualProbabilityByState matrix targets =
    S.vector
        <$> traverse
            (Hit.eventualProbabilityGivenInitialState matrix targets)
            finiteStates

hitRaceProbabilityByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    [state] ->
    [state] ->
    Either LinearSystemError (S.R (Cardinality state))
hitRaceProbabilityByState matrix successful competing =
    S.vector
        <$> traverse
            (Hit.raceProbabilityGivenInitialState matrix successful competing)
            finiteStates

hitExpectationByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    [state] ->
    Either LinearSystemError [Expectation]
hitExpectationByState matrix targets =
    traverse
        (Hit.expectationGivenInitialState matrix targets)
        finiteStates

returnProbabilityByState ::
    forall state.
    (FiniteState state) =>
    DiscreteEvent ->
    TransitionMatrix state ->
    S.R (Cardinality state)
returnProbabilityByState event matrix =
    S.vector
        [ Return.probabilityGivenInitialState event matrix initial
        | initial <- finiteStates
        ]

returnEventualProbabilityByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError (S.R (Cardinality state))
returnEventualProbabilityByState matrix =
    S.vector
        <$> traverse
            (Return.eventualProbabilityGivenInitialState matrix)
            finiteStates

returnExpectationByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError [Expectation]
returnExpectationByState matrix =
    traverse
        (Return.expectationGivenInitialState matrix)
        finiteStates

visitTotalProbabilityByState ::
    forall state.
    (FiniteState state) =>
    DiscreteEvent ->
    TransitionMatrix state ->
    state ->
    Either LinearSystemError (S.R (Cardinality state))
visitTotalProbabilityByState event matrix target =
    S.vector
        <$> traverse
            (Visit.totalProbabilityGivenInitialState event matrix target)
            finiteStates

visitInfiniteProbabilityByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError (S.R (Cardinality state))
visitInfiniteProbabilityByState matrix target =
    S.vector
        <$> traverse
            (Visit.infiniteProbabilityGivenInitialState matrix target)
            finiteStates

visitTotalExpectationByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError [Expectation]
visitTotalExpectationByState matrix target =
    traverse
        (Visit.totalExpectationGivenInitialState matrix target)
        finiteStates

absorptionProbabilityByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError (S.R (Cardinality state))
absorptionProbabilityByState matrix target =
    S.vector
        <$> traverse
            (Absorption.probabilityGivenInitialState matrix target)
            finiteStates

absorptionExpectationByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError [Expectation]
absorptionExpectationByState matrix =
    traverse
        (Absorption.expectationGivenInitialState matrix)
        finiteStates

{- | Absolute slack the tests use when comparing floating-point results. Kept
independent of the library's private validation threshold so a change there
cannot silently mask a regression here; the two happen to share a value.
-}
testTolerance :: Double
testTolerance = 1e-9

{- | Absolute-tolerance comparison of two scalar 'Double' results, matching the
@abs (x - y) <= tolerance@ convention of the vector and matrix helpers.
-}
approxEq :: Double -> Double -> Double -> Bool
approxEq tolerance left right =
    abs (left - right) <= tolerance

genSimplexPoint :: Int -> Gen [Double]
genSimplexPoint dimension = do
    entries <- vectorOf dimension genEntry
    let total = sum entries
    if total == 0
        then genSimplexPoint dimension
        else pure (map (/ total) entries)
  where
    genEntry =
        frequency
            [ (3, pure 0)
            , (7, choose (0, 1000))
            ]

genTransitionMatrix ::
    forall n.
    (KnownNat n) =>
    Gen (S.Sq n)
genTransitionMatrix = do
    rows <- vectorOf dimension (genSimplexPoint dimension)
    pure (S.matrix (concat rows))
  where
    dimension = fromIntegral (natVal (Proxy @n))

modifyMatrixRows ::
    (KnownNat n) =>
    ([[Double]] -> [[Double]]) ->
    S.Sq n ->
    S.Sq n
modifyMatrixRows transform =
    S.matrix
        . concat
        . transform
        . LA.toLists
        . S.extract

bumpSmallest :: Double -> [Double] -> [Double]
bumpSmallest _ [] = []
bumpSmallest amount entries =
    zipWith bump [0 :: Int ..] entries
  where
    smallestIndex =
        snd (minimum (zip entries [0 :: Int ..]))

    bump index entry
        | index == smallestIndex = entry + amount
        | otherwise = entry

bumpSmallestInFirstRow ::
    Double ->
    [[Double]] ->
    [[Double]]
bumpSmallestInFirstRow _ [] = []
bumpSmallestInFirstRow amount (row : rows) =
    bumpSmallest amount row : rows

setFirstEntry ::
    Double ->
    [[Double]] ->
    [[Double]]
setFirstEntry value ((_ : rest) : rows) =
    (value : rest) : rows
setFirstEntry _ rows = rows

approxTransitionMatrixEq ::
    (FiniteState state) =>
    Double ->
    TransitionMatrix state ->
    TransitionMatrix state ->
    Bool
approxTransitionMatrixEq tolerance left right =
    and (zipWith close (entries left) (entries right))
  where
    entries = LA.toList . LA.flatten . S.extract . unTransitionMatrix
    close x y = abs (x - y) <= tolerance

approxDistributionEq ::
    (FiniteState state) =>
    Double ->
    DistributionVector state ->
    DistributionVector state ->
    Bool
approxDistributionEq tolerance left right =
    and (zipWith close (entries left) (entries right))
  where
    entries = LA.toList . S.extract . unDistributionVector
    close x y = abs (x - y) <= tolerance
