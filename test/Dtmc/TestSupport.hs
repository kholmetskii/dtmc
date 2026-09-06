{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.TestSupport (
    testTolerance,
    approxEq,
    approxDistributionEq,
    approxTransitionMatrixEq,
    genSimplexPoint,
    genTransitionRows,
    chunksOf,
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
    toList,
 )
import Dtmc.State (
    FiniteState,
    finiteStates,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    toRows,
 )
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
    [Double]
hitProbabilityByState event matrix targets =
    [ Hit.probabilityGivenInitialState event matrix (`elem` targets) initial
    | initial <- finiteStates
    ]

hitEventualProbabilityByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    [state] ->
    Either LinearSystemError [Double]
hitEventualProbabilityByState matrix targets =
    traverse
        (Hit.eventualProbabilityGivenInitialState matrix targets)
        finiteStates

hitRaceProbabilityByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    [state] ->
    [state] ->
    Either LinearSystemError [Double]
hitRaceProbabilityByState matrix successful competing =
    traverse
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
    [Double]
returnProbabilityByState event matrix =
    [ Return.probabilityGivenInitialState event matrix initial
    | initial <- finiteStates
    ]

returnEventualProbabilityByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    Either LinearSystemError [Double]
returnEventualProbabilityByState matrix =
    traverse
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
    Either LinearSystemError [Double]
visitTotalProbabilityByState event matrix target =
    traverse
        (Visit.totalProbabilityGivenInitialState event matrix target)
        finiteStates

visitInfiniteProbabilityByState ::
    forall state.
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    Either LinearSystemError [Double]
visitInfiniteProbabilityByState matrix target =
    traverse
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
    Either LinearSystemError [Double]
absorptionProbabilityByState matrix target =
    traverse
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

{- | Generate a square grid of weights whose rows are probability vectors,
ready for 'Dtmc.Transition.Matrix.fromRows'.
-}
genTransitionRows :: Int -> Gen [[Double]]
genTransitionRows dimension =
    vectorOf dimension (genSimplexPoint dimension)

-- | Split a flat row-major list into rows of the given width.
chunksOf :: Int -> [value] -> [[value]]
chunksOf width values
    | width <= 0 || null values = []
    | otherwise = row : chunksOf width rest
  where
    (row, rest) = splitAt width values

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
    Double ->
    TransitionMatrix state ->
    TransitionMatrix state ->
    Bool
approxTransitionMatrixEq tolerance left right =
    and (zipWith close (entries left) (entries right))
  where
    entries = concat . toRows
    close x y = abs (x - y) <= tolerance

approxDistributionEq ::
    Double ->
    DistributionVector state ->
    DistributionVector state ->
    Bool
approxDistributionEq tolerance left right =
    and (zipWith close (entries left) (entries right))
  where
    entries = toList
    close x y = abs (x - y) <= tolerance
