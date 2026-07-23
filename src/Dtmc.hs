-- |
-- Module      : Dtmc
-- Description : Public facade re-exporting the library's curated API.
--
-- Single entry point for users of the library. It gathers the intended public
-- surface -- the 'Distribution' and 'TransitionMatrix' types with their
-- validating constructors and error types, the analytic forward dynamics
-- ('evolve' and 'matrixPower'), the qualitative structure
-- theory ('communicatingClasses', 'irreducible', 'period', 'classify', 'recurrentStates'), the
-- hitting theory ('hittingProbabilities', 'hittingProbability',
-- 'expectedHittingTimes', 'expectedHittingTime', 'returnProbabilities',
-- 'returnProbability', 'expectedReturnTimes', and 'expectedReturnTime'), the
-- random simulation primitives, and the shared numeric tolerance with its
-- approximate-equality helpers -- while hiding the "Dtmc.Internal" modules.
-- Import this module to build, analyse, and run chains.
module Dtmc (
    Distribution,
    DistributionError (..),
    mkDistribution,
    unDistribution,
    TransitionMatrix,
    TransitionMatrixError (..),
    mkTransitionMatrix,
    unTransitionMatrix,
    mulTransitionMatrix,
    rowAt,
    SimplexError (..),
    sampleFrom,
    step,
    evolve,
    evolveN,
    identityMatrix,
    matrixPower,
    supportEdge,
    accessible,
    communicates,
    communicatingClasses,
    irreducible,
    period,
    aperiodic,
    cyclicClasses,
    recurrentState,
    transientState,
    recurrentStates,
    transientStates,
    CommClass (..),
    Classification,
    classesOf,
    isIrreducible,
    isAperiodic,
    isErgodic,
    chainPeriod,
    recurrentStatesOf,
    transientStatesOf,
    absorbingStates,
    classify,
    Irreducible,
    witnessIrreducible,
    unIrreducible,
    MeanTime (..),
    hittingProbabilities,
    hittingProbability,
    expectedHittingTimes,
    expectedHittingTime,
    returnProbabilities,
    returnProbability,
    expectedReturnTimes,
    expectedReturnTime,
    tolerance,
    approxEq,
    approxEqDist,
    approxEqMatrix,
) where

import Dtmc.Distribution (
    Distribution,
    DistributionError (..),
    mkDistribution,
    unDistribution,
 )
import Dtmc.Simplex (
    SimplexError (..),
 )
import Dtmc.Simulation (
    sampleFrom,
    step,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrixError (..),
    TransitionMatrix,
    mkTransitionMatrix,
    mulTransitionMatrix,
    rowAt,
    unTransitionMatrix,
 )

import Dtmc.Dynamics (
    evolve,
    evolveN,
    identityMatrix,
    matrixPower,
 )

import Dtmc.Classification (
    Classification,
    CommClass (..),
    Irreducible,
    absorbingStates,
    accessible,
    aperiodic,
    chainPeriod,
    classesOf,
    classify,
    communicates,
    communicatingClasses,
    cyclicClasses,
    irreducible,
    isAperiodic,
    isErgodic,
    isIrreducible,
    period,
    recurrentState,
    recurrentStates,
    recurrentStatesOf,
    supportEdge,
    transientState,
    transientStates,
    transientStatesOf,
    unIrreducible,
    witnessIrreducible,
 )

import Dtmc.Approx (
    approxEq,
    approxEqDist,
    approxEqMatrix,
    tolerance,
 )

import Dtmc.Hitting (
    MeanTime (..),
    expectedHittingTime,
    expectedHittingTimes,
    expectedReturnTime,
    expectedReturnTimes,
    hittingProbability,
    hittingProbabilities,
    returnProbability,
    returnProbabilities,
 )
