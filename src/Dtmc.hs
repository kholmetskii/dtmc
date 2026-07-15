-- |
-- Module      : Dtmc
-- Description : Public facade re-exporting the library's curated API.
--
-- Single entry point for users of the library. It gathers the intended public
-- surface -- the 'Distribution' and 'TransitionMatrix' types with their
-- validating constructors and error types, the analytic forward dynamics
-- ('evolve' and 'matrixPower'), the qualitative structure
-- theory ('communicatingClasses', 'irreducible', 'period', 'classify'), the
-- random simulation primitives, and the shared numeric tolerance with its
-- approximate-equality helpers -- while hiding the "Dtmc.Internal" modules.
-- Import this module to build, analyse, and run chains.
module Dtmc (
    Distribution,
    DistributionError (..),
    mkDistribution,
    unDistribution,
    TransitionMatrix,
    TransitionError (..),
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
    SupportGraph,
    supportGraphOf,
    supportEdge,
    supportEdgeIn,
    accessible,
    accessibleIn,
    communicates,
    communicatesIn,
    communicatingClasses,
    communicatingClassesIn,
    irreducible,
    irreducibleIn,
    period,
    periodIn,
    aperiodic,
    aperiodicIn,
    CommClass (..),
    Classification,
    classesOf,
    classify,
    classifyIn,
    Irreducible,
    witnessIrreducible,
    irreducibleMatrix,
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
import Dtmc.Internal.Simplex (
    SimplexError (..),
 )
import Dtmc.Simulation (
    sampleFrom,
    step,
 )
import Dtmc.TransitionMatrix (
    TransitionError (..),
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
    SupportGraph,
    accessible,
    accessibleIn,
    aperiodic,
    aperiodicIn,
    classesOf,
    classify,
    classifyIn,
    communicates,
    communicatesIn,
    communicatingClasses,
    communicatingClassesIn,
    irreducible,
    irreducibleIn,
    irreducibleMatrix,
    period,
    periodIn,
    supportEdge,
    supportEdgeIn,
    supportGraphOf,
    witnessIrreducible,
 )

import Dtmc.Approx (
    approxEq,
    approxEqDist,
    approxEqMatrix,
    tolerance,
 )