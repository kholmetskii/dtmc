-- |
-- Module      : Dtmc
-- Description : Public facade re-exporting the library's curated API.
--
-- Single entry point for users of the library. It gathers the intended public
-- surface -- the 'Distribution' and 'TransitionMatrix' types with their
-- validating constructors and error types, the analytic forward dynamics
-- ('evolve' and 'matrixPower'), the qualitative structure
-- theory ('communicatingClasses', 'irreducible', 'period', 'classify'), and the
-- random simulation primitives -- while hiding the "Dtmc.Internal" modules.
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
    supportEdge,
    accessible,
    communicates,
    communicatingClasses,
    irreducible,
    period,
    aperiodic,
    CommClass (..),
    Classification,
    classesOf,
    classify,
    Irreducible,
    witnessIrreducible,
    irreducibleMatrix,
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
    accessible,
    aperiodic,
    classesOf,
    classify,
    communicatingClasses,
    communicates,
    irreducible,
    irreducibleMatrix,
    period,
    supportEdge,
    witnessIrreducible,
 )
