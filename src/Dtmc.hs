{-# LANGUAGE ExplicitNamespaces #-}

{- |
Module      : Dtmc
Description : Public facade re-exporting the library's curated API.

Curated entry point for the library's shared model-building API. It gathers
states, distributions, transitions, forward dynamics, simulation,
classification, stationary analysis, and common analysis result and error
types while hiding implementation modules.

Probability queries are intentionally imported from focused modules and used
qualified: "Dtmc.Analysis.FiniteTime", "Dtmc.Analysis.HittingTime",
"Dtmc.Analysis.ReturnTime", and "Dtmc.Analysis.VisitCount". This keeps their
short canonical names unambiguous and makes the analyzed random quantity
visible at each call site.
-}
module Dtmc (
    type Cardinality,
    LinearSystemError (..),
    FiniteState,
    finiteStates,
    stateIndex,
    stateAt,
    Distribution (..),
    DistributionVector,
    DistributionError (..),
    mkDistributionVector,
    unDistributionVector,
    DistributionMap,
    mkDistributionMap,
    unDistributionMap,
    toDistributionMap,
    pointMass,
    TransitionMatrix,
    TransitionMatrixError (..),
    mkTransitionMatrix,
    unTransitionMatrix,
    mulTransitionMatrix,
    rowAt,
    SimplexError (..),
    sample,
    step,
    simulateN,
    evolve,
    evolveN,
    evolveVector,
    evolveVectorN,
    Observation (..),
    ConditionalProbabilityError (..),
    identityMatrix,
    matrixPower,
    supportEdge,
    accessible,
    reachesAny,
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
    type CommClass (..),
    type Classification,
    classesOf,
    isIrreducible,
    isAperiodic,
    isErgodic,
    chainPeriod,
    recurrentStatesOf,
    transientStatesOf,
    absorbingStates,
    classify,
    type Irreducible,
    witnessIrreducible,
    unIrreducible,
    Expectation (..),
    stationaryDistribution,
    Transition (..),
    TransitionKernel,
    transitionKernel,
    deterministicKernel,
) where

import Dtmc.Analysis.Expectation (
    Expectation (..),
 )
import Dtmc.Analysis.FiniteTime (
    ConditionalProbabilityError (..),
    Observation (..),
 )
import Dtmc.Analysis.Stationary (
    stationaryDistribution,
 )
import Dtmc.Distribution (
    Distribution (..),
    DistributionError (..),
 )
import Dtmc.Distribution.Map (
    DistributionMap,
    mkDistributionMap,
    pointMass,
    toDistributionMap,
    unDistributionMap,
 )
import Dtmc.Distribution.Vector (
    DistributionVector,
    mkDistributionVector,
    unDistributionVector,
 )
import Dtmc.Simplex (
    SimplexError (..),
 )
import Dtmc.Simulation (
    sample,
    simulateN,
    step,
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
    finiteStates,
    stateAt,
    stateIndex,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    TransitionMatrixError (..),
    identityMatrix,
    matrixPower,
    mkTransitionMatrix,
    mulTransitionMatrix,
    rowAt,
    unTransitionMatrix,
 )

import Dtmc.Dynamics (
    evolve,
    evolveN,
    evolveVector,
    evolveVectorN,
 )

import Dtmc.Transition (
    Transition (..),
 )
import Dtmc.Transition.Kernel (
    TransitionKernel,
    deterministicKernel,
    transitionKernel,
 )

import Dtmc.Analysis.Classification (
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
    reachesAny,
    recurrentState,
    recurrentStates,
    recurrentStatesOf,
    supportEdge,
    transientState,
    transientStates,
    transientStatesOf,
    unIrreducible,
    witnessIrreducible,
    type Classification,
    type CommClass (..),
    type Irreducible,
 )

import Dtmc.Analysis.LinearSystem (
    LinearSystemError (..),
 )
