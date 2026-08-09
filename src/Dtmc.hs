{-# LANGUAGE ExplicitNamespaces #-}

{- |
Module      : Dtmc
Description : Public facade re-exporting the library's curated API.

Single entry point for users of the library. It gathers the intended public
surface -- the 'Distribution' abstraction, its 'DistributionVector' and
'DistributionMap' representations, the 'Transition' abstraction and its
'TransitionMatrix' and 'TransitionKernel' representations, their validating
constructors and error types, and the forward dynamics ('evolve' and
'matrixPower'), the scalar probability queries ('probabilityAt',
'transitionProbability', 'transitionProbabilityN', 'probabilityAtTime', and
'pathProbability'), the timed event queries ('probability' and
'conditionalProbability' over 'Observation' values), the qualitative structure
theory ('communicatingClasses', 'irreducible', 'period', 'classify', 'recurrentStates'), the
bounded, eventual, and expected hitting and first-return queries in
"Dtmc.Hitting", and the
random simulation primitives -- while hiding internal modules. Import this
module to build, analyse, and run chains. Focused imports are available through
"Dtmc.Distribution", "Dtmc.Distribution.Vector", "Dtmc.Distribution.Map",
"Dtmc.Transition", "Dtmc.Transition.Matrix", and "Dtmc.Transition.Kernel".
-}
module Dtmc (
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
    transitionProbability,
    transitionProbabilityN,
    SimplexError (..),
    sample,
    step,
    simulateN,
    evolve,
    evolveN,
    evolveVector,
    evolveVectorN,
    probabilityAtTime,
    pathProbability,
    probability,
    conditionalProbability,
    Observation (..),
    FiniteObservation,
    ProbabilityError (..),
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
    MeanTime (..),
    hittingTimeProbabilitiesAt,
    hittingTimeProbabilityAt,
    hittingTimeProbabilitiesBefore,
    hittingTimeProbabilityBefore,
    hittingProbabilities,
    hittingProbability,
    hittingBeforeProbabilities,
    hittingBeforeProbability,
    expectedHittingTimes,
    expectedHittingTime,
    returnTimeProbabilitiesAt,
    returnTimeProbabilityAt,
    returnTimeProbabilitiesBefore,
    returnTimeProbabilityBefore,
    returnProbabilities,
    returnProbability,
    expectedReturnTimes,
    expectedReturnTime,
    Transition (..),
    TransitionKernel,
    transitionKernel,
    deterministicKernel,
) where

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
import Dtmc.Probability (
    FiniteObservation,
    Observation (..),
    ProbabilityError (..),
    conditionalProbability,
    pathProbability,
    probability,
    probabilityAtTime,
    transitionProbability,
    transitionProbabilityN,
 )
import Dtmc.Simplex (
    SimplexError (..),
 )
import Dtmc.Simulation (
    sample,
    simulateN,
    step,
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

import Dtmc.Classification (
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
    type Classification,
    type CommClass (..),
    type Irreducible,
 )

import Dtmc.Hitting (
    MeanTime (..),
    expectedHittingTime,
    expectedHittingTimes,
    expectedReturnTime,
    expectedReturnTimes,
    hittingBeforeProbabilities,
    hittingBeforeProbability,
    hittingProbabilities,
    hittingProbability,
    hittingTimeProbabilitiesAt,
    hittingTimeProbabilitiesBefore,
    hittingTimeProbabilityAt,
    hittingTimeProbabilityBefore,
    returnProbabilities,
    returnProbability,
    returnTimeProbabilitiesAt,
    returnTimeProbabilitiesBefore,
    returnTimeProbabilityAt,
    returnTimeProbabilityBefore,
 )
