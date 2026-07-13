module Dtmc (
    Distribution,
    DistributionError (..),
    mkDistribution,
    unDistribution,
    approxDistributionEq,
    TransitionMatrix,
    TransitionError (..),
    mkTransitionMatrix,
    unTransitionMatrix,
    mulTransitionMatrix,
    rowAt,
    approxTransitionMatrixEq,
    SimplexError (..),
    sampleFrom,
    step,
    evolve,
    evolveN,
    identityMatrix,
    matrixPower,
    chapmanKolmogorov,
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
    approxDistributionEq,
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
    approxTransitionMatrixEq,
    mkTransitionMatrix,
    mulTransitionMatrix,
    rowAt,
    unTransitionMatrix,
 )

import Dtmc.Dynamics (
    chapmanKolmogorov,
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
