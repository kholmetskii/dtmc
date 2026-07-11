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
    chapmanKolmogorov
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
    evolve, 
    evolveN, 
    identityMatrix, 
    matrixPower, 
    chapmanKolmogorov
    )