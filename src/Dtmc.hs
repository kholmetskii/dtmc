module Dtmc
  ( Distribution
  , DistributionError (..)
  , mkDistribution
  , unDistribution
  , approxDistributionEq
  , TransitionMatrix
  , TransitionError (..)
  , mkTransitionMatrix
  , unTransitionMatrix
  , mulTransitionMatrix
  , rowAt
  , approxTransitionMatrixEq
  , SimplexError (..)
  , sampleFrom
  , step
  ) where

import Dtmc.Distribution
  ( Distribution
  , DistributionError (..)
  , mkDistribution
  , unDistribution
  , approxDistributionEq
  )
import Dtmc.TransitionMatrix
  ( TransitionMatrix
  , TransitionError (..)
  , mkTransitionMatrix
  , unTransitionMatrix
  , mulTransitionMatrix
  , rowAt
  , approxTransitionMatrixEq
  )
import Dtmc.Internal.Simplex
  ( SimplexError (..)
  )
import Dtmc.Simulation
  ( sampleFrom
  , step
  )