module Dtmc
  ( Distribution
  , DistributionError (..)
  , mkDistribution
  , unDistribution
  , TransitionMatrix
  , TransitionError (..)
  , mkTransitionMatrix
  , unTransitionMatrix
  , rowAt
  , mulTransitionMatrix
  , sampleFrom
  , step
  ) where

import Dtmc.Distribution
    ( DistributionError(..), Distribution(..), mkDistribution )
import Dtmc.Simulation ( sampleFrom, step )
import Dtmc.TransitionMatrix
    ( TransitionError(..),
      TransitionMatrix(..),
      mkTransitionMatrix,
      mulTransitionMatrix,
      rowAt )
