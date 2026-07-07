module Dtmc.Stochastic where

import Numeric.LinearAlgebra

newtype TransitionMatrix = TransitionMatrix (Matrix Double)

stationaryDistribution :: TransitionMatrix -> Vector Double
stationaryDistribution = undefined

