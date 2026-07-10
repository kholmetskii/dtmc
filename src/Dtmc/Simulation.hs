module Dtmc.Simulation
  ( sampleFrom
  , step
  ) where

import Control.Monad.Primitive (PrimMonad, PrimState)
import Data.Finite (Finite, finite)
import Dtmc.Distribution
  ( Distribution
  , simplexTolerance
  , unDistribution
  )
import Dtmc.TransitionMatrix
  ( TransitionMatrix
  , rowAt
  )
import GHC.TypeNats (KnownNat)
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Static as S
import qualified System.Random.MWC as MWC
import qualified System.Random.MWC.Distributions as MWCD

sampleFrom :: (KnownNat n, PrimMonad m) => Distribution n -> MWC.Gen (PrimState m) -> m (Finite n)
sampleFrom distribution generator = do
  index <- MWCD.categorical weights generator
  pure (finite (fromIntegral index))
  where
    weights =
      clampToleratedNegatives
        (S.extract (unDistribution distribution))

step :: (KnownNat n, PrimMonad m) => TransitionMatrix n -> Finite n -> MWC.Gen (PrimState m) -> m (Finite n)
step matrix state =
  sampleFrom (rowAt matrix state)

clampToleratedNegatives :: LA.Vector Double -> LA.Vector Double
clampToleratedNegatives =
  LA.cmap clamp
  where
    clamp value
      | value >= 0 = value
      | value >= negate simplexTolerance = 0
      | otherwise =
          error
            ( "Dtmc.Simulation: coordinate "
                <> show value
                <> " is below -simplexTolerance"
            )
