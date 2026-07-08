{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Simulation
  ( step
  ) where

import Control.Monad.Primitive (PrimMonad, PrimState)
import Data.Finite (Finite, finite, getFinite)
import Dtmc.StochasticMatrix (StochasticMatrix, unStochasticMatrix)
import GHC.TypeNats (KnownNat)
import Numeric.LinearAlgebra (toRows)
import qualified Numeric.LinearAlgebra as LA
import qualified System.Random.MWC as MWC
import qualified System.Random.MWC.Distributions as MWC

-- | Sample one transition of a discrete-time Markov chain.
--
-- Given a stochastic matrix @P@ and a current state @i@, this samples the next
-- state from row @i@ of @P@.
--
-- The function is polymorphic in 'PrimMonad', so it can run in 'IO' for real
-- simulations and in 'ST' for deterministic tests with a seeded generator.
--
-- The transition row is passed to categorical sampling as a vector of weights.
-- We do not re-normalise it: stochastic rows have already been validated up to
-- tolerance, and categorical sampling uses relative weights internally.
step
  :: forall n m
   . (KnownNat n, PrimMonad m)
  => StochasticMatrix n
  -> Finite n
  -> MWC.Gen (PrimState m)
  -> m (Finite n)
step matrix currentState gen = do
  let rowIndex =
        fromInteger (getFinite currentState)

      weights =
        transitionRow rowIndex matrix

  nextIndex <- MWC.categorical weights gen

  pure (finite @n (fromIntegral nextIndex))

transitionRow :: Int -> StochasticMatrix n -> LA.Vector Double
transitionRow rowIndex matrix =
  clampTinyNegativeEntries $
    toRows (unStochasticMatrix matrix) !! rowIndex

clampTinyNegativeEntries :: LA.Vector Double -> LA.Vector Double
clampTinyNegativeEntries =
  LA.cmap (\x -> if x < 0.0 then 0.0 else x)