-- |
-- Module      : Dtmc.Simulation
--
-- The simulation primitive. Does not import 'Dtmc.Internal': it constructs
-- only through 'Dtmc.Kernel.rowAt'.
module Dtmc.Simulation
  ( sampleFrom
  , step
  ) where

import Control.Monad.Primitive (PrimMonad, PrimState)
import Data.Finite (Finite, finite)
import Dtmc.Internal ( Distribution, unDistribution, StochasticMatrix )
import Dtmc.Kernel (rowAt)
import Dtmc.Simplex (simplexTolerance)
import GHC.TypeNats (KnownNat)
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Static as S
import qualified System.Random.MWC as MWC
import qualified System.Random.MWC.Distributions as MWCD

-- | Sample a state from a distribution.
--
-- Compare with a hypothetical @sampleFrom :: R n -> ...@: that would mean
-- "sample from an arbitrary vector of numbers" and would require the caller to
-- remember about normalisation. With 'Distribution' the signature means
-- "sample from a distribution", and that is TRUE BY TYPE.
--
-- Polymorphic in 'PrimMonad' rather than fixed to 'IO': tests run in
-- 'Control.Monad.ST.ST' with a fixed seed and are therefore deterministic
-- (docs/TESTING.md T5).
--
-- Proof (totality of 'finite'):
--   @MWCD.categorical@ returns an index in [0, length weights).
--   @weights@ is the extraction of an @R n@, so its length is exactly n.
--   Hence 0 ≤ j < n, and the partial function @finite@ cannot fail. ∎
--
-- @categorical@ normalises the weights internally, so a sum of 1 ± ε is
-- harmless and must NOT be renormalised (NUMERICS N4).
--
-- Verified: Dtmc.SimulationSpec
sampleFrom
  :: (KnownNat n, PrimMonad m)
  => Distribution n -> MWC.Gen (PrimState m) -> m (Finite n)
sampleFrom d gen = do
  j <- MWCD.categorical weights gen
  pure (finite (fromIntegral j))
  where
    weights = clampToleratedNegatives (S.extract (unDistribution d))

-- | One step of the chain. Reads as the definition: take row i as a
-- distribution, sample from it.
step
  :: (KnownNat n, PrimMonad m)
  => StochasticMatrix n -> Finite n -> MWC.Gen (PrimState m) -> m (Finite n)
step p i = sampleFrom (rowAt p i)

-- | Zeroes ONLY those negative coordinates that
-- 'Dtmc.Simplex.validateSimplexPoint' tolerates, i.e. those in [-ε, 0).
--
-- Anything below is a violated invariant, i.e. a bug, not data. The previous
-- version was named @clampTinyNegativeEntries@ and silently clamped ANY
-- negative: the name lied, and a real bug would have been masked.
clampToleratedNegatives :: LA.Vector Double -> LA.Vector Double
clampToleratedNegatives = LA.cmap clamp
  where
    clamp x
      | x >= 0 = x
      | x >= negate simplexTolerance = 0
      | otherwise =
          error $
            "Dtmc.Simulation: coordinate " <> show x
              <> " < -simplexTolerance; Distribution invariant violated"
