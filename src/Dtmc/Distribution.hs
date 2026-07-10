-- |
-- Module      : Dtmc.Distribution
--
-- A probability distribution over a finite state space.
--
-- This is a validation boundary: the 'Distribution' constructor is applied
-- here but not re-exported, so 'mkDistribution' is the only door in.
--
-- This module does NOT import 'Dtmc.StochasticMatrix', nor is it imported by
-- it. The shared predicate lives below, in 'Dtmc.Simplex'. 'Dtmc.Kernel' is
-- what connects them. See docs/DECISIONS.md D3.
module Dtmc.Distribution
  ( Distribution
  , unDistribution
  , DistributionError (..)
  , mkDistribution
  , approxDistributionEq
  ) where

import Data.Bifunctor (first)
import Dtmc.Internal (Distribution (..))
import Dtmc.Simplex (SimplexError, validateSimplexPoint)
import GHC.TypeNats (KnownNat)
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Static as S

-- | A wrapper over 'SimplexError'. No row index: a standalone distribution has
-- no rows.
newtype DistributionError = DistributionError SimplexError
  deriving (Eq, Show)

-- | Validates but does NOT mutate: no normalising, no clamping. That is why
-- the round-trip property is exact equality rather than approximate.
--
-- Verified: Dtmc.DistributionSpec
mkDistribution :: KnownNat n => S.R n -> Either DistributionError (Distribution n)
mkDistribution v =
  Distribution { unDistribution = v } <$ first DistributionError (validateSimplexPoint v)

-- | Coordinatewise comparison with an EXPLICIT tolerance.
--
-- The tolerance is a parameter, not 'Dtmc.Simplex.simplexTolerance': the
-- project has two coexisting regimes six orders of magnitude apart
-- (NUMERICS N6). An explicit argument forces the caller to name which one.
--
-- @and (zipWith ...)@ rather than @LA.maxElement (abs (a - b))@: the latter
-- throws on an empty vector (n = 0).
approxDistributionEq
  :: KnownNat n => Double -> Distribution n -> Distribution n -> Bool
approxDistributionEq tol (Distribution a) (Distribution b) =
  and (zipWith close (LA.toList (S.extract a)) (LA.toList (S.extract b)))
  where
    close x y = abs (x - y) <= tol
