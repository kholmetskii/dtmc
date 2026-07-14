-- |
-- Module      : Dtmc.Distribution
-- Description : Probability distributions over a finite state space.
--
-- Public interface for 'Distribution', the initial or marginal state
-- distribution of a chain on @n@ states. Values are built through the validating
-- 'mkDistribution', so anything of type @Distribution n@ is a genuine point on
-- the standard @(n-1)@-simplex.
module Dtmc.Distribution (
    Distribution,
    DistributionError (..),
    mkDistribution,
    unDistribution,
    approxDistributionEq,
) where

import Data.Bifunctor (
    first,
 )
import Dtmc.Internal.Simplex (
    SimplexError,
    validateSimplex,
 )
import Dtmc.Internal.Types (
    Distribution (..),
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

-- | A distribution-specific wrapper around 'SimplexError', keeping these
-- failures distinct in the type from the row failures of "Dtmc.TransitionMatrix".
newtype DistributionError = DistributionError SimplexError
    deriving (Eq, Show)

-- | Smart constructor: accept a raw vector only if it satisfies the simplex
-- invariant, otherwise return the reason. The single sanctioned way to obtain a
-- 'Distribution'. (@<$@ keeps the validated vector on success, the error on failure.)
mkDistribution :: (KnownNat n) => S.R n -> Either DistributionError (Distribution n)
mkDistribution vector =
    Distribution vector <$ first DistributionError (validateSimplex vector)

-- | Entrywise comparison within an absolute @tolerance@. Distributions hold
-- floating-point data, so this -- rather than exact @(==)@ -- is the right notion
-- of equality for tests and numerical results.
approxDistributionEq :: (KnownNat n) => Double -> Distribution n -> Distribution n -> Bool
approxDistributionEq tolerance (Distribution left) (Distribution right) =
    and (zipWith close (entries left) (entries right))
  where
    entries = LA.toList . S.extract
    close x y = abs (x - y) <= tolerance
