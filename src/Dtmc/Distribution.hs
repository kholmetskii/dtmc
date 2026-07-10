module Dtmc.Distribution
  ( Distribution
  , unDistribution
  , SimplexError (..)
  , DistributionError (..)
  , simplexTolerance
  , validateSimplex
  , mkDistribution
  , approxDistributionEq
  ) where

import Data.Bifunctor (first)
import Dtmc.Internal.Simplex
  ( SimplexError (..)
  , simplexTolerance
  , validateSimplex
  )
import Dtmc.Internal.Types (Distribution (..))
import GHC.TypeNats (KnownNat)
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Static as S

newtype DistributionError = DistributionError SimplexError
  deriving (Eq, Show)

mkDistribution :: KnownNat n => S.R n -> Either DistributionError (Distribution n)
mkDistribution vector =
  Distribution vector <$ first DistributionError (validateSimplex vector)

approxDistributionEq :: KnownNat n => Double -> Distribution n -> Distribution n -> Bool
approxDistributionEq tolerance (Distribution left) (Distribution right) =
  and (zipWith close (entries left) (entries right))
  where
    entries = LA.toList . S.extract
    close x y = abs (x - y) <= tolerance
