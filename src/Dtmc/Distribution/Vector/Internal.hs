{- |
Module      : Dtmc.Distribution.Vector.Internal
Description : Unsafe carrier for dense distribution vectors.

Public smart constructors validate the simplex invariant. Internal callers may
use the constructor only when their operation preserves that invariant up to
floating-point error.
-}
module Dtmc.Distribution.Vector.Internal (
    DistributionVector (DistributionVector),
    unDistributionVector,
) where

import Data.Finite (
    Finite,
    finites,
    getFinite,
 )
import Dtmc.Distribution (
    Distribution (..),
 )
import GHC.TypeNats (
    KnownNat,
    Nat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

{- | A state distribution vector over @{0 .. n-1}@. The internal constructor
performs no validation.
-}
newtype DistributionVector (n :: Nat)
    = DistributionVector (S.R n)

-- Nominal role prevents changing the static dimension through 'Data.Coerce'.
type role DistributionVector nominal

deriving instance (KnownNat n) => Show (DistributionVector n)

{- | Return the stored probability vector unchanged. This is an @O(1)@
projection and performs no copy, validation, clamping, or renormalisation.
-}
unDistributionVector :: DistributionVector n -> S.R n
unDistributionVector (DistributionVector vector) = vector

instance (KnownNat n) => Distribution (DistributionVector n) where
    type DistributionState (DistributionVector n) = Finite n

    probabilityAt distribution index =
        S.extract (unDistributionVector distribution)
            `LA.atIndex` fromIntegral (getFinite index)

    distributionWeights distribution =
        [ (state, weight)
        | (state, weight) <- zip finites weights
        , weight /= 0
        ]
      where
        weights = LA.toList (S.extract (unDistributionVector distribution))
