{- |
Module      : Dtmc.Distribution.Internal
Description : Raw carrier for probability distributions (unsafe underbelly).

Raw carriers behind t'Dtmc.Distribution.DistributionVector' and
t'Dtmc.Distribution.SparseDistribution'. Public smart constructors validate
the simplex invariant; internal callers may use the constructors only when
their operation preserves that invariant up to floating-point error.

The constructor is positional so the public 'unDistributionVector' projection
cannot also act as a record-update setter.
-}
module Dtmc.Distribution.Internal (
    DistributionVector (DistributionVector),
    unDistributionVector,
    SparseDistribution (SparseDistribution),
    unSparseDistribution,
) where

import Data.Map.Strict (
    Map,
 )
import GHC.TypeNats (
    KnownNat,
    Nat,
 )
import Numeric.LinearAlgebra.Static qualified as S

{- | A state distribution vector over @{0 .. n-1}@. Values built by
'Dtmc.Distribution.mkDistributionVector' satisfy its tolerant simplex check
and retain the original coordinates. The internal @DistributionVector@
constructor performs no validation.
-}
newtype DistributionVector (n :: Nat)
    = -- | Build without simplex validation.
      DistributionVector (S.R n)

-- Nominal role: forbids 'Data.Coerce.coerce' from changing @n@. 'S.R' is
-- phantom in its size at runtime, so without this a @DistributionVector 3@
-- could be coerced to a @DistributionVector 5@, breaking the dimension
-- invariant.
type role DistributionVector nominal

deriving instance (KnownNat n) => Show (DistributionVector n)

{- | Return the stored probability vector unchanged. This is an @O(1)@
projection and performs no copy, validation, clamping, or renormalisation.
-}
unDistributionVector :: DistributionVector n -> S.R n
unDistributionVector (DistributionVector vector) = vector

{- | A finite-support probability distribution over an otherwise unrestricted
state type. The internal constructor performs no validation.
-}
newtype SparseDistribution state
    = SparseDistribution (Map state Double)

type role SparseDistribution nominal

deriving instance (Eq state) => Eq (SparseDistribution state)
deriving instance (Show state) => Show (SparseDistribution state)

-- | Return the canonical state-to-weight map without copying or validation.
unSparseDistribution :: SparseDistribution state -> Map state Double
unSparseDistribution (SparseDistribution weights) = weights
