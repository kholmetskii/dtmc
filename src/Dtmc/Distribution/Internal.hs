{- |
Module      : Dtmc.Distribution.Internal
Description : Raw carrier for probability distributions (unsafe underbelly).

Raw dimension-indexed carrier behind t'Dtmc.Distribution.Distribution'. The
public smart constructor validates the simplex invariant; internal callers may
use the constructor only when their operation preserves that invariant up to
floating-point error.

The constructor is positional so the public 'unDistribution' projection
cannot also act as a record-update setter.
-}
module Dtmc.Distribution.Internal (
    Distribution (Distribution),
    unDistribution,
) where

import GHC.TypeNats (
    KnownNat,
    Nat,
 )
import Numeric.LinearAlgebra.Static qualified as S

{- | A state distribution over @{0 .. n-1}@. Values built by
'Dtmc.Distribution.mkDistribution' satisfy its tolerant simplex check and
retain the original coordinates. The internal @Distribution@ constructor
performs no validation.
-}
newtype Distribution (n :: Nat)
    = -- | Build without simplex validation.
      Distribution (S.R n)

-- Nominal role: forbids 'Data.Coerce.coerce' from changing @n@. 'S.R' is
-- phantom in its size at runtime, so without this a @Distribution 3@ could be
-- coerced to a @Distribution 5@, breaking the dimension invariant.
type role Distribution nominal

deriving instance (KnownNat n) => Show (Distribution n)

{- | Return the stored probability vector unchanged. This is an @O(1)@
projection and performs no copy, validation, clamping, or renormalisation.
-}
unDistribution :: Distribution n -> S.R n
unDistribution (Distribution vector) = vector
