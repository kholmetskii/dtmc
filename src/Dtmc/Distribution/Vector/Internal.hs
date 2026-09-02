{- |
Module      : Dtmc.Distribution.Vector.Internal
Description : Unsafe carrier for dense distribution vectors.

The public smart constructor validates and canonicalises the simplex
invariant. Internal callers may use the constructor only when their operation
preserves that invariant up to floating-point error.
-}
module Dtmc.Distribution.Vector.Internal (
    DistributionVector (DistributionVector),
    unDistributionVector,
) where

import Dtmc.Distribution (
    Distribution (..),
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
    finiteStates,
 )
import Dtmc.State.Internal (
    stateIndexInt,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

{- | A state distribution vector whose coordinates follow the canonical order
of its finite state type. The internal constructor performs no validation.
-}
newtype DistributionVector state
    = DistributionVector (S.R (Cardinality state))

-- Nominal role prevents coercion between distinct state types, including
-- state types with the same cardinality.
type role DistributionVector nominal

deriving instance (FiniteState state) => Show (DistributionVector state)

{- | Return the stored probability vector unchanged. This is an @O(1)@
projection and performs no copy, validation, clamping, or renormalisation.
-}
unDistributionVector ::
    DistributionVector state ->
    S.R (Cardinality state)
unDistributionVector (DistributionVector vector) = vector

instance (FiniteState state) => Distribution (DistributionVector state) where
    type DistributionState (DistributionVector state) = state

    probabilityAt distribution state =
        S.extract (unDistributionVector distribution)
            `LA.atIndex` stateIndexInt state

    distributionWeights distribution =
        [ (state, weight)
        | (state, weight) <- zip finiteStates weights
        , weight /= 0
        ]
      where
        weights = LA.toList (S.extract (unDistributionVector distribution))
