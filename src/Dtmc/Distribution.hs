{- |
Module      : Dtmc.Distribution
Description : Dense finite and sparse finite-support probability distributions.

Probability distribution vectors over @n@ states, used as initial and marginal
laws of a DTMC. 'mkDistributionVector' checks the simplex invariant with the
@1e-9@ tolerance documented by 'SimplexError'. t'SparseDistribution' provides
the same validated invariant for finite support over an unrestricted state
type.
-}
module Dtmc.Distribution (
    Distribution (..),
    DistributionVector,
    DistributionError (..),
    mkDistributionVector,
    unDistributionVector,
    SparseDistribution,
    mkSparseDistribution,
    pointMass,
) where

import Data.Bifunctor (
    first,
 )
import Data.Finite (
    Finite,
    finites,
    getFinite,
 )
import Data.Map.Strict qualified as Map
import Dtmc.Distribution.Internal (
    DistributionVector (DistributionVector),
    SparseDistribution (SparseDistribution),
    unDistributionVector,
    unSparseDistribution,
 )
import Dtmc.Simplex (
    SimplexError,
 )
import Dtmc.Simplex.Internal (
    validateSimplex,
    validateSimplexEntries,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

{- | A simplex failure while constructing either distribution representation.
For a sparse law, coordinate indices refer to ascending state order after
duplicate states have been combined. The wrapper keeps distribution failures
distinct from transition-matrix row failures.
-}
newtype DistributionError
    = -- | Wrap the underlying simplex failure.
      DistributionError SimplexError
    deriving (Eq, Show)

{- | Validate a raw state distribution vector. On success it is preserved
exactly; tolerated coordinates are not clamped and the total is not
renormalised. For @n = 0@, returns
@Left (DistributionError (SumOffBy 0))@.

Time: @O(n)@. Space: @O(n)@ for validation.
-}
mkDistributionVector ::
    (KnownNat n) =>
    S.R n ->
    Either DistributionError (DistributionVector n)
mkDistributionVector vector =
    DistributionVector vector <$ first DistributionError (validateSimplex vector)

{- | Combine duplicate states, remove entries whose combined weight is exactly
zero, and validate the resulting finite-support probability law. Input order
is irrelevant; accepted weights are otherwise preserved without clamping or
renormalisation.

Time: @O(m log m)@ for @m@ supplied entries. Space: @O(s)@ for @s@ distinct
states.
-}
mkSparseDistribution ::
    (Ord state) =>
    [(state, Double)] ->
    Either DistributionError (SparseDistribution state)
mkSparseDistribution entries =
    SparseDistribution canonical
        <$ first DistributionError (validateSimplexEntries (Map.elems canonical))
  where
    canonical = Map.filter (/= 0) (Map.fromListWith (+) entries)

-- | The point mass concentrated on one state. This is valid by construction.
pointMass :: state -> SparseDistribution state
pointMass state = SparseDistribution (Map.singleton state 1)

{- | A representation of a discrete probability distribution with finite
stored support. A dense t'DistributionVector' and a t'SparseDistribution'
expose the same mathematical operations through this class.
-}
class Distribution distribution where
    -- | State type carried by the distribution representation.
    type DistributionState distribution

    {- | Read the stored probability of one state, returning exactly zero when
    the representation does not store that state. The value is returned
    without clamping or revalidation.
    -}
    probabilityAt ::
        (Ord (DistributionState distribution)) =>
        distribution ->
        DistributionState distribution ->
        Double

    {- | Return canonical ascending state weights. Exact-zero weights are
    omitted; tolerated negative values remain visible.
    -}
    distributionWeights ::
        distribution ->
        [(DistributionState distribution, Double)]

    {- | States with strictly positive stored weight, in ascending order.
    Tolerated negative coordinates are not mathematical support.
    -}
    support :: distribution -> [DistributionState distribution]
    support distribution =
        [ state
        | (state, weight) <- distributionWeights distribution
        , weight > 0
        ]

    {- | Convert to finite support without revalidation or renormalisation.
    Exact zero coordinates may be omitted.
    -}
    toSparseDistribution ::
        distribution ->
        SparseDistribution (DistributionState distribution)

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

    toSparseDistribution ::
        DistributionVector n ->
        SparseDistribution (DistributionState (DistributionVector n))
    toSparseDistribution distribution =
        SparseDistribution (Map.fromDistinctAscList (distributionWeights distribution))

instance Distribution (SparseDistribution state) where
    type DistributionState (SparseDistribution state) = state

    probabilityAt distribution state =
        Map.findWithDefault 0 state (unSparseDistribution distribution)

    distributionWeights = Map.toAscList . unSparseDistribution

    toSparseDistribution ::
        SparseDistribution state ->
        SparseDistribution
            (DistributionState (SparseDistribution state))
    toSparseDistribution = id
