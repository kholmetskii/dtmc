{- |
Module      : Dtmc.Distribution
Description : Dense finite and sparse finite-support probability distributions.

Probability distributions over @n@ states, used as initial and marginal laws
of a DTMC. 'mkDistribution' checks the simplex invariant with the @1e-9@
tolerance documented by 'SimplexError'. t'SparseDistribution' provides the same
validated invariant for finite support over an unrestricted state type.
-}
module Dtmc.Distribution (
    Distribution,
    DistributionError (..),
    mkDistribution,
    unDistribution,
    probabilityAt,
    SparseDistribution,
    SparseDistributionError (..),
    ToSparseDistribution (..),
    mkSparseDistribution,
    pointMass,
    sparseEntries,
    sparseProbabilityAt,
    sparseSupport,
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
    Distribution (Distribution),
    SparseDistribution (SparseDistribution),
    unDistribution,
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

{- | A simplex failure while constructing a distribution. The wrapper keeps it
distinct from a transition-matrix row failure.
-}
newtype DistributionError
    = -- | Wrap the underlying simplex failure.
      DistributionError SimplexError
    deriving (Eq, Show)

{- | Validate a raw state distribution. On success the vector is preserved
exactly; tolerated coordinates are not clamped and the total is not
renormalised. For @n = 0@, returns
@Left (DistributionError (SumOffBy 0))@.

Time: @O(n)@. Space: @O(n)@ for validation.
-}
mkDistribution :: (KnownNat n) => S.R n -> Either DistributionError (Distribution n)
mkDistribution vector =
    Distribution vector <$ first DistributionError (validateSimplex vector)

{- | Read the stored probability @mu(j)@ of a single state.

Mathematically this is the coordinate of the law at index @j@: for an initial
or marginal distribution @mu@, @probabilityAt mu j = mu(j) = P(X = j)@. The
bounded 'Finite' @n@ index makes the lookup total, so there is no
out-of-range or partial case.

The stored coordinate is returned exactly as held. No clamping to @[0, 1]@,
renormalisation, or revalidation is performed, so a value carried through
tolerated construction or floating-point arithmetic is reported unchanged and
may lie slightly outside @[0, 1]@.

Time: @O(n)@ to materialise the underlying vector; the index itself is @O(1)@.
-}
probabilityAt :: (KnownNat n) => Distribution n -> Finite n -> Double
probabilityAt distribution index =
    S.extract (unDistribution distribution)
        `LA.atIndex` fromIntegral (getFinite index)

{- | A simplex failure while constructing a sparse finite-support law. Entry
indices in the wrapped error refer to ascending state order after duplicate
states have been combined.
-}
newtype SparseDistributionError
    = -- | Wrap the shared simplex failure.
      SparseDistributionError SimplexError
    deriving (Eq, Show)

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
    Either SparseDistributionError (SparseDistribution state)
mkSparseDistribution entries =
    SparseDistribution canonical
        <$ first SparseDistributionError (validateSimplexEntries (Map.elems canonical))
  where
    canonical = Map.filter (/= 0) (Map.fromListWith (+) entries)

-- | The point mass concentrated on one state. This is valid by construction.
pointMass :: state -> SparseDistribution state
pointMass state = SparseDistribution (Map.singleton state 1)

{- | Initial-law representations accepted by shared sparse probability
algorithms. Dense finite distributions are converted once; sparse
distributions pass through unchanged.
-}
class ToSparseDistribution distribution where
    -- | State type carried by the distribution representation.
    type DistributionState distribution

    {- | Convert to finite support without revalidation or renormalisation.
    Exact zero coordinates may be omitted.
    -}
    toSparseDistribution ::
        distribution ->
        SparseDistribution (DistributionState distribution)

instance (KnownNat n) => ToSparseDistribution (Distribution n) where
    type DistributionState (Distribution n) = Finite n

    toSparseDistribution :: Distribution n -> SparseDistribution (DistributionState (Distribution n))
    toSparseDistribution distribution =
        SparseDistribution $
            Map.fromDistinctAscList
                [ (state, weight)
                | (state, weight) <- zip finites weights
                , weight /= 0
                ]
      where
        weights = LA.toList (S.extract (unDistribution distribution))

instance ToSparseDistribution (SparseDistribution state) where
    type DistributionState (SparseDistribution state) = state

    toSparseDistribution ::
        SparseDistribution state ->
        SparseDistribution
            (DistributionState (SparseDistribution state))
    toSparseDistribution = id

{- | Return the canonical ascending list of stored state weights. Exact zero
entries are absent; tolerated negative entries remain visible.

Time and result space: @O(s)@ for stored support size @s@.
-}
sparseEntries :: SparseDistribution state -> [(state, Double)]
sparseEntries = Map.toAscList . unSparseDistribution

{- | Read a state's stored probability, returning exactly zero when absent.
The stored value is returned without clamping or revalidation.

Time: @O(log s)@ for stored support size @s@.
-}
sparseProbabilityAt ::
    (Ord state) =>
    SparseDistribution state ->
    state ->
    Double
sparseProbabilityAt distribution state =
    Map.findWithDefault 0 state (unSparseDistribution distribution)

{- | States with strictly positive stored weight, in ascending order.
Tolerated negative coordinates are not members of mathematical support.

Time: @O(s)@; result space: at most @O(s)@.
-}
sparseSupport :: SparseDistribution state -> [state]
sparseSupport distribution =
    [state | (state, weight) <- sparseEntries distribution, weight > 0]
