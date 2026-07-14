-- |
-- Module      : Dtmc.Approx
-- Description : Floating-point tolerance policy and approximate equality.
--
-- The single home for the library's numeric slack. Distributions and transition
-- matrices are backed by 'Double's and deliberately have no 'Eq' instance
-- (bit-exact equality is the wrong notion for probabilities), so this module
-- provides the approximate comparisons users and tests actually want, all
-- referred to one shared 'tolerance'. It also hosts 'snapToSimplex', the
-- floating-point repair used before categorical sampling.
module Dtmc.Approx (
    tolerance,
    approxEq,
    approxEqR,
    approxEqDist,
    approxEqMatrix,
    snapToSimplex,
) where

import Dtmc.Internal.Types (
    Distribution,
    TransitionMatrix,
    unDistribution,
    unTransitionMatrix,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S

-- | Absolute slack allowed on every numeric comparison in the library
-- (simplex validation, approximate equality, and the snap in 'snapToSimplex').
-- Small enough to catch real modelling errors while tolerating rounding noise.
tolerance :: Double
tolerance = 1e-9

-- | Approximate equality of two scalars within 'tolerance'.
approxEq :: Double -> Double -> Bool
approxEq x y = abs (x - y) <= tolerance

-- | Coordinatewise approximate equality of two static vectors.
approxEqR :: (KnownNat n) => S.R n -> S.R n -> Bool
approxEqR a b =
    and (zipWith approxEq (entries a) (entries b))
  where
    entries = LA.toList . S.extract

-- | Approximate equality of two distributions, coordinatewise.
approxEqDist :: (KnownNat n) => Distribution n -> Distribution n -> Bool
approxEqDist a b = approxEqR (unDistribution a) (unDistribution b)

-- | Approximate equality of two transition matrices, entrywise.
approxEqMatrix :: (KnownNat n) => TransitionMatrix n -> TransitionMatrix n -> Bool
approxEqMatrix a b =
    and (zipWith approxEq (entries a) (entries b))
  where
    entries = LA.toList . LA.flatten . S.extract . unTransitionMatrix

-- | Snap small negative coordinates -- those within 'tolerance' of zero, an
-- artefact of floating-point arithmetic -- to exactly zero, so a probability
-- vector is accepted by a categorical sampler. A coordinate more negative than
-- that signals a real invariant violation (a programmer error) and fails loudly.
snapToSimplex :: LA.Vector Double -> LA.Vector Double
snapToSimplex =
    LA.cmap snap
  where
    snap value
        | value >= 0 = value
        | value >= negate tolerance = 0
        | otherwise =
            error
                ( "Dtmc.Approx.snapToSimplex: probability coordinate "
                    <> show value
                    <> " is below -tolerance"
                )
