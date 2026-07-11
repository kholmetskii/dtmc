module Dtmc.Internal.Types (
    Distribution (..),
    TransitionMatrix (..),
) where

import GHC.TypeNats (
    KnownNat,
    Nat,
 )
import Numeric.LinearAlgebra.Static qualified as S

newtype Distribution (n :: Nat) = Distribution
    { unDistribution :: S.R n
    }

type role Distribution nominal

deriving instance (KnownNat n) => Show (Distribution n)

newtype TransitionMatrix (n :: Nat) = TransitionMatrix
    { unTransitionMatrix :: S.Sq n
    }

type role TransitionMatrix nominal

deriving instance (KnownNat n) => Show (TransitionMatrix n)

instance (KnownNat n) => Semigroup (TransitionMatrix n) where
    (<>) :: TransitionMatrix n -> TransitionMatrix n -> TransitionMatrix n
    TransitionMatrix a <> TransitionMatrix b = TransitionMatrix{unTransitionMatrix = a S.<> b}

instance (KnownNat n) => Monoid (TransitionMatrix n) where
    mempty :: TransitionMatrix n
    mempty = TransitionMatrix{unTransitionMatrix = S.eye}
