module Dtmc.Internal.Types
  ( Distribution (..)
  , TransitionMatrix (..)
  ) where

import GHC.TypeNats 
  (KnownNat
  , Nat
  )
import qualified Numeric.LinearAlgebra.Static as S

newtype Distribution (n :: Nat) = Distribution
  { unDistribution :: S.R n
  }

type role Distribution nominal

deriving instance KnownNat n => Show (Distribution n)


newtype TransitionMatrix (n :: Nat) = TransitionMatrix
  { unTransitionMatrix :: S.Sq n
  }

type role TransitionMatrix nominal

deriving instance KnownNat n => Show (TransitionMatrix n)

