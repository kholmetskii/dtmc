-- |
-- Module      : Dtmc.Internal
--
-- Carriers and nothing else. No validation, no logic, no paired spec.
-- Listed under @other-modules@: invisible outside the library.
--
-- The raw constructors build a value BYPASSING validation. The rule:
--
--   * 'Dtmc.Distribution' and 'Dtmc.StochasticMatrix' are the validation
--     boundary; constructing is their job.
--   * Any OTHER module importing 'Dtmc.Internal' must apply the constructor
--     under a @-- Proof:@ block (see 'Dtmc.Kernel.rowAt').
--
-- Enforced by the @discipline@ job in CI, not by memory.
module Dtmc.Internal
  ( Distribution (..)
  , StochasticMatrix (..)
  ) where

import GHC.TypeNats (KnownNat, Nat)
import qualified Numeric.LinearAlgebra.Static as S

-- | A point of the standard simplex Δ^{n-1}.
newtype Distribution (n :: Nat) = Distribution
  { unDistribution :: S.R n
  }

-- | A row-stochastic transition matrix.
newtype StochasticMatrix (n :: Nat) = StochasticMatrix
  { unStochasticMatrix :: S.Sq n -- Sq n = L n n
  }

-- CRITICAL. See docs/DECISIONS.md D2.
--
-- Without these, GHC infers a phantom role for n, because in hmatrix
--
-- > newtype Dim (n :: Nat) t = Dim t          -- n does not occur on the right
-- > newtype R n   = R (Dim n (Vector ℝ))
-- > newtype L m n = L (Dim m (Dim n (Matrix ℝ)))
--
-- and `grep -rn "type role" hmatrix/` finds nothing. Coercible between the
-- same type constructor at a phantom parameter is solved by the lifting rule,
-- WITHOUT unwrapping the newtype, so a hidden constructor does not save us:
--
-- > coerce :: StochasticMatrix 2 -> StochasticMatrix 3   -- compiles without type role
--
-- A role annotation can only strengthen the inferred role. Regression test:
-- check/Role.hs must FAIL to compile.
type role Distribution nominal

type role StochasticMatrix nominal

-- hmatrix declares Show (R n) and Show (L m n) with a KnownNat context: you
-- need the size at runtime to print. Standalone deriving makes that
-- requirement visible in the source.
deriving instance KnownNat n => Show (Distribution n)

deriving instance KnownNat n => Show (StochasticMatrix n)

-- Eq is not derived and cannot be: hmatrix's R n and L m n have no Eq
-- instances. The rule "no structural equality on floating-point carriers" is
-- now held by the compiler rather than by discipline.
