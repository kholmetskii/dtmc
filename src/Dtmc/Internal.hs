-- |
-- Module      : Dtmc.Internal
--
-- Carriers and nothing else. No validation, no logic, no paired spec.
-- Listed under @other-modules@: invisible outside the library.
--
-- Haskell has no visibility modifier on a data constructor. Visibility is the
-- module export list, i.e. a file. So "the constructor is visible to a couple
-- of modules but not to the user" has exactly one encoding: a module kept out
-- of @exposed-modules@. This is @internal@ in C#, package-private in Java,
-- @pub(crate)@ in Rust. Precedent: @Data.Text.Internal@, and @Internal.Static@
-- inside hmatrix itself.
--
-- The raw constructors build a value BYPASSING validation. The rule:
--
--   * 'Dtmc.Distribution' is the validation boundary for Δ^{n-1}; constructing
--     is its job.
--   * Any other module importing 'Dtmc.Internal' must apply a raw constructor
--     under a @-- Proof:@ block. Today that is 'Dtmc.TransitionMatrix'
--     ('Dtmc.TransitionMatrix.mulTransitionMatrix', 'Dtmc.TransitionMatrix.rowAt');
--     M2 will add @Dtmc.Dynamics@.
--
-- The list of importers of this module is the library's audit list:
--
-- > grep -rl '^import Dtmc.Internal' src/
--
-- Enforced by the @discipline@ job in CI, not by memory.
module Dtmc.Internal
  ( Distribution (..)
  , TransitionMatrix (..)
  ) where

import GHC.TypeNats (KnownNat, Nat)
import qualified Numeric.LinearAlgebra.Static as S

-- | A point of the standard simplex Δ^{n-1}.
--
-- A row of a transition matrix IS a distribution — not "resembles one".
-- There is therefore no separate @StochasticRow@ type: it would be a second
-- name for the same set. See docs/DECISIONS.md D9.
newtype Distribution (n :: Nat) = Distribution
  { unDistribution :: S.R n
  }

-- | A row-stochastic transition matrix.
--
-- The representation is a dense @Sq n@ — a pointer to a contiguous buffer that
-- LAPACK multiplies in one BLAS call. It does NOT hold n 'Distribution' values;
-- it holds n² doubles. A row is not stored, it is carved out
-- ('Dtmc.TransitionMatrix.rowAt').
newtype TransitionMatrix (n :: Nat) = TransitionMatrix
  { unTransitionMatrix :: S.Sq n -- Sq n = L n n
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
-- > coerce :: TransitionMatrix 2 -> TransitionMatrix 3   -- compiles without type role
--
-- A role annotation can only strengthen the inferred role. Regression test:
-- check/Role.hs must FAIL to compile.
type role Distribution nominal

type role TransitionMatrix nominal

-- hmatrix declares Show (R n) and Show (L m n) with a KnownNat context: you
-- need the size at runtime to print. Standalone deriving makes that
-- requirement visible in the source.
deriving instance KnownNat n => Show (Distribution n)

deriving instance KnownNat n => Show (TransitionMatrix n)

-- Eq is not derived and cannot be: hmatrix's R n and L m n have no Eq
-- instances. The rule "no structural equality on floating-point carriers" is
-- now held by the compiler rather than by discipline.
