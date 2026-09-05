{- |
Module      : Dtmc.Analysis.LinearSystem
Description : Numerical errors shared by finite-state linear-system analyses.

The explicit failure type shared by eventual hitting, return-time expectation,
and stationary-distribution calculations. The solvers themselves remain an
implementation detail; this module owns only their public error contract.
-}
module Dtmc.Analysis.LinearSystem (
    LinearSystemError (..),
) where

{- | Why a numerical linear-system result could not be accepted safely.

The reciprocal condition estimate and relative residual are dimensionless.
Smaller reciprocal condition estimates indicate greater sensitivity; smaller
relative residuals indicate a better computed solution.
-}
data LinearSystemError
    = -- | A required solve or decomposition has no usable unique result.
      SingularSystem
    | -- | The coefficient matrix is too sensitive for the numerical contract.
      IllConditionedSystem
        { reciprocalConditionEstimate :: Double
        -- ^ The backend's estimated reciprocal condition number.
        }
    | -- | A coefficient or right-hand-side entry was @NaN@ or infinite.
      NonFiniteSystem
    | -- | The solver produced a @NaN@ or infinite result.
      NonFiniteSolution
    | -- | The computed solution did not satisfy the equations closely enough.
      ResidualTooLarge
        { relativeResidual :: Double
        -- ^ The scaled residual of the computed solution.
        , residualLimit :: Double
        -- ^ The largest scaled residual accepted by the solver.
        }
    deriving (Eq, Show)
