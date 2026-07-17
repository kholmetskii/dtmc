-- |
-- Module      : Dtmc.Internal.Block
-- Description : Sub-block extraction and @(I - Q)@ linear solves.
--
-- The numeric substrate for the hitting-time and absorption theory: extract
-- the block of a matrix indexed by chosen subsets of states, and solve the
-- linear systems of the shape @(I - Q) x = b@ that first-step decompositions
-- produce (hitting probabilities, expected hitting times, the fundamental
-- matrix).
--
-- Everything here is deliberately /dynamic/ ("Numeric.LinearAlgebra" rather
-- than the statically sized types): the sizes of the blocks -- how many
-- states are transient, how many lie outside a target set -- are runtime
-- values that the type-level dimension @n@ cannot express without existential
-- plumbing. Public modules translate 'Data.Finite.Finite' indices to raw
-- 'Int's on the way in and re-index the results on the way out, so no
-- dynamically sized value leaks into the public API.
--
-- Like "Dtmc.Internal.Graph", this module knows nothing about probabilities:
-- it is plain linear algebra, and the Markov-chain meaning of each system
-- lives in the caller. It is hidden like the other internal modules and, like
-- them, is verified through the public API: the hitting-time and absorption
-- specs exercise every function here against closed forms and simulations.
module Dtmc.Internal.Block (
    subMatrix,
    rowSums,
    solveIminusQ,
    solveIminusQVector,
    fundamental,
) where

import Numeric.LinearAlgebra qualified as LA

-- | The block of @m@ picked out by the given row and column indices: entry
-- @(a, b)@ of the result is @m (rows !! a, cols !! b)@. Order (and, if
-- present, multiplicity) of the index lists is preserved. Precondition:
-- indices in bounds and both lists non-empty -- callers short-circuit the
-- empty case (an empty transient set, a target covering every state) before
-- building a system.
subMatrix :: [Int] -> [Int] -> LA.Matrix Double -> LA.Matrix Double
subMatrix rowIdx colIdx m =
    m LA.?? (LA.Pos (LA.idxs rowIdx), LA.Pos (LA.idxs colIdx))

-- | The vector of row sums of @m@, i.e. @m@ applied to a vector of ones. For
-- a block @P[D, A]@ of a transition matrix this is the one-step probability
-- of stepping from each state of @D@ straight into @A@.
rowSums :: LA.Matrix Double -> LA.Vector Double
rowSums m = m LA.#> LA.konst 1 (LA.cols m)

-- | Solve @(I - Q) X = B@ by LU decomposition, one column of @X@ per column
-- of @B@. @Nothing@ when the factorisation finds the system singular. For
-- every system this library builds, the theory guarantees @I - Q@ invertible
-- (the spectral radius of @Q@ is below one), so @Nothing@ signals a caller
-- bug -- typically a state wrongly included in the solve set -- rather than
-- a property of the chain.
solveIminusQ :: LA.Matrix Double -> LA.Matrix Double -> Maybe (LA.Matrix Double)
solveIminusQ q =
    LA.linearSolve (LA.ident (LA.rows q) - q)

-- | 'solveIminusQ' specialised to a single right-hand side.
solveIminusQVector :: LA.Matrix Double -> LA.Vector Double -> Maybe (LA.Vector Double)
solveIminusQVector q b =
    LA.flatten <$> solveIminusQ q (LA.asColumn b)

-- | The fundamental matrix of @Q@: @(I - Q)^-1 = sum_k Q^k@, computed as the
-- solve @(I - Q) G = I@ rather than by explicit inversion. The Neumann-series
-- identity is the caller's theorem to state and prove; here it is only the
-- reason the name is apt.
fundamental :: LA.Matrix Double -> Maybe (LA.Matrix Double)
fundamental q =
    solveIminusQ q (LA.ident (LA.rows q))
