-- | Generators. Per the M1 brief these are "core, not test scaffolding" — for
-- now they live under test/, and should move into a @dtmc-testing@ sublibrary
-- the first time anything outside the test suite needs them.
module Dtmc.Generators
  ( genSimplexPointList
  , genStochasticSq
  , onLists
  , bumpSmallest
  , bumpSmallestInRow0
  , setEntry00
  ) where

import Data.Proxy (Proxy (..))
import GHC.TypeNats (KnownNat, natVal)
import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Static as S
import Test.QuickCheck (Gen, choose, frequency, vectorOf)

-- | A point of Δ^{k-1} as a list.
--
-- This is NOT Dirichlet(1,…,1): 30% of coordinates are exactly zero. See
-- docs/TESTING.md T3 — a dense Dirichlet yields strictly positive rows, so
-- EVERY matrix would be regular (irreducible and aperiodic), and the M3
-- classification tests would be vacuously green.
--
-- Zeros do occur here, but they occur BY CHANCE, not by construction. Covering
-- reducible chains in M3 still requires structured generators with a prescribed
-- zero pattern. Do not rely on this generator there.
genSimplexPointList :: Int -> Gen [Double]
genSimplexPointList k = do
  xs <- vectorOf k genEntry
  let total = sum xs
  if total == 0.0
    then genSimplexPointList k -- probability 0.3^k; retrying is safe
    else pure (map (/ total) xs)
  where
    genEntry = frequency [(3, pure 0.0), (7, choose (0.0, 1000.0))]

-- | A row-stochastic n×n matrix. The size comes from the type.
genStochasticSq :: forall n. KnownNat n => Gen (S.Sq n)
genStochasticSq = do
  rows <- vectorOf k (genSimplexPointList k)
  pure (S.matrix (concat rows))
  where
    k = fromIntegral (natVal (Proxy @n))

-- | Perturbations are easier to describe on lists. @S.matrix@ throws on a
-- wrong-length list, so @f@ must preserve the shape.
onLists :: KnownNat n => ([[Double]] -> [[Double]]) -> S.Sq n -> S.Sq n
onLists f = S.matrix . concat . f . LA.toLists . S.extract

-- | Add @d@ to the SMALLEST coordinate.
--
-- Not the first one. See docs/TESTING.md T1: QuickCheck found the counterexample
-- [1.0, 0.0, 0.0] — bumping the first coordinate pushes it past 1 + ε, and the
-- coordinatewise check (which by contract runs BEFORE the sum) returns
-- EntryAboveOne rather than SumOffBy.
--
-- For a simplex point min ≤ 1/n, so after +1e-6 the smallest coordinate stays
-- safely below 1 + ε and the error is isolated to SumOffBy.
bumpSmallest :: Double -> [Double] -> [Double]
bumpSmallest _ [] = []
bumpSmallest d xs = zipWith bump [0 :: Int ..] xs
  where
    j = snd (minimum (zip xs [0 :: Int ..]))
    bump i x = if i == j then x + d else x

bumpSmallestInRow0 :: Double -> [[Double]] -> [[Double]]
bumpSmallestInRow0 _ [] = []
bumpSmallestInRow0 d (r : rs) = bumpSmallest d r : rs

-- | Replace entry (0,0). The row sum breaks too, but the coordinatewise check
-- runs first, so the reported error is NegativeEntry 0.
setEntry00 :: Double -> [[Double]] -> [[Double]]
setEntry00 v ((_ : rest) : rs) = (v : rest) : rs
setEntry00 _ m = m
