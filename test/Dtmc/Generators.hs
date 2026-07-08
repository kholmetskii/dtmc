{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Generators
  ( genDenseStochasticRawMatrix2
  , genDenseStochasticMatrix2
  ) where

import Dtmc.StochasticMatrix
  ( StochasticMatrix
  , mkStochasticMatrix
  )
import Numeric.LinearAlgebra
  ( Matrix
  , fromLists
  )
import Test.QuickCheck
  ( Gen
  , choose
  , frequency
  , vectorOf
  )

-- | Generate a raw stochastic matrix of size 2.
--
-- The rows are generated as non-negative vectors and then normalised in
-- 'Double'. Therefore row sums may differ from @1@ by ordinary floating-point
-- normalisation error.
--
-- The validation tolerance in 'mkStochasticMatrix' must be larger than this
-- normalisation error. For small dimensions, this error is on the order of
-- @n * u@, where @u@ is double-precision machine epsilon, and the library's
-- tolerance @1e-9@ is safely larger.
--
-- This generator is appropriate for constructor round-trip and multiplication
-- closure properties.
--
-- Warning: dense Dirichlet-style generators produce strictly positive rows,
-- hence typically only regular/ergodic chains. Such generators are not suitable
-- for testing reducible chains, absorbing chains, communicating classes, or
-- prescribed zero patterns. Those require separate structured generators.
--
-- This current generator is not a true Dirichlet generator: it intentionally
-- allows zero entries, so it can also produce sparse matrices.
genDenseStochasticRawMatrix2 :: Gen (Matrix Double)
genDenseStochasticRawMatrix2 = do
  rawRows <- vectorOf 2 (genNonZeroRow 2)
  pure (fromLists (map normalise rawRows))

-- | Generate a validated stochastic matrix of size 2.
--
-- This should always succeed if the generator's floating-point normalisation
-- error is below the validation tolerance.
genDenseStochasticMatrix2 :: Gen (StochasticMatrix 2)
genDenseStochasticMatrix2 = do
  matrix <- genDenseStochasticRawMatrix2
  case mkStochasticMatrix @2 matrix of
    Right stochasticMatrix ->
      pure stochasticMatrix
    Left err ->
      error ("genDenseStochasticMatrix2 produced invalid matrix: " <> show err)

genNonZeroRow :: Int -> Gen [Double]
genNonZeroRow n = do
  row <- vectorOf n genEntry
  if sum row == 0.0
    then genNonZeroRow n
    else pure row

genEntry :: Gen Double
genEntry =
  frequency
    [ (3, pure 0.0)
    , (7, choose (0.0, 1000.0))
    ]

normalise :: [Double] -> [Double]
normalise row =
  let rowTotal = sum row
   in map (/ rowTotal) row