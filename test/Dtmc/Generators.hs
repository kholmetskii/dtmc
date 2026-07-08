{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Generators
  ( genDenseStochasticMatrix2
  ) where

import Dtmc.StochasticMatrix
  ( StochasticMatrix
  , mkStochasticMatrix
  )
import Numeric.LinearAlgebra (fromLists)
import Test.QuickCheck
  ( Gen
  , choose
  , frequency
  , vectorOf
  )

-- | Generate a dense/sparse-ish stochastic matrix of size 2.
--
-- This generator currently normalises randomly generated rows and may include
-- zero entries because 'genEntry' can return 0.
--
-- Important: a truly Dirichlet(1,...,1)-based generator would produce strictly
-- positive rows. Such dense generators are useful for basic invariants, but
-- they do not test reducible, absorbing, or prescribed-zero-pattern chains.
-- Classification tests need separate structured generators.
genDenseStochasticMatrix2 :: Gen (StochasticMatrix 2)
genDenseStochasticMatrix2 = do
  rawRows <- vectorOf 2 (genNonZeroRow 2)
  let matrix = fromLists (map normalise rawRows)

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