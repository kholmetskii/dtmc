{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.TestSupport (
    genSimplexPoint,
    genTransitionMatrix,
    modifyMatrixRows,
    bumpSmallest,
    bumpSmallestInFirstRow,
    setFirstEntry,
) where

import Data.Proxy (
    Proxy (..),
 )
import GHC.TypeNats (
    KnownNat,
    natVal,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Test.QuickCheck (
    Gen,
    choose,
    frequency,
    vectorOf,
 )

genSimplexPoint :: Int -> Gen [Double]
genSimplexPoint dimension = do
    entries <- vectorOf dimension genEntry
    let total = sum entries
    if total == 0
        then genSimplexPoint dimension
        else pure (map (/ total) entries)
  where
    genEntry =
        frequency
            [ (3, pure 0)
            , (7, choose (0, 1000))
            ]

genTransitionMatrix ::
    forall n.
    (KnownNat n) =>
    Gen (S.Sq n)
genTransitionMatrix = do
    rows <- vectorOf dimension (genSimplexPoint dimension)
    pure (S.matrix (concat rows))
  where
    dimension = fromIntegral (natVal (Proxy @n))

modifyMatrixRows ::
    (KnownNat n) =>
    ([[Double]] -> [[Double]]) ->
    S.Sq n ->
    S.Sq n
modifyMatrixRows transform =
    S.matrix
        . concat
        . transform
        . LA.toLists
        . S.extract

bumpSmallest :: Double -> [Double] -> [Double]
bumpSmallest _ [] = []
bumpSmallest amount entries =
    zipWith bump [0 :: Int ..] entries
  where
    smallestIndex =
        snd (minimum (zip entries [0 :: Int ..]))

    bump index entry
        | index == smallestIndex = entry + amount
        | otherwise = entry

bumpSmallestInFirstRow ::
    Double ->
    [[Double]] ->
    [[Double]]
bumpSmallestInFirstRow _ [] = []
bumpSmallestInFirstRow amount (row : rows) =
    bumpSmallest amount row : rows

setFirstEntry ::
    Double ->
    [[Double]] ->
    [[Double]]
setFirstEntry value ((_ : rest) : rows) =
    (value : rest) : rows
setFirstEntry _ rows = rows
