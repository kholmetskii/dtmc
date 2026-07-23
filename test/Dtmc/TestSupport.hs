{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.TestSupport (
    testTolerance,
    approxDistributionEq,
    approxTransitionMatrixEq,
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
import Dtmc.Distribution (
    Distribution,
    unDistribution,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    unTransitionMatrix,
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

-- | Absolute slack the tests use when comparing floating-point results. Kept
-- independent of the library's private validation threshold so a change there
-- cannot silently mask a regression here; the two happen to share a value.
testTolerance :: Double
testTolerance = 1e-9

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

approxTransitionMatrixEq ::
    (KnownNat n) => Double -> TransitionMatrix n -> TransitionMatrix n -> Bool
approxTransitionMatrixEq tolerance left right =
    and (zipWith close (entries left) (entries right))
  where
    entries = LA.toList . LA.flatten . S.extract . unTransitionMatrix
    close x y = abs (x - y) <= tolerance

approxDistributionEq ::
    (KnownNat n) => Double -> Distribution n -> Distribution n -> Bool
approxDistributionEq tolerance left right =
    and (zipWith close (entries left) (entries right))
  where
    entries = LA.toList . S.extract . unDistribution
    close x y = abs (x - y) <= tolerance
