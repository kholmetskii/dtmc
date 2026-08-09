{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.TestSupport (
    testTolerance,
    approxEq,
    approxDistributionEq,
    approxTransitionMatrixEq,
    genSimplexPoint,
    genTransitionMatrix,
    modifyMatrixRows,
    bumpSmallest,
    bumpSmallestInFirstRow,
    setFirstEntry,
    assignment1Matrix,
    assignment1Lambda,
) where

import Data.Proxy (
    Proxy (..),
 )
import Dtmc.Distribution.Vector (
    DistributionVector,
    mkDistributionVector,
    unDistributionVector,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    mkTransitionMatrix,
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

{- | Absolute slack the tests use when comparing floating-point results. Kept
independent of the library's private validation threshold so a change there
cannot silently mask a regression here; the two happen to share a value.
-}
testTolerance :: Double
testTolerance = 1e-9

{- | Absolute-tolerance comparison of two scalar 'Double' results, matching the
@abs (x - y) <= tolerance@ convention of the vector and matrix helpers.
-}
approxEq :: Double -> Double -> Double -> Bool
approxEq tolerance left right =
    abs (left - right) <= tolerance

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
    (KnownNat n) =>
    Double ->
    DistributionVector n ->
    DistributionVector n ->
    Bool
approxDistributionEq tolerance left right =
    and (zipWith close (entries left) (entries right))
  where
    entries = LA.toList . S.extract . unDistributionVector
    close x y = abs (x - y) <= tolerance

{- | Assignment 1 transition matrix over states @[A, B, C, D, E]@, shared by the
transition, dynamics, and probability specs.
-}
assignment1Matrix :: TransitionMatrix 5
assignment1Matrix =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 0
                , 0
                , 0
                , 1
                , 0
                , 1 / 3
                , 0
                , 0
                , 0
                , 2 / 3
                , 0
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 1 / 3
                , 2 / 3
                , 0
                , 1 / 4
                , 1 / 4
                , 0
                , 0
                , 1 / 2
                ] ::
                S.Sq 5
            )

-- | Assignment 1 initial law @lambda = [1/4, 1/2, 0, 1/4, 0]@.
assignment1Lambda :: DistributionVector 5
assignment1Lambda =
    either (error . show) id $
        mkDistributionVector (S.vector [1 / 4, 1 / 2, 0, 1 / 4, 0] :: S.R 5)
