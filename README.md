# dtmc

A small Haskell library for finite discrete-time Markov chains with type-safe dimensions.

## Current features

- Validated finite probability distributions
- Validated stochastic transition matrices
- Type-level vector and matrix dimensions using `hmatrix`
- Transition-matrix multiplication
- Access to the transition distribution from a given state
- Random sampling from a distribution
- Single-step Markov-chain simulation
- Approximate equality helpers for numerical comparisons

## Quick start

```haskell
{-# LANGUAGE DataKinds #-}

import Data.Finite (finite)
import Dtmc
import qualified Numeric.LinearAlgebra.Static as S
import qualified System.Random.MWC as MWC

main :: IO ()
main =
  case mkTransitionMatrix transitionMatrix of
    Left err ->
      print err

    Right matrix -> do
      generator <- MWC.createSystemRandom
      nextState <- step matrix (finite 0) generator
      print nextState
  where
    transitionMatrix :: S.Sq 2
    transitionMatrix =
      S.matrix
        [ 0.9, 0.1
        , 0.4, 0.6
        ]
```

## Building

The project requires GHC and Cabal.

```bash
cabal update
cabal build all --enable-tests
cabal test all --test-show-details=direct
```

`hmatrix` requires BLAS and LAPACK system libraries.

On Ubuntu or Debian:

```bash
sudo apt-get install libblas-dev liblapack-dev
```

On macOS, `hmatrix` can use Apple Accelerate.

## Status

The library is in early development. The current API covers the basic finite-state objects and single-step simulation. Multi-step dynamics and further Markov-chain analysis are planned.