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
- Distribution evolution, matrix powers, and the Chapman–Kolmogorov law
- Scalar, path, timed-event, and conditional probability queries
- Support-graph classification: accessibility, communication, communicating
  classes, irreducibility, and per-class periods (exact, combinatorial), with
  an `Irreducible` witness type
- Exact-time, bounded, eventual, competing, and expected hitting and
  first-return quantities

## Numerical contract

The public constructors validate probability distributions and transition
matrices. Analysis functions then assume those invariants and use ordinary
`Double` arithmetic: computed values are not clamped, snapped, or
renormalised. Support-graph algorithms treat a stored entry as an edge exactly
when it is strictly positive.

Construction and mathematically undefined query results use explicit error
values. The current `0.x` hitting routines can still raise an exception if the
numerical backend rejects a linear system that is nonsingular in exact
arithmetic; their module documentation identifies those cases. New numerical
analyses must return an explicit error value instead, and the existing solver
failure path will be migrated before `1.0`.

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

The library is in early development. The current API covers the finite-state objects, single- and multi-step dynamics (Chapman–Kolmogorov), the combinatorial structure theory of the support graph (communicating classes, irreducibility, periodicity, and recurrence/transience), and exact-time, bounded, eventual, and expected hitting and return quantities. Absorbing chains via the fundamental matrix, and stationary/limiting behaviour, are planned.
