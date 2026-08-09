# dtmc

A Haskell library for finite discrete-time Markov chains with type-safe
dimensions and locally finite countable-state chains with exact finite-horizon
algorithms.

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
- Validated finite-support laws and locally finite kernels over countable state
  types
- Exact sparse countable-state evolution, timed probability queries, bounded
  hitting/first-return probabilities, and finite simulation
- A shared `MarkovKernel` interface for applying the same sparse finite-horizon
  functions to finite matrices and locally finite infinite kernels

## Shared kernel interface

`MarkovKernel` captures exactly the operation shared by finite and infinite
chains: obtaining the validated finite-support law for one transition from a
given state. Both `TransitionMatrix n` and `TransitionKernel state` implement
it. The associated `KernelState` type keeps each kernel tied to its state type.
Likewise, `Distribution` is the common initial-law abstraction implemented by
`DistributionVector n` and `SparseDistribution state`.

Functions live in their mathematical subject modules and use `MarkovKernel`
where the finite and infinite signatures genuinely agree:

```haskell
import qualified Dtmc.Probability as Probability

-- matrix  :: TransitionMatrix n
-- initial :: DistributionVector n
-- state   :: Finite n

result =
  Probability.probabilityAtTime
    10
    initial
    matrix
    state
```

`Dtmc.Dynamics` keeps dense `evolve`/`evolveN` and explicitly named
`evolveSparse`/`evolveSparseN` because their result representations differ.
Probability, scalar bounded hitting/return, and simulation functions use the
shared kernel abstraction directly.

## Countable-state boundary

For an infinite state type, construct a `TransitionKernel` whose every row has
finite support. This guarantees that every finite-time calculation terminates,
although reachable support can still grow quickly. Shared functions perform no
state-space enumeration, truncation, clamping, or hidden approximation.

The shared API intentionally stops at finite horizons. It does not
offer eventual hitting probabilities, expected hitting/return times,
classification, or stationary distributions for an arbitrary infinite chain:
those questions need additional structure and will live in specialised
modules. This keeps the library a collection of DTMC algorithms rather than a
general probabilistic query interpreter.

For example, a simple random walk on all integers is locally finite:

```haskell
import qualified Dtmc.Distribution as Distribution
import qualified Dtmc.Hitting as Hitting
import qualified Dtmc.Kernel as Kernel

randomWalk :: Kernel.TransitionKernel Integer
randomWalk = Kernel.transitionKernel $ \i ->
  either (error . show) id $
    Distribution.mkSparseDistribution [(i - 1, 0.5), (i + 1, 0.5)]

-- P_0(H_{2} < 3) = 1/4
hitTwoBeforeThree :: Double
hitTwoBeforeThree =
  Hitting.hittingTimeProbabilityBefore randomWalk (== 2) 0 3
```

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

The library is in early development. The finite API covers type-safe objects,
single- and multi-step dynamics, support-graph classification, and exact-time,
bounded, eventual, and expected hitting and return quantities. The shared
kernel abstractions cover finite-time queries without tying algorithms to a
finite or infinite representation. Specialised infinite-chain solvers,
absorbing-chain summaries, and stationary/limiting behaviour are planned.
