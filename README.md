# dtmc

A Haskell library for discrete-time Markov chains with type-safe finite models
and locally finite countable-state models.

## Chain representations

### Finite chains

Finite chains use a `FiniteState` type and a validated `TransitionMatrix`.
State cardinality is known at compile time, so dense vectors and matrices have
type-safe dimensions.

Initial distributions can use either:

- `DistributionVector` for a dense, statically sized representation;
- `DistributionMap` for a sparse finite-support representation.

### Countable-state chains

Countable-state chains use `TransitionKernel`. A kernel may have an infinite
state type, but the transition law from each state must have finite support.
Distributions use the sparse `DistributionMap` representation.

`Distribution` and `Transition` provide shared interfaces, allowing the same
finite-horizon algorithms to work with both finite matrices and locally finite
kernels.

## Analysis

The library supports:

- distribution evolution and multi-step transition probabilities;
- path, joint, and conditional probabilities;
- exact-time and bounded hitting and first-return probabilities;
- eventual and competing hitting probabilities for finite chains;
- expected hitting and return times for finite chains;
- finite-horizon and infinite-horizon visit-count analysis;
- accessibility, communicating classes, recurrence, and periodicity;
- stationary distributions for finite irreducible chains;
- random sampling and finite simulation.

Finite-horizon analysis works with both chain representations.
Infinite-horizon analysis and structural classification currently require a
finite `TransitionMatrix`.

## Example

```haskell
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}

import Dtmc.Analysis.FiniteTime qualified as FT
import Dtmc.State (FiniteState)
import Dtmc.Transition.Matrix (TransitionMatrix, mkTransitionMatrix)
import GHC.Generics (Generic)
import qualified Numeric.LinearAlgebra.Static as S

data Weather = Dry | Wet
  deriving (Eq, Ord, Show, Generic, FiniteState)

weather :: TransitionMatrix Weather
weather =
  either (error . show) id $
    mkTransitionMatrix
      (S.matrix [0.9, 0.1, 0.4, 0.6] :: S.Sq 2)

wetTomorrow :: Double
wetTomorrow = FT.stepProbability weather Dry Wet
```

## Building

The project requires GHC, Cabal, and BLAS/LAPACK for `hmatrix`.

```bash
cabal build all --enable-tests
cabal test all --test-show-details=direct
```

On Ubuntu or Debian, install `libblas-dev` and `liblapack-dev`. On macOS,
`hmatrix` can use Apple Accelerate.
