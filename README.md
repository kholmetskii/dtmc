# dtmc

A Haskell library for discrete-time Markov chains with type-safe finite models
and locally finite countable-state models.

Add the library to your Cabal package:

```cabal
build-depends: dtmc ^>=0.2.0.0
```

## Chain representations

### Finite chains

Finite chains use a `FiniteState` type and a validated `TransitionMatrix`.
State cardinality is known at compile time, so dense vectors and matrices have
type-safe dimensions.

Initial distributions can use either:

- `DistributionVector` for a dense, statically sized representation;
- `DistributionMap` for a sparse finite-support representation.

Distribution smart constructors accept floating-point error within `1e-9`,
then clamp tolerated coordinate error and normalise before storing the value.
The explicit `hmatrix` matrix constructor applies the same policy row by row.
`NaN` and infinite coordinates are rejected explicitly.

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
- occupation and fundamental matrices;
- canonical decomposition and absorption probabilities and expectations;
- accessibility, communicating classes, recurrence, and periodicity;
- extremal stationary distributions for every recurrent class;
- ordinary limiting matrices and cyclic subsequential limits;
- random sampling and finite simulation.

Finite-horizon analysis works with both chain representations.
Infinite-horizon analysis and structural classification currently require a
finite `TransitionMatrix`.

Simulation functions return `Either SimulationError ...`. Invalid unchecked
weights are reported without consuming randomness; only successfully validated
sampling passes the supplied generator to the categorical backend.

## Example

```haskell
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}

import Dtmc.Analysis.FiniteTime qualified as FT
import Dtmc.Distribution.Map (DistributionMap, mkDistributionMap)
import Dtmc.State (FiniteState)
import Dtmc.Transition.Kernel (transitionKernel)
import Dtmc.Transition.Matrix (TransitionMatrix, fromKernel)
import GHC.Generics (Generic)

data Weather = Dry | Wet
  deriving (Eq, Ord, Show, Generic, FiniteState)

law :: [(Weather, Double)] -> DistributionMap Weather
law = either (error . show) id . mkDistributionMap

weather :: TransitionMatrix Weather
weather =
  fromKernel $
    transitionKernel $ \state ->
      case state of
        Dry -> law [(Dry, 0.9), (Wet, 0.1)]
        Wet -> law [(Dry, 0.4), (Wet, 0.6)]

wetTomorrow :: Double
wetTomorrow = FT.stepProbability weather Dry Wet
```

See the `Dtmc` module on Hackage for a map of the construction, simulation,
and analysis modules. Analysis modules are intended to be imported qualified
because several subjects expose concise names such as `probability` and
`expectation`.

## hmatrix interoperability

The ordinary construction and inspection APIs use state-labelled weights and
plain lists. Projects already using `hmatrix` can opt into
`Dtmc.Distribution.Vector.HMatrix` and `Dtmc.Transition.Matrix.HMatrix` for
statically sized construction and zero-copy projection.

`dtmc` still uses `hmatrix` internally, so building the package requires a
BLAS/LAPACK implementation even when an application only uses kernels.

## Building

The project requires GHC 9.10 or later, Cabal, and BLAS/LAPACK for `hmatrix`.

```bash
cabal build all --enable-tests
cabal test all --test-show-details=direct
```

On Ubuntu or Debian, install `libblas-dev` and `liblapack-dev`. On macOS,
`hmatrix` can use Apple Accelerate.
