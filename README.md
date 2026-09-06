# dtmc

`dtmc` is a Haskell library for discrete-time Markov chains. It provides
type-safe finite models, locally finite models over potentially infinite state
spaces, validated probability distributions, exact structural analysis, and
checked numerical analysis.

Use it to evolve distributions, calculate path probabilities and hitting
times, classify finite chains, find stationary and limiting distributions, or
simulate trajectories.

## Installation

Add the library to your Cabal package:

```cabal
build-depends: dtmc ^>=0.2.0.0
```

The package requires GHC 9.10 or later. It uses `hmatrix` internally, so a
BLAS/LAPACK implementation is required even when an application uses only
transition kernels.

On Ubuntu or Debian, install `libblas-dev` and `liblapack-dev`. On macOS,
`hmatrix` can use Apple Accelerate.

## Quick start

This example defines a finite weather chain, constructs it without unchecked
probabilities, and performs finite-time, hitting-time, and stationary
analyses.

```haskell
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}

module Main (main) where

import qualified Dtmc.Analysis.FiniteTime as FiniteTime
import qualified Dtmc.Analysis.HittingTime as HittingTime
import qualified Dtmc.Analysis.Stationary as Stationary
import Dtmc.Distribution (DistributionError)
import qualified Dtmc.Distribution.Map as DistributionMap
import Dtmc.Distribution.Vector (DistributionVector, DistributionVectorError)
import qualified Dtmc.Distribution.Vector as DistributionVector
import Dtmc.State (FiniteState)
import Dtmc.Transition.Kernel (TransitionKernel, fromLaws)
import Dtmc.Transition.Matrix (TransitionMatrix, fromKernel)
import GHC.Generics (Generic)

data Weather = Dry | Wet
  deriving (Eq, Ord, Show, Generic, FiniteState)

weatherKernel :: Either DistributionError (TransitionKernel Weather)
weatherKernel = do
  dryLaw <- DistributionMap.fromList [(Dry, 0.9), (Wet, 0.1)]
  wetLaw <- DistributionMap.fromList [(Dry, 0.4), (Wet, 0.6)]
  pure $
    fromLaws $ \state ->
      case state of
        Dry -> dryLaw
        Wet -> wetLaw

weatherMatrix :: Either DistributionError (TransitionMatrix Weather)
weatherMatrix = fromKernel <$> weatherKernel

initialWeather :: Either DistributionVectorError (DistributionVector Weather)
initialWeather = DistributionVector.fromList [1, 0]

main :: IO ()
main =
  case (weatherMatrix, initialWeather) of
    (Left problem, _) -> print problem
    (_, Left problem) -> print problem
    (Right matrix, Right initial) -> do
      -- P(X_2 = Wet | X_0 = Dry)
      print (FiniteTime.nStepProbability 2 matrix Dry Wet)

      -- P(eventually visit Wet) under the initial distribution
      print (HittingTime.eventualProbability matrix [Wet] initial)

      -- One extremal stationary distribution per recurrent class
      case Stationary.stationaryDistributions matrix of
        Left problem -> print problem
        Right classes ->
          print
            [ (states, DistributionVector.toList distribution)
            | (states, distribution) <- classes
            ]
```

Analysis modules are intended to be imported qualified because several of
them expose concise names such as `probability` and `expectation`.

## Choosing a representation

The state-space size and storage format are separate choices. A finite chain
can use either a functional kernel or a dense matrix. A potentially infinite
chain must use a locally finite kernel.

| State space | Transition representation | Transition laws | Initial distributions | Available operations |
|---|---|---|---|---|
| Finite | `TransitionKernel` | Sparse `DistributionMap` | `DistributionMap` or `DistributionVector` | Finite-horizon analysis and simulation |
| Finite | `TransitionMatrix` | Dense matrix rows | `DistributionMap` or `DistributionVector` | All analysis and simulation |
| Potentially infinite | `TransitionKernel` | Finite-support `DistributionMap` | Finite-support `DistributionMap` | Finite-horizon analysis and simulation |

Structural and infinite-horizon analyses—including classification,
absorption, stationarity, and limiting behaviour—require a finite
`TransitionMatrix`.

### Finite state types

Derive `FiniteState` for a fieldless enumeration with `Generic`, as in the
quick-start example. Constructor declaration order determines vector and
matrix order. A stock-derived `Ord` instance has the same order and is the
intended companion.

Use `Finite n` directly when named constructors are unnecessary. The library
also provides `FiniteState` instances for `()`, `Bool`, and `Ordering`.

### Potentially infinite state types

A `TransitionKernel` does not enumerate its state space. The following
examples define deterministic and stochastic chains over the integers:

```haskell
import Dtmc.Distribution (DistributionError)
import qualified Dtmc.Distribution.Map as DistributionMap
import Dtmc.Transition.Kernel
  ( TransitionKernel
  , fromLaws
  )

countUp :: TransitionKernel Integer
countUp = fromLaws (DistributionMap.pointMass . (+ 1))

randomWalk :: Either DistributionError (TransitionKernel Integer)
randomWalk = do
  steps <- DistributionMap.fromList
    [ (-1, 0.5)
    , (1, 0.5)
    ]

  pure $
    fromLaws $ \position ->
      DistributionMap.mapStates (+ position) steps
```

The random walk validates its relative step distribution once. `mapStates`
then translates that law to the current position while preserving its
probability mass. Every transition still has finite support, although the
states reachable over the whole lifetime of the chain are infinite.

Each transition law and every representable initial distribution has finite
support; no global enumeration or truncation is performed.

There is no `DistributionVector` or `TransitionMatrix` for an infinite state
space.

## Construction reference

These are the built-in public construction paths:

| Value | Recommended construction | Alternative construction |
|---|---|---|
| Sparse distribution | `Dtmc.Distribution.Map.fromList` or `pointMass` | `fromDistribution` |
| Dense finite distribution | `Dtmc.Distribution.Vector.fromList` | — |
| Functional transition kernel | `Dtmc.Transition.Kernel.fromLaws` | — |
| Dense finite transition matrix | `Dtmc.Transition.Matrix.fromRows` or `fromKernel` | — |

`Dtmc.Distribution.Map.mapStates` transforms a validated sparse distribution
and combines weights when several source states map to the same target.
`identity`, `compose`, and `power` construct new transition matrices from
existing ones.

Advanced users may define custom `Distribution` and `Transition` instances.
Generic finite-horizon analysis and simulation can use them when they satisfy
the documented finite-support contracts. Finite-chain structural analysis
still requires `TransitionMatrix`.

## Analysis guide

| Task | Module |
|---|---|
| Evolve a distribution | `Dtmc.Dynamics` |
| Calculate transition, joint, or conditional probabilities | `Dtmc.Analysis.FiniteTime` |
| Analyse hitting times or competing targets | `Dtmc.Analysis.HittingTime` |
| Analyse first-return times | `Dtmc.Analysis.ReturnTime` |
| Analyse finite or total visit counts | `Dtmc.Analysis.VisitCount` |
| Find communicating classes, periods, and recurrent states | `Dtmc.Analysis.Classification` |
| Calculate fundamental matrices and absorption quantities | `Dtmc.Analysis.Absorption` |
| Find extremal stationary distributions | `Dtmc.Analysis.Stationary` |
| Find ordinary or cyclic long-run limits | `Dtmc.Analysis.Limiting` |
| Sample states or simulate trajectories | `Dtmc.Simulation` |

Finite-horizon analysis works with both finite matrices and locally finite
kernels. Structural and infinite-horizon analysis currently requires a finite
`TransitionMatrix`.

## Validation and numerical behaviour

Probability-distribution smart constructors return `Either DistributionError`
instead of storing invalid values. They accept coordinate and total-mass error
within `1e-9`, clamp tolerated coordinate error to `[0, 1]`, and normalise the
repaired weights. `NaN` and infinite coordinates are rejected.

The explicit `hmatrix` transition-matrix constructor validates and repairs
each row under the same policy, returning `TransitionMatrixError` on failure.

Linear-system-based hitting, return, absorption, stationary, and limiting
analyses use checked `Double` calculations. Numerical failures are returned
as `LinearSystemError`; computed results are not silently clamped or
renormalised. Expectations that are mathematically infinite use
`InfiniteExpectation` rather than floating-point infinity.

Classification uses the exact positive support of the stored matrix: every
entry greater than zero is an edge, with no tolerance. Simulation validates
stored weights before consuming randomness and returns `SimulationError` on
failure.

See each function's Haddock documentation for its edge cases, numerical
contract, and time and space complexity.

## `hmatrix`

No `hmatrix` type appears anywhere in the public API. Matrices are built from
and read back as plain lists of weights, and distributions as state-labelled
weights or plain lists, so an application never mentions `hmatrix` in its own
signatures. The package does use `hmatrix` internally for the linear algebra,
so building it needs a BLAS/LAPACK implementation.

## Building and testing

```bash
cabal build all --enable-tests
cabal test all --test-show-details=direct
```

## Documentation and support

- The `Dtmc` module provides a complete module map.
- API documentation is available on
  [Hackage](https://hackage.haskell.org/package/dtmc/docs/Dtmc.html).
- Report problems through the
  [GitHub issue tracker](https://github.com/kholmetskii/dtmc/issues).

## License

`dtmc` is distributed under the BSD 3-Clause licence. See `LICENSE`.
