# dtmc

A Haskell library for finite discrete-time Markov chains with type-safe
dimensions and locally finite countable-state chains with exact finite-horizon
algorithms.

## Current features

- Validated dense and finite-support probability distributions
- Validated stochastic transition matrices
- Compile-time vector and matrix dimensions derived from finite state types
- Transition-matrix multiplication
- Access to the transition distribution from a given state
- Random sampling from a distribution
- Single-step Markov-chain simulation
- Distribution evolution, matrix powers, and the Chapman–Kolmogorov law
- Scalar, path, timed-event, and conditional probability queries
- Support-graph classification: accessibility, communication, communicating
  classes, irreducibility, and per-class periods (exact, combinatorial), with
  an `Irreducible` witness type
- Exact-time, bounded, eventual, competing, and expected hitting and
  first-return quantities
- Exact finite-horizon visit-count distributions and infinite-horizon total
  visit-count probabilities and expectations
- Unique stationary distributions for finite irreducible chains
- Validated finite-support laws and locally finite kernels over countable state
  types
- Exact sparse countable-state evolution, timed probability queries, bounded
  hitting/first-return probabilities, and finite simulation
- A shared `Transition` interface for applying the same sparse finite-horizon
  functions to finite matrices and locally finite infinite kernels

## Shared transition interface

`Transition` captures exactly the operation shared by finite and infinite
chains: obtaining the validated finite-support law for one transition from a
given state. Both `TransitionMatrix state` and `TransitionKernel state` implement
it. The associated `TransitionState` type keeps each transition representation
tied to its state type.
Likewise, `Distribution` is the common initial-law abstraction implemented by
`DistributionVector state` and `DistributionMap state`.

## Module layout

The abstract capabilities and concrete representations are separated:

```text
Dtmc                              complete public facade

Dtmc.State
Dtmc.Simplex

Dtmc.Distribution
├── Dtmc.Distribution.Vector
└── Dtmc.Distribution.Map

Dtmc.Transition
├── Dtmc.Transition.Matrix
└── Dtmc.Transition.Kernel

Dtmc.Dynamics
Dtmc.Simulation

Dtmc.Analysis.*                    focused analysis namespace
├── Dtmc.Analysis.FiniteTime
├── Dtmc.Analysis.LinearSystem
├── Dtmc.Analysis.HittingTime
├── Dtmc.Analysis.ReturnTime
├── Dtmc.Analysis.VisitCount
├── Dtmc.Analysis.Classification
└── Dtmc.Analysis.Stationary
```

Import `Dtmc` for the curated complete API. Use the focused modules when a
library component should depend only on an abstraction or one representation.
There is intentionally no second `Dtmc.Analysis` facade: the
`Dtmc.Analysis.*` names group focused mathematical subjects without duplicating
the top-level export surface.

Functions live in their mathematical subject modules and use `Transition`
where the finite and infinite signatures genuinely agree:

```haskell
import qualified Dtmc.Analysis.FiniteTime as FiniteTime

-- matrix  :: TransitionMatrix Weather
-- initial :: DistributionVector Weather
-- state   :: Weather

result =
  FiniteTime.probabilityAtTime
    10
    initial
    matrix
    state
```

`Dtmc.Dynamics` exposes shared `evolve`/`evolveN` operations that always return
a `DistributionMap`. The optimized `evolveVector`/`evolveVectorN` operations
preserve `DistributionVector` when both inputs use the dense finite
representation. `sample` accepts either distribution representation directly.
Probability, scalar bounded hitting/return, and simulation functions use the
shared abstractions directly.

## Analysis API

Finite-time functions work with any `Distribution` and `Transition`, including
locally finite kernels over countably infinite state types:

| Function | Quantity |
| --- | --- |
| `transitionProbability` | One-step transition probability |
| `transitionProbabilityN` | Transition probability after exactly `n` steps |
| `probabilityAtTime` | Marginal state probability at time `n` |
| `pathProbability` | Probability of a consecutive non-empty path |
| `jointProbability` | Joint probability of timed observations |
| `conditionalProbability` | Conditional probability of timed observations |

Singular function names return the result for one initial state. Plural names
return results for every finite state in canonical state order.

| Quantity | All states | One state |
| --- | --- | --- |
| First hit at exactly time `n` | `hittingTimeProbabilities` | `hittingTimeProbability` |
| First hit strictly before time `n` | `hittingTimeBeforeProbabilities` | `hittingTimeBeforeProbability` |
| Eventual hit | `hittingProbabilities` | `hittingProbability` |
| Hit one set before another | `hittingBeforeProbabilities` | `hittingBeforeProbability` |
| Expected hitting time | `hittingTimeExpectations` | `hittingTimeExpectation` |
| First return at exactly time `n` | `returnTimeProbabilities` | `returnTimeProbability` |
| First return strictly before time `n` | `returnTimeBeforeProbabilities` | `returnTimeBeforeProbability` |
| Eventual return | `returnProbabilities` | `returnProbability` |
| Expected return time | `returnTimeExpectations` | `returnTimeExpectation` |
| Total visit-count outcome | `visitCountProbabilities` | `visitCountProbability` |
| Expected total visits | `visitCountExpectations` | `visitCountExpectation` |

`hittingTimeBeforeProbability n` is a finite-time query for
`P(H_A < n)`. In contrast, `hittingBeforeProbability` is the competing-event
query `P(H_A < H_B)` and has no time-bound argument.

Finite-horizon visit-count analysis starts from an initial distribution and
uses a predicate to identify the visited set:

| Function | Result before time `n` |
| --- | --- |
| `visitCountDistributionBefore` | Complete distribution of the visit count |
| `visitCountProbabilityBefore` | Probability of exactly a supplied number of visits |
| `visitCountExpectationBefore` | Expected number of visits |

Every `Before` bound is strict. Visit counts before `n` inspect times
`0, ..., n - 1`, so the initial state is counted when `n > 0`.

Infinite-horizon analysis instead fixes one target state `i` and studies the
total count
`V_i = sum (t = 0 .. infinity) 1_{X_t = i}` for each possible initial state.
It also includes time zero. Query an individual outcome with:

```haskell
visitCountProbability matrix target (FiniteVisits 3) initial
visitCountProbability matrix target InfiniteVisits initial
visitCountExpectation matrix target initial
```

`VisitCountOutcome` distinguishes `FiniteVisits n` from `InfiniteVisits`, and
`MeanCount` distinguishes `FiniteMeanCount x` from `InfiniteMeanCount`. The
library intentionally exposes outcome probabilities rather than a finite map
called a “distribution”: a transient target can have infinitely many positive
finite-count probabilities, while a recurrent target can carry mass at
infinity. These finite-matrix calculations use hitting probabilities, return
probabilities, and exact support-graph classification; they do not simulate or
truncate the path.

`Dtmc.Analysis.LinearSystem` owns the `LinearSystemError` contract shared by
eventual hitting and return, total visit-count, and stationary-distribution
calculations. Numerical solver functions remain internal.

## Countable-state boundary

For an infinite state type, construct a `TransitionKernel` whose every row has
finite support. This guarantees that every finite-time calculation terminates,
although reachable support can still grow quickly. Shared functions perform no
state-space enumeration, truncation, clamping, or hidden approximation.

For arbitrary infinite kernels, the shared API intentionally stops at finite
horizons. It does not offer eventual hitting probabilities, expected
hitting/return times, infinite-horizon total visits, classification, or
stationary distributions without additional structure. Those calculations
currently require a finite `TransitionMatrix`. This keeps the library a
collection of DTMC algorithms rather than a general probabilistic query
interpreter.

For example, a simple random walk on all integers is locally finite:

```haskell
import qualified Dtmc.Distribution.Map as DistributionMap
import qualified Dtmc.Analysis.HittingTime as HittingTime
import qualified Dtmc.Transition.Kernel as Kernel

randomWalk :: Kernel.TransitionKernel Integer
randomWalk = Kernel.transitionKernel $ \i ->
  either (error . show) id $
    DistributionMap.mkDistributionMap [(i - 1, 0.5), (i + 1, 0.5)]

-- P_0(H_{2} < 3) = 1/4
hitTwoBeforeThree :: Double
hitTwoBeforeThree =
  HittingTime.hittingTimeBeforeProbability 3 randomWalk (== 2) 0
```

## Numerical contract

The public constructors validate probability distributions and transition
matrices. Analysis functions then assume those invariants and use ordinary
`Double` arithmetic: computed values are not clamped, snapped, or
renormalised. Support-graph algorithms treat a stored entry as an edge exactly
when it is strictly positive.

Construction and mathematically undefined query results use explicit error
values. Eventual and expected hitting/return analyses and finite stationary
distribution calculations return
`Either LinearSystemError result`. The shared solver rejects non-finite systems or
solutions, reciprocal condition estimates below `1e-12`, and scaled residuals
above `1e-9`; no numerical failure is converted into a runtime exception.

## Quick start

```haskell
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}

import Dtmc
import GHC.Generics (Generic)
import qualified Numeric.LinearAlgebra.Static as S
import qualified System.Random.MWC as MWC

data Weather = Dry | Wet
  deriving (Eq, Ord, Show, Generic, FiniteState)

weatherMatrix :: S.Sq 2
weatherMatrix =
  S.matrix
    [ 0.9, 0.1
    , 0.4, 0.6
    ]

main :: IO ()
main =
  case mkTransitionMatrix weatherMatrix :: Either TransitionMatrixError (TransitionMatrix Weather) of
    Left err ->
      print err

    Right matrix -> do
      generator <- MWC.createSystemRandom
      nextState <- step matrix Dry generator
      print nextState
```

Constructor declaration order defines the dense vector and matrix coordinate
order. `Finite n` remains available as the low-level indexed state type when
named constructors are not useful.

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
bounded, eventual, and expected hitting and return quantities, finite-horizon
visit counts, infinite-horizon total visit-count probabilities and
expectations, and stationary distributions of finite irreducible chains. The
shared kernel abstractions cover finite-time queries without tying algorithms
to a finite or infinite representation. Specialised infinite-chain solvers,
absorbing-chain summaries, and limiting behaviour are planned.
