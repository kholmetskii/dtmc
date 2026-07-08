# Design Ledger

This file records design decisions that affect the public API, mathematical invariants, numerical assumptions, and testing strategy.

The goal of this ledger is not to document every implementation detail. It records decisions that future contributors should not accidentally undo.

---

## M1 — Core stochastic guarantee

### Decision

The core validated carriers are:

```haskell
ProbabilityVector n
StochasticMatrix n
```

Both types are abstract. Their constructors are not exported.

Values can only be created through smart constructors:

```haskell
mkProbabilityVector
mkStochasticMatrix
```

The intended invariant is:

```text
ProbabilityVector n
  = a vector of length n whose entries are numerically probabilities

StochasticMatrix n
  = an n-by-n matrix whose rows are ProbabilityVector n values
```

### Reason

A row-stochastic matrix is mathematically a square matrix whose rows are probability distributions.

Therefore, matrix validation should not duplicate probability-vector validation. Instead, `mkStochasticMatrix` extracts each row and delegates row validation to `mkProbabilityVectorAt`.

This gives a single source of truth for the statement:

```text
a row is a probability distribution
```

### Consequence

A value of type:

```haskell
StochasticMatrix n
```

means:

```text
This matrix was checked at a constructor boundary to be n-by-n and numerically row-stochastic.
```

It does not mean that the rows sum to exactly `1` in the mathematical real-number sense.

---

## M2 — Abstract constructors

### Decision

The constructors for `ProbabilityVector n` and `StochasticMatrix n` are not exported.

The modules should export the type names and accessors, but not the constructors.

Good:

```haskell
module Dtmc.ProbabilityVector
  ( ProbabilityVector
  , unProbabilityVector
  , mkProbabilityVector
  , mkProbabilityVectorAt
  , probabilityTolerance
  ) where
```

Good:

```haskell
module Dtmc.StochasticMatrix
  ( StochasticMatrix
  , unStochasticMatrix
  , mkStochasticMatrix
  , mulStochasticMatrix
  , approxStochasticMatrixEq
  ) where
```

Bad:

```haskell
module Dtmc.StochasticMatrix
  ( StochasticMatrix(..)
  ) where
```

### Reason

If the constructor is exported, users can bypass validation:

```haskell
StochasticMatrix rawInvalidMatrix
```

That would destroy the invariant.

### Consequence

The phrase “invalid stochastic matrices are unrepresentable” is only true because the constructors are hidden and all construction goes through smart constructors.

---

## N1 — Numerical stochasticity

### Decision

The library uses `Double` and numerical validation with a fixed tolerance:

```haskell
probabilityTolerance :: Double
probabilityTolerance = 1e-9
```

Probability-vector validation checks:

```haskell
x >= -probabilityTolerance
x <= 1.0 + probabilityTolerance
abs (sum row - 1.0) <= probabilityTolerance
```

### Reason

The library is intended to support spectral computations such as stationary distributions, power iteration, eigensolvers, and Perron-Frobenius-style calculations.

For a general stochastic matrix, eigenvalues and eigenvectors are not generally rational. They are typically algebraic or irrational. Therefore, an exact-rational representation would still need to be converted to floating point for the spectral layer.

Using `Rational` would also make repeated matrix multiplication expensive because denominators can grow quickly.

So the core carrier is `Double`.

### Tolerance constraint

The tolerance must satisfy:

```text
floating-point normalisation error < probabilityTolerance << meaningful modelling error
```

For a row of length `n`, ordinary floating-point normalisation error is roughly on the order of:

```text
n * u
```

where `u` is double-precision machine epsilon, approximately:

```text
1.1e-16
```

For dimensions up to around `n = 1000`, this is around:

```text
1e-13
```

So:

```haskell
probabilityTolerance = 1e-9
```

is safely above expected floating-point noise while remaining far below meaningful modelling error such as a row summing to `1.01`.

### Consequence

A value of type:

```haskell
ProbabilityVector n
```

or:

```haskell
StochasticMatrix n
```

means numerical validity up to `probabilityTolerance`, not exact symbolic validity.

---

## N2 — No derived Eq for probability carriers

### Decision

Do not derive `Eq` for `ProbabilityVector n` or `StochasticMatrix n`.

Good:

```haskell
newtype StochasticMatrix n = StochasticMatrix
  { unStochasticMatrix :: Matrix Double
  }
  deriving (Show)
```

Bad:

```haskell
newtype StochasticMatrix n = StochasticMatrix
  { unStochasticMatrix :: Matrix Double
  }
  deriving (Eq, Show)
```

### Reason

Structural equality on `Double` probability data is misleading.

For example, two values may be numerically equivalent for the library’s purposes but not exactly equal as `Double`s.

### Consequence

Tests and users should use an explicit approximate comparison function:

```haskell
approxStochasticMatrixEq
  :: Double
  -> StochasticMatrix n
  -> StochasticMatrix n
  -> Bool
```

This makes the tolerance visible at the call site.

---

## M3 — Informative validation errors

### Decision

Smart constructors return:

```haskell
Either ValidationError a
```

rather than:

```haskell
Maybe a
```

The validation error type records the reason for failure.

Examples:

```haskell
NegativeEntry
EntryAboveOne
RowSumOffBy
NonSquareMatrix
MatrixDimensionMismatch
VectorDimensionMismatch
```

### Reason

`Maybe` only says whether validation succeeded or failed.

`Either ValidationError` says why validation failed.

This matters for testing and debugging. A rejection test should assert the exact validation error, not merely that validation failed.

### Consequence

A test like this is too weak:

```haskell
mkStochasticMatrix @2 matrix `shouldSatisfy` isLeft
```

A stronger test checks the exact constructor:

```haskell
mkStochasticMatrix @2 matrix
  `shouldFailWith` NegativeEntry
    { row = 1
    , col = 0
    , val = -0.1
    }
```

This ensures the matrix is rejected for the correct reason.

---

## M4 — Matrix validation delegates to vector validation

### Decision

`mkStochasticMatrix` validates each row using `mkProbabilityVectorAt`.

Conceptually:

```haskell
mkStochasticMatrix matrix =
  check matrix is n-by-n
  validate each row with mkProbabilityVectorAt
  wrap as StochasticMatrix
```

### Reason

A stochastic matrix is an `n`-tuple of probability vectors over the same state space.

Duplicating row validation inside `mkStochasticMatrix` would create two sources of truth.

### Consequence

If the probability-vector invariant changes, matrix validation automatically follows.

For example, if the tolerance policy changes in `ProbabilityVector`, stochastic matrix validation uses the same policy.

---

## M5 — Simulation step

### Decision

The simulation primitive has the shape:

```haskell
step
  :: (KnownNat n, PrimMonad m)
  => StochasticMatrix n
  -> Finite n
  -> MWC.Gen (PrimState m)
  -> m (Finite n)
```

### Reason

A DTMC step samples the next state from the transition row indexed by the current state.

The current state has type:

```haskell
Finite n
```

so it is guaranteed to be a valid state index in:

```text
0, 1, ..., n - 1
```

The result also has type:

```haskell
Finite n
```

so `step` cannot return an invalid state index.

The function is polymorphic over `PrimMonad`, rather than fixed to `IO`, so it can run in:

```haskell
IO
```

for real simulations and in:

```haskell
ST s
```

for deterministic tests with seeded random generators.

### Consequence

The simulation API does not force `IO`.

This makes simulation properties easier to test reproducibly.

---

## N3 — Sampling from stochastic rows

### Decision

`step` uses categorical sampling from `mwc-random`.

The transition row is passed as a vector of weights.

We do not re-normalise the row inside `step`.

### Reason

Rows have already been validated up to `probabilityTolerance`.

Categorical sampling works with relative weights. Therefore a row summing to:

```text
1 ± probabilityTolerance
```

is harmless for sampling.

Re-normalising inside `step` would duplicate responsibility and introduce additional floating-point operations.

### Consequence

The constructor is responsible for validation.

The simulation function is responsible for sampling.

---

## T1 — Rejection properties

### Decision

Constructor rejection tests should assert the exact validation error.

Examples:

```haskell
NegativeEntry
EntryAboveOne
RowSumOffBy
NonSquareMatrix
MatrixDimensionMismatch
```

### Reason

A test that only checks rejection is too weak.

For example, this is not enough:

```haskell
mkStochasticMatrix @2 badMatrix `shouldSatisfy` isLeft
```

A broken implementation could reject everything with the wrong error and still pass.

### Required property

If a generated stochastic matrix is perturbed off the simplex by a change larger than `probabilityTolerance`, then `mkStochasticMatrix` must return `Left` with the matching constructor.

For example:

```text
sign flip below -probabilityTolerance
  -> NegativeEntry

row-sum bump above probabilityTolerance
  -> RowSumOffBy
```

### Consequence

The tests verify not only that invalid matrices are rejected, but that they are rejected for the mathematically correct reason.

---

## T2 — Constructor round-trip property

### Decision

Generated valid stochastic matrices should round-trip through the constructor.

The property is:

```haskell
case mkStochasticMatrix @n rawMatrix of
  Right stochasticMatrix ->
    unStochasticMatrix stochasticMatrix == rawMatrix
  Left _ ->
    False
```

### Reason

The constructor should validate input. It should not mutate, repair, or normalise input.

For example, an invalid row like:

```haskell
[0.5, 0.4]
```

should be rejected, not silently normalised to:

```haskell
[0.555..., 0.444...]
```

### Consequence

The constructor boundary is honest:

```text
valid input is accepted unchanged
invalid input is rejected
```

This property also connects the generator to the tolerance policy. If generator normalisation error exceeds `probabilityTolerance`, generated matrices may be rejected. Therefore the generator and tolerance policy must be compatible.

---

## G1 — Generator scope

### Decision

The current QuickCheck stochastic matrix generator is suitable for basic constructor and multiplication properties.

It is appropriate for:

```text
constructor round-trip
multiplication closure
simulation smoke tests
```

It is not sufficient for future classification algorithms.

### Reason

Dense Dirichlet-style generators produce strictly positive rows.

A strictly positive stochastic matrix is typically regular, irreducible, and aperiodic. Such generators do not cover reducible chains, absorbing chains, closed communicating classes, or prescribed zero patterns.

The current generator intentionally allows zero entries, so it is not a true Dirichlet generator. However, the warning still applies to future dense generators.

### Consequence

Classification algorithms must use separate structured generators.

Examples of future structured generators:

```haskell
genAbsorbingStochasticMatrix
genReducibleStochasticMatrix
genStochasticMatrixWithZeroPattern
```

Do not use a dense generator as evidence that classification algorithms work on reducible or absorbing chains.

---

## S1 — Scope discipline

### Decision

Do not add generic container-style instances for stochastic carriers.

For example, do not add:

```haskell
Functor
Foldable
Traversable
```

for `StochasticMatrix n`.

### Reason

A stochastic matrix is not merely a container of numbers. It is a validated mathematical object.

Mapping an arbitrary function over entries can break the invariant.

For example:

```haskell
fmap (+ 1) matrix
```

would make row sums invalid.

Similarly:

```haskell
fmap negate matrix
```

would create negative probabilities.

### Consequence

The API should stay narrow and mathematically meaningful.

Good operations are ones that preserve the invariant by construction or revalidate at the boundary.

Examples:

```haskell
mkStochasticMatrix
mulStochasticMatrix
step
approxStochasticMatrixEq
```

Suspicious operations are ones that expose arbitrary entrywise mutation without validation.

---

## M1 status

M1 is considered complete when the following are true:

```text
- ProbabilityVector n exists.
- StochasticMatrix n exists.
- Constructors are hidden.
- Smart constructors validate dimensions and stochasticity.
- mkStochasticMatrix delegates row validation to mkProbabilityVectorAt.
- Validation errors are informative.
- Double + probabilityTolerance policy is documented.
- Eq is not derived for probability carriers.
- approxStochasticMatrixEq exists.
- Finite n is used for state indices.
- step is implemented.
- Valid stochastic matrix generators exist.
- Rejection properties assert exact error constructors.
- Round-trip property confirms constructors do not mutate input.
```

Current status:

```text
M1: complete, except any future replacement of the current generator with a true Dirichlet-backed generator should update G1.
```