# Numerical honesty

Every trade of accuracy for tractability gets an entry: where, what was chosen,
what was rejected, and how it is checked.

Architectural decisions live in [docs/DECISIONS.md](docs/DECISIONS.md) (`D*`).
Test design lives in [docs/TESTING.md](docs/TESTING.md) (`T*`).

---

## N1 — `mulStochasticMatrix` does not re-validate its result

**Where.** `Dtmc.StochasticMatrix.mulStochasticMatrix`.

**Choice.** The result is returned through the raw constructor, without
`mkStochasticMatrix`.

**Why.** The closure theorem is exact. Float row sums of the product equal 1 to
within `O(n·u)`, `u ≈ 1.1e-16`. Running a provably correct result through the
strict constructor trades a proof for a tolerance gamble that can only lose
(spurious rejection).

**Consequence.** The signature is total: no `Either`. The proof pays rent. If a
function carries a `-- Proof:` block and still returns `Maybe`/`Either` for the
property it supposedly proves, the proof carries nothing.

**Checked by.** `Dtmc.StochasticMatrixSpec.prop_productRowStochastic`, which
*does* re-validate — appropriate in a test, not in the implementation.

---

## N2 — `simplexTolerance = 1e-9`, absolute

**Where.** `Dtmc.Simplex`.

**Choice.** `Double` plus an ε.

**Why not `Rational`.** Eigenvalues of a stochastic matrix are algebraic and
generically irrational; there is no exact rational eigendecomposition. By M6 we
would have to coerce to `Double` anyway, having paid for `Rational` with
denominator blow-up under repeated multiplication in `matPow`. The spectral
layer is intrinsically floating-point.

**Why 1e-9.** It is bounded from both sides:

```
generator normalisation error  <  ε  <<  meaningful invalidity
        ~ n·u                              (what a model would notice)
```

The lower bound is forced by the round-trip property: the generator normalises
rows in `Double`, so their sums land in `1 ± O(n·u)`. Any tighter and "generated
matrices are always accepted" fails. At `n ≤ 10³` that is `≈ 1e-13` — four orders
of margin. Absolute equals relative here because the target is exactly 1.

**What this means.** `StochasticMatrix n` does NOT mean "the rows sum to 1". It
means "no row deviates from 1 by more than ε, and this was checked once at a
constructor boundary".

**Open.** Revisit past `n > 10³`.

---

## N3 — naive summation

**Where.** `Dtmc.Simplex.validateSimplexPoint`, `total = sum entries`.

**Choice.** Naive summation, error `O(n·u)`. Rejected: Kahan/Neumaier, `O(u)`.

**Why.** At `n ≤ 10³` the naive error is `≈ 1e-13`, four orders below `ε`. Kahan
would buy accuracy that ε immediately eats. Revisit together with N2.

---

## N4 — `categorical` normalises the weights itself

**Where.** `Dtmc.Simulation.sampleFrom`.

`System.Random.MWC.Distributions.categorical` normalises the weight vector
internally, so a sum of `1 ± ε` is harmless and the row must NOT be renormalised
before sampling. Recorded so that nobody later "fixes" this by dividing by the sum.

`clampToleratedNegatives` zeroes coordinates in `[-ε, 0)` (rounding noise) and
throws on anything below: that is a violated invariant, i.e. a bug, not data.

---

## N5 — the BLAS backend is not pinned

**Where.** The whole library: hmatrix links against the system LAPACK/BLAS.

Results are reproducible only up to the BLAS implementation. Linux gives
OpenBLAS/reference, macOS gives Accelerate.

**Checked by.** CI: Linux on every push, macOS/Accelerate on a weekly schedule.
The old CI ran on `macos-14` only — that is, on exactly one, atypical backend, so
the claim was asserted but never tested.

**Open.** `cabal freeze` plus `index-state`, so CI cannot go red from an upstream
release.

---

## N6 — two tolerance regimes coexist; never conflate them

**Where.** `simplexTolerance` versus the Monte Carlo tolerance arriving in M2.

- Floating point: `ε = 1e-9`.
- Monte Carlo: `O(1/√N)`; at `N = 2·10⁵` that is `≈ 1.1e-3` — **six orders of
  magnitude looser**.

That is why `approxDistributionEq` and `approxStochasticMatrixEq` take the
tolerance as an **explicit parameter** rather than reading `simplexTolerance`.
Reach for `simplexTolerance` in the empirical-marginal convergence test and CI is
red forever; reach for the MC tolerance in a constructor and the constructor
accepts garbage.

**Naming rule for M2.** `epsFloat` and `mcTol`, never a bare `tol`.

---

## Deferred to M6

Naive matrix powering, `O(n³ log k)`, is numerically BENIGN: for a row-stochastic
`P` we have `‖P‖_∞ = 1` exactly, the norm is submultiplicative, so `‖P^k‖_∞ = 1`
for every `k`. Powering is non-expansive; errors accumulate additively as
`O(k·n·u)` rather than geometrically as they do for a general matrix (`‖A‖^k`).

The spectral route buys speed (`O(n³)` once, then `O(n²)` per power) and pays in
the conditioning `κ(V)` of the eigenvector matrix, which blows up for non-normal
`P` and near-defective spectra — chains close to periodic, where eigenvalues on
the unit circle nearly collide.

Static provides a general non-symmetric solver:
`instance KnownNat n => Eigen (Sq n) (C n) (M n n)`, i.e.
`eigensystem :: Sq n -> (C n, M n n)`. No detour into the dynamic API is needed.

**Note for M2.** `Numeric.LinearAlgebra.Static` has NO `<#`. It offers only
`(<>)`, `(#>)`, `(<.>)`. So `λP` must be written `tr p #> λ`. This is precisely
why `pushforward` must be the only exported operation: a bare `#>` would compute
`Pλ`, and π would come out a right eigenvector instead of a left one.
