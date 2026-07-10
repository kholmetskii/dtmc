# Test and generator design

Entries `T*`. Falsified by a test that should have failed and did not.

Numerical decisions live in [../NUMERICS.md](../NUMERICS.md) (`N*`).
Architectural decisions live in [DECISIONS.md](DECISIONS.md) (`D*`).

---

## T1 — bump the SMALLEST coordinate, not the first

**Where.** `Dtmc.Generators.bumpSmallest`.

**History.** QuickCheck (seed 180161876) failed two properties on the
counterexample `[1.0, 0.0, 0.0]`. Bumping the first coordinate gives
`1.000001 > 1 + ε`, and the coordinatewise check — which by contract runs BEFORE
the sum — returns `EntryAboveOne`, not the expected `SumOffBy`.

**Fix.** For a simplex point `min ≤ 1/n ≤ 0.5`, so after `+1e-6` the smallest
coordinate stays safely below `1 + ε`. The error is isolated to `SumOffBy`.

**Side conclusion.** Two failures from one seed, in `SimplexSpec` and
`DistributionSpec`, are direct evidence that the predicate is shared (D3).

---

## T2 — a 3-cycle, not `I₂`

**Where.** `Dtmc.KernelSpec`, `Dtmc.SimulationSpec`.

**Problem.** The old test used the identity matrix `I₂`. It is **symmetric**:
row `i` equals column `i`. Had `rowAt` used `toColumns` instead of `toRows`, both
tests would still have passed.

At `n = 2` even the permutation matrix `[[0,1],[1,0]]` is symmetric. You need
`n = 3`:

```
P = [[0,1,0],[0,0,1],[1,0,0]]     orbit under P : 0 → 1 → 2 → 0
                                  orbit under Pᵀ: 0 → 2 → 1 → 0
```

**Why this matters more than it looks.** The λ ↦ λP convention is why π is a LEFT
eigenvector. The right eigenvector also exists, is also pretty, and is the vector
of ones. Confusing them breaks all of M5/M6 without producing a single red build.

---

## T3 — a generator with zeros, not Dirichlet

**Where.** `Dtmc.Generators.genSimplexPointList`.

**Trap.** Dirichlet(1,…,1) yields strictly positive rows. A strictly positive
stochastic matrix is regular: irreducible AND aperiodic, hence ergodic. So
**every** generated matrix would be ergodic, and the M3 classification tests
would be vacuously green — they would never see a reducible or absorbing chain.

**Fix.** 30% of coordinates are exactly zero
(`frequency [(3, pure 0), (7, choose (0,1000))]`).

**Remaining limitation.** The zeros appear BY CHANCE, not by construction.
Covering reducible chains in M3 still requires structured generators with a
prescribed zero pattern. Do not rely on this generator there.

---

## T4 — `counterexample`, not `_ -> False`

**Where.** every `*Spec.hs`.

**Problem.** If a property asserts the **exact error constructor**, the
"wrong constructor" branch must print what it got. Otherwise QuickCheck shows the
input (`[1.0, 0.0, 0.0]`) but stays silent about the result
(`EntryAboveOne 0 1.000001`), and the diagnosis takes ten minutes instead of one.

**Rule.**

```haskell
case validateSimplexPoint v of
  Left (SumOffBy _) -> property True
  other -> counterexample ("expected SumOffBy, got " <> show other) False
```

`isLeft` is not a test. A property that checks only the fact of rejection stays
green under any broken check.

---

## T5 — a fixed seed, not a random one

**Where.** `Dtmc.SimulationSpec` today; M2's `empiricalMarginal` shortly.

**Rule.** A statistical test with a random seed has a nonzero failure rate by
construction. Flaky tests destroy trust in the whole prove/verify claim faster
than any real bug does.

So simulation tests are `it`, not `prop`, and use `MWC.create` (fixed seed) inside
`runST` rather than `createSystemRandom`. This is also the only place where
`step`'s `PrimMonad` polymorphism is exercised: production path `IO`, test path `ST`.

**For M2.** The simulation test is the only NON-circular check of
Chapman–Kolmogorov: `P^(m+n) ≈ P^m · P^n` merely tests the associativity of
floating-point multiplication, because both sides are built from the same
`matPow` operation. The temptation to reach for a random-seeded `prop` will be
strong. Resist it.

Tolerance: a `5σ` band, `|p̂_j − (λPⁿ)_j| < 5·√(p̂_j(1−p̂_j)/N)`. Changing the seed
to make it pass is falsification, not debugging.

**Do not take `n` too large.** As `n` grows, `λPⁿ → π`, and the test starts
verifying the stationary distribution rather than Chapman–Kolmogorov: it stops
discriminating between correct and slightly wrong dynamics.
