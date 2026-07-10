# Architectural decisions

Format: where, choice, rejected, why, evidence, checked by, open.

Prefix `D*`, so decision numbers cannot collide with milestone numbers `M*`.
Numerical decisions live in [../NUMERICS.md](../NUMERICS.md) (`N*`).
Test design lives in [TESTING.md](TESTING.md) (`T*`).

Membership rule: **what would falsify this entry?**
A number → `N`. A counterexample in the types → `D`. A failing test → `T`.

---

## D1 — carriers from `Numeric.LinearAlgebra.Static`, not a phantom `Nat`

- **Where.** `Dtmc.Internal`.
- **Choice.** `Distribution n = Distribution (R n)`, `StochasticMatrix n = StochasticMatrix (Sq n)`.
- **Rejected.** `newtype StochasticMatrix (n :: Nat) = StochasticMatrix (Matrix Double)` —
  `n` never occurs on the right, i.e. a sticker on a box.
- **Why.** Three of the six constructors of the old `ValidationError` were
  dimension errors (`VectorDimensionMismatch`, `NonSquareMatrix`,
  `MatrixDimensionMismatch`). With `Sq n` they are unrepresentable: their
  existence was a confession that the types were not working. `LA.<>` on a
  dynamic matrix is a runtime dimension check; `S.<>` is a compile-time one.
- **Evidence.** Three error constructors and all the `Proxy`/`natVal` plumbing
  deleted. Two tests (`rejects a non-square matrix`, `rejects a wrong
  type-level dimension`) no longer compile.
- **Cost.** The dynamic→static boundary is
  `S.create :: Matrix Double -> Maybe (Sq n)`, returning `Maybe`. A size mismatch
  and a stochasticity failure are different kinds of failure; mixing them in one
  type was the original mistake.
- **Checked by.** The compiler.

---

## D2 — `type role ... nominal` is mandatory; Static does not provide it

- **Where.** `Dtmc.Internal`.
- **Observation.** In hmatrix:
  ```haskell
  newtype Dim (n :: Nat) t = Dim t          -- n does not occur on the right
  newtype R n   = R (Dim n (Vector ℝ))
  newtype L m n = L (Dim m (Dim n (Matrix ℝ)))
  ```
  and `grep -rn "type role" hmatrix/` finds **nothing**. Static's dimension
  indices have phantom roles.
- **Consequence 1.** Migrating to Static does NOT by itself close
  `coerce :: StochasticMatrix 2 -> StochasticMatrix 3`. `Coercible` between the
  same type constructor at a phantom parameter is solved by the lifting rule,
  without unwrapping the newtype, so the hidden constructor does not help.
- **Decision.** `type role Distribution nominal`, `type role StochasticMatrix nominal`.
  A role annotation can only strengthen the inferred role.
- **Residual assumption.** At the hmatrix level, `coerce :: Sq 2 -> Sq 3` remains
  possible. Our guarantees hold under the assumption that users do not coerce
  `Sq`/`R` directly. Documented rather than hidden.
- **Checked by.** `check/Role.hs` must FAIL to compile: `ghc -isrc check/Role.hs`.

---

## D3 — one shared predicate Δ^{n-1}, independent error types

- **Where.** `Dtmc.Simplex`, `Dtmc.Distribution`, `Dtmc.StochasticMatrix`.
- **Choice.** `validateSimplexPoint` and `simplexTolerance` live in one module.
  The error types differ: `DistributionError` carries no row index,
  `StochasticError` does (`InRow`).
- **Rejected.** Independent validation sharing an `ε` from a `Config` module.
- **Why.** What gets duplicated is not the number but three decisions: `x < -ε`
  rather than `x < 0`; coordinatewise BEFORE the sum; left to right. The copies
  drift, and rejection tests start exercising different copies. The definition
  (ST227 §2.4) — a matrix is row-stochastic iff every row is a distribution — is
  one predicate, not a similar property.
- **Evidence.** The counterexample `[1.0, 0.0, 0.0]` failed `SimplexSpec` and
  `DistributionSpec` on the same seed (180161876), because the predicate is one.
- **Why not a `Config` module.** `ε` is a conclusion drawn from a proof (N2), not
  a setting. A separate file turns the conclusion into a knob that gets tuned to
  make a test green — and it attracts the MC tolerance, six orders away (N6).
- **Checked by.** CI, job `discipline`, step "No import edge".
- **Open.** `SimplexError` is deliberately transparent: rejection tests match on
  the exact constructor. Close the constructors if an external consumer appears.

---

## D4 — `Dtmc.Internal`: the raw constructor only under a `-- Proof:`

- **Where.** `Dtmc.Internal` (in `other-modules`).
- **Problem.** A hidden constructor is the validation boundary. But theorems
  sometimes establish the invariant without a check, and they need the
  constructor. `rowAt` returning `Either DistributionError (Distribution n)`
  would be a lie: failure is impossible. `either (error "impossible") id . mkDistribution`
  is worse: an `O(n)` re-check on every simulation step, plus an `error` in a
  place where we have a proof.
- **Choice.** Both carriers in one internal module. `Distribution` and
  `StochasticMatrix` are the validation boundary and do construct. Every other
  importer must carry a `-- Proof:`.
- **Rejected.** A subdirectory `src/Dtmc/Distribution/Internal.hs`. It existed so
  that the glob `src/Dtmc/*.hs` would not see it — i.e. the file hierarchy was
  doing the `.cabal` file's job. It was also asymmetric (`StochasticMatrix` had
  no such directory) and promised two directories of one file each.
- **Precedent.** `Data.Text.Internal`; `Internal.Static` inside hmatrix itself.
- **Checked by.** CI, step "Dtmc.Internal is imported only at the boundary or
  under a Proof".
- **Limitation.** The check looks for a `-- Proof:` anywhere in the file, not
  directly above the constructor application. A strict check would require
  parsing Haskell. The job catches the real risk: a NEW module quietly acquiring
  access.

---

## D5 — `Dtmc.Kernel`: the matrix as a Markov kernel

- **Where.** `Dtmc.Kernel.rowAt`.
- **Claim.** `P : S → Δ(S)`, `i ↦ P_{i·}`. So `StochasticMatrix n` is morally
  isomorphic to `Finite n -> Distribution n`.
- **Why a separate module.** `Distribution` and `StochasticMatrix` know nothing
  of each other (D3), so the connecting operation must live below both. A module
  holding one function is not a smell when the function *is* the concept.
- **What it buys.** `sampleFrom :: Distribution n -> ...` is true by type —
  "sample from a distribution", not "from an arbitrary vector of numbers".
  `step p i = sampleFrom (rowAt p i)` reads as a definition and is total.
- **Convention.** Distributions are ROWS acting on the left (λ ↦ λP). Hence π is
  a LEFT eigenvector, `πP = π`.
- **Checked by.** `Dtmc.KernelSpec` (3-cycle, T2).

---

## D6 — `rowAt` builds every row to get one

- **Where.** `Dtmc.Kernel.rowAt`, `S.toRows (...) !! i`.
- **Cost.** `O(n²)` per call instead of `O(n)`. Irrelevant at `n = 3`; visible in
  M2's `empiricalMarginal` (2·10⁵ paths × 10 steps).
- **Decision.** Deferred until M2 profiling. Do not optimise blind.

---

## D7 — a flat module tree; the source of truth is `.cabal`

- **Choice.** Zero subdirectories under `src/Dtmc/`.
- **Why.** The "every module has a paired spec" check must read
  `exposed-modules` from `dtmc.cabal`, not glob the filesystem. A glob is a
  source of truth that lies: it cannot distinguish `other-modules` from the
  public API.
- **Consequence.** `Dtmc.Internal` needs no subdirectory in order to "hide from
  the glob".
- **Checked by.** CI, step "Every exposed module has a paired spec".

---

## D8 — the guarantee is stated in two tiers, explicitly

- **Where.** README, and the Haddock of `Dtmc.StochasticMatrix`.
- **Statement.**
  - *Proved by the compiler:* dimension, squareness, dimension agreement.
  - *Checked once at a boundary, then carried by the type:* stochasticity.
- **Why it matters.** "The compiler proves stochasticity" is false:
  `S.matrix [1,2,-3,0.5] :: Sq 2` is perfectly valid. This is the one objection
  that can be raised against the project in earnest, and it is answered by an
  honest statement rather than by silence.
