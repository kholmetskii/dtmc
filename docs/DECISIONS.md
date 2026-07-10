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
- **Choice.** `Distribution n = Distribution (R n)`, `TransitionMatrix n = TransitionMatrix (Sq n)`.
- **Rejected.** `newtype TransitionMatrix (n :: Nat) = TransitionMatrix (Matrix Double)` —
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
- **Consequence.** Migrating to Static does NOT by itself close
  `coerce :: TransitionMatrix 2 -> TransitionMatrix 3`. `Coercible` between the
  same type constructor at a phantom parameter is solved by the lifting rule,
  without unwrapping the newtype, so the hidden constructor does not help.
- **Decision.** `type role Distribution nominal`, `type role TransitionMatrix nominal`.
  A role annotation can only strengthen the inferred role.
- **Residual assumption.** At the hmatrix level, `coerce :: Sq 2 -> Sq 3` remains
  possible. Our guarantees hold under the assumption that users do not coerce
  `Sq`/`R` directly. Documented rather than hidden.
- **Checked by.** `check/Role.hs` must FAIL to compile: `ghc -isrc check/Role.hs`.

---

## D3 — one predicate, owned by `Distribution`; a one-way dependency

- **Where.** `Dtmc.Distribution.validateSimplex`, `Dtmc.TransitionMatrix.mkTransitionMatrix`.
- **Choice.** The Δ^{n-1} predicate and `simplexTolerance` live in
  `Dtmc.Distribution`. `mkTransitionMatrix` *calls* it row by row and wraps the
  error as `InRow i`.
- **Rejected.** (a) A separate `Dtmc.Simplex` module — an extra file for a
  predicate that *is* the definition of a distribution. (b) Independent
  validation in each carrier, sharing only an ε from a `Config` module.
- **Why (b) fails.** What gets duplicated is not the number but three decisions:
  `x < -ε` rather than `x < 0`; coordinatewise BEFORE the sum; left to right.
  The copies drift, and rejection tests start exercising different copies. The
  definition (ST227 §2.4) — a matrix is row-stochastic iff every row is a
  distribution — is one predicate, not a similar property.
- **Evidence.** The counterexample `[1.0, 0.0, 0.0]` failed two properties on the
  same seed (180161876), because the predicate is one.
- **Why not a `Config` module.** `ε` is a conclusion drawn from a proof (N2), not
  a setting. A separate file turns the conclusion into a knob that gets tuned to
  make a test green — and it attracts the MC tolerance, six orders away (N6).
- **Direction.** `TransitionMatrix → Distribution` is legitimate (`has-a`).
  `Distribution → TransitionMatrix` is forbidden: Δ(S) knows nothing of kernels.
- **Checked by.** CI, job `discipline`, step "Distribution does not depend on
  TransitionMatrix".
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
- **Why a module at all.** Haskell has no visibility modifier on a data
  constructor. Visibility is the module export list, i.e. a file. So "the
  constructor is visible to a couple of modules but not to the user" has exactly
  one encoding. This is `internal` in C#, package-private in Java, `pub(crate)`
  in Rust. Precedent: `Data.Text.Internal`; `Internal.Static` inside hmatrix.
- **Why not fold the carriers into `Distribution.hs`.** The set of theorems that
  may construct without checking grows; the carrier definition does not. Today:
  `mulTransitionMatrix`, `rowAt`. M2 adds `pushforward`, `matPow`; M4
  `absorptionProbabilities`; M5 `stationary`; M7 the MH kernel. Without
  `Internal` they all must live in the file that defines the newtype, and within
  six months `Distribution.hs` is a two-thousand-line god module where spectral
  theory sits next to simplex validation.
- **What it buys beyond `private`.** An audit list:
  `grep -rl '^import Dtmc.Internal' src/` is exactly the set of places where the
  invariant is held by a proof rather than by a check.
- **Boundary.** `Dtmc.Distribution` owns Δ^{n-1} and constructs freely. Every
  other importer must carry a `-- Proof:`.
- **Rejected.** A subdirectory `src/Dtmc/Distribution/Internal.hs`. It existed so
  that the glob `src/Dtmc/*.hs` would not see it — i.e. the file hierarchy was
  doing the `.cabal` file's job. It was also asymmetric and promised two
  directories of one file each.
- **Checked by.** CI, step "Dtmc.Internal is imported only at the boundary or
  under a Proof".
- **Limitation.** The check looks for a `-- Proof:` anywhere in the file, not
  directly above the constructor application. A strict check would require
  parsing Haskell. The job catches the real risk: a NEW module quietly acquiring
  access.

---

## D5 — `rowAt` is a method on the matrix, not a module

- **Where.** `Dtmc.TransitionMatrix.rowAt`.
- **Claim.** `P : S → Δ(S)`, `i ↦ P_{i·}`. So `TransitionMatrix n` is morally
  isomorphic to `Finite n -> Distribution n`. Without `rowAt` the library has no
  chain — only a validated square array.
- **Rejected.** A separate `Dtmc.Kernel` module holding `rowAt`. It existed to
  keep `Distribution` and `TransitionMatrix` mutually independent — but mutual
  independence was never the goal, and it is mathematically wrong: a transition
  matrix *is* a family of distributions. The dependency is one-way (D3), so
  `rowAt` belongs on the matrix, exactly as `matrix.getRow(i)` would.
- **What it buys.** `sampleFrom :: Distribution n -> ...` is true by type —
  "sample from a distribution", not "from an arbitrary vector of numbers".
  `step p i = sampleFrom (rowAt p i)` reads as a definition and is total.
- **Convention.** Distributions are ROWS acting on the left (λ ↦ λP). Hence π is
  a LEFT eigenvector, `πP = π`.
- **Checked by.** `Dtmc.TransitionMatrixSpec` (3-cycle, T2).

---

## D6 — `rowAt` builds every row to get one

- **Where.** `Dtmc.TransitionMatrix.rowAt`, `S.toRows (...) !! i`.
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
- **Consequence.** `Dtmc.Internal` and `Dtmc.Errors` need no subdirectory in
  order to "hide from the glob".
- **Open.** Revisit at M3, which adds `Chain`, `Property`, `Witness`. Five or six
  types earn a `Types/` directory; two do not. `structures/`, `data/`, `objects/`
  were rejected: these are not containers, not data blobs, and not objects.
- **Checked by.** CI, step "Every exposed module has a paired spec".

---

## D8 — the guarantee is stated in two tiers, explicitly

- **Where.** README, and the Haddock of `Dtmc.TransitionMatrix`.
- **Statement.**
  - *Proved by the compiler:* dimension, squareness, dimension agreement.
  - *Checked once at a boundary, then carried by the type:* stochasticity.
- **Why it matters.** "The compiler proves stochasticity" is false:
  `S.matrix [1,2,-3,0.5] :: Sq 2` is perfectly valid. This is the one objection
  that can be raised against the project in earnest, and it is answered by an
  honest statement rather than by silence.

---

## D9 — no `StochasticRow` type

- **Considered.** A hidden `StochasticRow n` carrier, with `Distribution` and
  the rows of `TransitionMatrix` both built on it.
- **Rejected.** A row of a transition matrix *is* a distribution. Both denote the
  same set Δ^{n-1}. A second name for the same set is a synonym, not an
  abstraction.
- **The decisive argument.** `TransitionMatrix n` is a dense `Sq n` — a pointer to
  a contiguous buffer that LAPACK multiplies in one BLAS call. It does not
  *contain* n rows; it contains n² doubles. `S.toRows` **constructs** a fresh
  `R n`. So a row is not what the matrix is made of, it is what gets carved out
  of it — and the result of carving is already a `Distribution`.
- **The alternative that would justify it.**
  `newtype TransitionMatrix n = TransitionMatrix (Vec n (StochasticRow n))` makes
  "every row is a distribution" true by construction. But then
  `mulTransitionMatrix` rebuilds a dense matrix on every call, and M2 (`matPow`)
  and M6 (`eigensystem :: Sq n -> (C n, M n n)`) die. Representation is dense for
  BLAS; the invariant is a boundary check because in floating point `Σ = 1`
  cannot be structural anyway (it is not even transitive under addition).
- **Falsification test.** What breaks if `StochasticRow` is removed?
  `mkDistribution` works. `mkTransitionMatrix` works, delegating.
  `rowAt` works. No ε is duplicated. Nothing breaks — so the type carries nothing.

---

## D10 — `Dtmc.Errors` as shared vocabulary, kept internal

- **Where.** `Dtmc.Errors` (in `other-modules`).
- **Choice.** `SimplexError` is the common core; `DistributionError` and
  `TransitionError` wrap it with their own context. Users obtain the constructors
  by re-export from the module that raises them.
- **Why internal.** Keeping it out of `exposed-modules` means it needs no paired
  spec (it has no logic), and the CI rule stays free of an exemption list.
- **Why the wrappers differ.** A standalone distribution has no rows, so
  `DistributionError` has no row index. `TransitionError` does: `InRow 2 (SumOffBy 1.03)`
  reads as "row 2 violates the distribution invariant". The old
  `RowSumOffBy { row = 0 }` for a standalone vector was simply a lie.