# dtmc

Type-safe discrete-time Markov chains in Haskell.

## What the type system guarantees

The claim is two-tiered, and the difference matters.

**Proved by the compiler:** dimension, squareness, dimension agreement under
multiplication. That is `Numeric.LinearAlgebra.Static` — `R n`, `Sq n`, `KnownNat`.

**Checked once at a boundary, then carried by the type:** stochasticity. The
`TransitionMatrix` constructor is hidden; the only way to obtain a value is
`mkTransitionMatrix`, which runs every row through `Distribution`'s `validateSimplex`.

The compiler does **not** prove stochasticity: `S.matrix [1,2,-3,0.5] :: Sq 2` is
a perfectly valid `Sq 2`. And `TransitionMatrix n` does not mean "the rows sum to
1"; it means "no row deviates from 1 by more than ε, and this was checked once at
a single audited boundary". The justification for ε is [NUMERICS.md](NUMERICS.md), N2.

hmatrix's dimension indices have phantom roles, so migrating to `Static` does not
by itself close `coerce :: TransitionMatrix 2 -> TransitionMatrix 3`. The carriers
are annotated `type role ... nominal`; the regression test is `check/Role.hs`,
which **must fail to compile** ([docs/DECISIONS.md](docs/DECISIONS.md), D2).

## Prove / verify

Every theorem carries a `-- Proof:` block (a prose derivation) and a paired
property or simulation. The proof pays rent: it licenses a total signature.

```haskell
mulTransitionMatrix :: KnownNat n => TransitionMatrix n -> TransitionMatrix n -> TransitionMatrix n
--                                                                              ^ no Either: closure is proved
rowAt :: KnownNat n => TransitionMatrix n -> Finite n -> Distribution n
--                                                       ^ no Either: the row is already validated
```

If a function carries a `-- Proof:` and still returns `Either` for the property it
supposedly proves, the proof carries nothing.

The discipline is enforced by CI, not by memory: every exposed module has a paired
spec; `Dtmc.Internal` is imported only at the validation boundary or under a
`-- Proof:`; the `Distribution` → `TransitionMatrix` dependency direction never
reverses.

## Numerical honesty

Every trade of accuracy for tractability gets an entry: where, what was chosen,
what was rejected, how it is checked.

- [NUMERICS.md](NUMERICS.md) — numerical decisions (`N*`)
- [docs/DECISIONS.md](docs/DECISIONS.md) — architecture and types (`D*`)
- [docs/TESTING.md](docs/TESTING.md) — test and generator design (`T*`)

## Building

```bash
cabal build all
cabal test all --test-show-details=direct
ghc -isrc check/Role.hs      # must FAIL to compile
```

hmatrix links against the system LAPACK/BLAS. Linux: `libblas-dev liblapack-dev`.
macOS: Accelerate out of the box. Results are reproducible only up to the backend,
so CI exercises both ([NUMERICS.md](NUMERICS.md), N5).