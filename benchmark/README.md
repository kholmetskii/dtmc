# `dtmc` versus PyDTMC benchmark

This directory contains a public-API benchmark of `dtmc` against PyDTMC
9.0.0. It deliberately measures equivalent mathematical results rather than
similarly named functions.

## Quick start

Python 3.12, GHC 9.10.3, Cabal, and a BLAS/LAPACK installation are required.
Using another compiler supported by `dtmc` requires regenerating the Cabal
freeze file and records a different environment.

Haskell benchmark dependencies are pinned by `cabal.project.freeze`; direct
and transitive Python package versions are pinned by
`python/requirements.lock` and captured with the results.

```console
python3 benchmark/run.py all --mode smoke
```

The smoke workflow creates an isolated Python environment, builds the Haskell
benchmark, generates five 10-state fixtures, verifies deterministic outputs,
runs short Criterion and pyperf measurements, aggregates them, and creates
plots. The full reproducible experiment is:

```console
python3 benchmark/run.py all --mode full
```

The full suite is intentionally expensive. Individual phases are available as
`bootstrap`, `generate`, `verify`, `benchmark`, `profile`, and `analyse`
commands. `profile` records Haskell RTS allocation/GC statistics and Python
`cProfile` data for stationary, hitting, committor, mean-recurrence, bounded
visit, and long-simulation representatives.

## Fairness contract

- Both implementations read the same little-endian IEEE-754 payload. SHA-256
  hashes live in the generated `data/generated/manifest.json`; input decoding
  is not timed. The committed experiment definition is `data/spec.json`.
- Inputs are valid stochastic matrices away from each library's validation
  tolerance. Invalid-input behavior is not compared.
- Finite matrices are dense in both libraries. The `low-outdegree` family has
  exact zeros but is not described as sparse storage.
- Construction measures each public constructor as defined. PyDTMC eagerly
  builds its NetworkX graph, while `dtmc` leaves its support graph lazy.
- `structure/classes-lifecycle` measures construction through the first SCC
  answer. `classes-cold` excludes transition-matrix construction and graph
  materialization. `classes-warm/access` measures cached property access;
  `classes-warm/consumed` traverses every returned state in both libraries.
- Cached PyDTMC properties and lazy `dtmc` graph/classification fields are
  recreated for each cold sample. Fixture setup is excluded from those
  samples.
- `hitting-*/cold-all-states` computes and consumes the complete all-state
  result. The `dtmc`-only `hitting-*/warm-lookup` rows measure one scalar
  lookup through the same partially applied function after its shared solve
  has been forced; PyDTMC has no equivalent reusable lookup object.
- Bounded return cases ask for the probability of returning to state zero by
  step 10 or 100. PyDTMC exposes the individual first-passage masses, so its
  timed adapter sums them. Bounded visit cases ask for the expected visits to
  the first manifest target, starting in state zero. PyDTMC's reward horizon
  is one less because it includes both time zero and the named final step.
- Pyperf records 30 fresh full-operation values per normal cell (three worker
  processes, ten values each). Warm-cache lookups are internally batched;
  Criterion uses its calibrated sample schedule. Smoke mode intentionally uses
  fewer values and is not performance evidence.
- Numeric results are consumed with deterministic checksums. This prevents
  Haskell laziness from turning a benchmark into wrapper construction alone.
- BLAS thread environment variables are fixed to one. The actual Haskell and
  NumPy linkage is captured in `results/environment.json`.
- Simulation compares throughput and validates behavior statistically; equal
  seeds do not imply equal paths because the RNG algorithms differ. PyDTMC's
  public simulation call creates its generator during the timed operation;
  `dtmc` receives the generator required by its public API before timing. The
  `dtmc` case uses its public matrix-specialized API; its lazy row preparation
  and cache are created inside the timed operation.

## Dataset families

- `dense`: every entry is positive; irreducible and aperiodic.
- `low-outdegree`: four or fewer exact-positive entries per row, including a
  self-loop and the next ring state; irreducible and aperiodic.
- `absorbing`: two singleton absorbing states and a transient dense block from
  which absorption is certain.
- `reducible`: a transient communicating block feeding two closed irreducible
  classes, with probabilities kept away from either library's near-zero
  threshold.
- `periodic`: a deterministic directed cycle, included so period and cyclic
  decomposition do nontrivial work.

Full runs use sizes 10, 25, 50, 100, 250, 500, and 1000 with seeds 1729, 2718,
and 31415. Cubic operations stop at 500, occupation matrices at 250, bounded
return and visit cases at 100, and dense NetworkX structural cases at 500.
These caps are symmetric.

## Measured mappings

| Benchmark | `dtmc` | PyDTMC |
|---|---|---|
| `construction/public-consumed` | `fromRows` | `MarkovChain(P)` |
| `evolution/one-step` | `evolveVector` | `redistribute(1)` |
| `evolution/k-*` | `evolveVectorN` | `redistribute(k)` |
| `power/*` | `power` | `to_nth_order` |
| `structure/classes-*` | `communicatingClasses` | `communicating_classes` |
| `structure/irreducible-*` | `irreducible` | `is_irreducible` |
| `structure/period-*` | `chainPeriod` | `period` |
| `structure/cyclic-classes-*` | `cyclicClasses` | `cyclic_classes` |
| `stationary` | `stationaryDistributions` | `pi` |
| `hitting-probability/cold-all-states` | all-state eventual hitting | `hitting_probabilities` |
| `hitting-time/cold-all-states` | all-state expected hitting | `hitting_times` |
| `race/forward-committor` | `raceProbabilityGivenInitialState` | `committor_probabilities("forward", ...)` |
| `return/mean-recurrence` | `expectationGivenInitialState` | `mean_recurrence_times` |
| `return/bounded/k-*` | `probabilityGivenInitialState (AtMost k)` | summed `first_passage_probabilities` |
| `visits/bounded-expectation/k-*` | `boundedExpectationGivenInitialState` | `expected_rewards(k - 1, ...)` |
| `fundamental-matrix` | `fundamentalMatrix` | `fundamental_matrix` |
| `absorption-time` | all-state expectation, transient part | `mean_absorption_times` |
| `absorption/probabilities` | per-absorber transient probabilities | `absorption_probabilities` |
| `occupation-matrix` | `occupationMatrix` | `mean_number_visits` |
| `simulation/*` | `simulateMatrix` | `simulate` |

`evolveVectorN` chooses between `k` matrix-vector steps and matrix powering
using the matrix dimension and requested step count; PyDTMC performs `k`
matrix-vector steps. The comparison is a public end-to-end workload, not a
claim that their internal primitives are identical.

Race/committor and mean-recurrence comparisons use only the dense and
low-outdegree aperiodic irreducible families because PyDTMC returns `None` for
non-ergodic chains. Absorption-probability rows use the manifest absorber order
and canonical transient-state order.

PyDTMC's `mean_number_visits` excludes the visit at time zero, while
`occupationMatrix` includes it. The Python case therefore adds the identity
matrix to the PyDTMC result. That semantic adapter is included in the public
workload timing and in correctness verification.

## Correctness

`benchmark/python/verify.py` compares complete outputs through size 100 by
default. It checks structural equality, exact infinity masks, and numeric
values using `atol=1e-10` and `rtol=1e-8`. It records maximum absolute and
relative discrepancies in `results/verification.json` and exits unsuccessfully
on a mismatch.

## Interpreting output

Raw engine output is retained under `results/raw/`. `summary.csv` contains
median, MAD, quartiles, and sample counts. `ratios.csv` reports
`PyDTMC time / dtmc time`, so values above one favor `dtmc`. Do not compare
cached and cold rows or infer sparse-storage performance from low-outdegree
fixtures. The generated SVGs include both runtime scaling and speed-ratio
plots; ratio confidence bands summarize the three deterministic seeds.
Analysis reads only raw JSON files corresponding to datasets in the current
generated manifest. Raw files from another mode are retained but ignored, so a
smoke run can safely follow a full run without deleting the full results.
