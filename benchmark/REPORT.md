# Benchmark report: `dtmc` versus PyDTMC 9.0.0

## Status

The benchmark implementation builds successfully with GHC 9.10.3. On
2026-09-08, the untimed cross-language gate passed all 360 comparisons for
every family, seed, and size through 100. The largest observed absolute
discrepancy was `3.525e-12`; the largest relative discrepancy was `2.709e-14`.
Thirty-six larger datasets were deliberately skipped by the complete-output
verification limit.

These checks are not performance results. Full-machine timing results must
only be added after running:

```console
python3 benchmark/run.py all --mode full
```

No performance winner is asserted before that run completes. Generated
measurements are accompanied by `results/environment.json`, raw Criterion and
pyperf output, correctness discrepancies, robust summaries, paired speed
ratios, and log-log scaling plots.

## Scope

The comparison covers public construction, distribution evolution, matrix
powers, communicating classes, irreducibility, stationary distributions,
eventual and expected hitting, conventional absorbing-chain calculations,
occupation matrices, and simulation throughput.

Locally finite kernels, multi-target first-passage distributions, and general
per-recurrent-state absorption probabilities are excluded because PyDTMC does
not expose semantically equivalent operations. Low-outdegree fixtures still use
dense matrix storage in both implementations.

## Interpretation checklist

1. Confirm `results/verification.json` reports `passed: true`.
2. Check the recorded BLAS/LAPACK libraries and single-thread settings.
3. Use medians and dispersion, not a single run or arithmetic mean alone.
4. Compare only like-named cold, warm, or lifecycle rows.
5. Treat timeouts and memory limits as censored observations.
6. Use profiles to explain large gaps before proposing library changes.

Likely profiling areas, based on source inspection rather than measurements,
are the checked LU-plus-SVD solve contract, support-graph materialization,
simulation row conversions, and the GTH stationary solver.
