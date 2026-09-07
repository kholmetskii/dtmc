# Changelog

## 0.2.0.0

First Hackage release.

- Added validated dense and sparse probability distributions.
- Added type-safe finite transition matrices and locally finite transition
  kernels, with representation-independent finite-horizon evolution and
  simulation.
- Added finite-time joint and conditional probabilities.
- Added exact-time, bounded, eventual, competing, and expected hitting and
  return quantities.
- Added finite- and infinite-horizon visit-count analysis, including the
  occupation matrix.
- Added communicating-class, recurrence, transience, periodicity, and cyclic
  class analysis.
- Added canonical decomposition, fundamental matrices, and absorption
  probabilities and expectations.
- Added extremal stationary distributions for every recurrent class, ordinary
  limiting matrices, and cyclic subsequential limits.
- Added state-labelled and list-based construction and inspection. No
  `hmatrix` type appears in the public API.
- Changed the internal dense storage of `DistributionVector` and
  `TransitionMatrix` from statically sized values to ordinary `hmatrix`
  vectors and matrices. The public types remain state-indexed and abstract,
  and their smart constructors continue to validate dimensions against the
  finite state cardinality.
- Added `Dtmc.Distribution.Map.mapStates` for transforming sparse
  distributions, combining the weights of states that share a target.
- Added `Dtmc.Transition.Matrix.fromRows`, which builds a matrix from a grid
  of weights and reports shape mismatches as typed errors.
- Made GTH stationary-distribution normalisation robust when finite weights
  have a sum that overflows `Double`.
- Reduced dense transition-row lookup from quadratic to linear time and
  space.
- Made `Dtmc.Distribution.Vector.fromList` positional: it now takes one weight
  per state in canonical state order, so it is the exact inverse of `toList`,
  and reports a length mismatch through the new `DistributionVectorError`.
  Labelled construction, where duplicates combine and missing states default
  to zero, remains `Dtmc.Distribution.Map.fromList`.
- Supports GHC 9.6 through 9.14.
- Narrowed `Dtmc.Analysis.Classification` to the queries themselves. The
  `Classification` report and `classify` are no longer exported, and with them
  the `Of` suffixes that existed only to keep record fields from colliding
  with the standalone functions. `absorbingStates`, `chainPeriod` and
  `ergodic` are now functions on a matrix, and `communicatingClasses` returns
  `[CommClass state]`, carrying each class's period and closedness rather than
  its members alone.
