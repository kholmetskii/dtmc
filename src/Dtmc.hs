{- |
Module      : Dtmc
Description : Orientation and module map for discrete-time Markov chains.

This package models time-homogeneous discrete-time Markov chains (DTMCs).
It supports two complementary representations:

* finite chains over a 'Dtmc.State.FiniteState' use
  'Dtmc.Transition.Matrix.TransitionMatrix' and may use dense
  'Dtmc.Distribution.Vector.DistributionVector' values;
* locally finite chains over unrestricted state types use
  'Dtmc.Transition.Kernel.TransitionKernel' and sparse
  'Dtmc.Distribution.Map.DistributionMap' values.

Start with these modules:

* "Dtmc.State" for finite named state types;
* "Dtmc.Distribution.Map" and "Dtmc.Distribution.Vector" for validated
  probability laws;
* "Dtmc.Transition.Kernel" and "Dtmc.Transition.Matrix" for transition
  models;
* "Dtmc.Dynamics" and "Dtmc.Simulation" for evolution and sampling.

Analysis is organised by mathematical subject:

* "Dtmc.Analysis.FiniteTime" for transition, joint, and conditional
  probabilities;
* "Dtmc.Analysis.HittingTime", "Dtmc.Analysis.ReturnTime", and
  "Dtmc.Analysis.VisitCount" for path-time and occupation quantities;
* "Dtmc.Analysis.Classification" for communication, recurrence, and
  periodicity;
* "Dtmc.Analysis.Absorption", "Dtmc.Analysis.Stationary", and
  "Dtmc.Analysis.Limiting" for finite-chain long-run behaviour.

The package deliberately has no broad facade of re-exports because several
analysis modules use the same concise names, such as @probability@ and
@expectation@, for their subject-specific operations. Import analysis modules
qualified.

The ordinary construction and inspection APIs do not expose @hmatrix@ types.
Callers already using @hmatrix@ can opt into
"Dtmc.Distribution.Vector.HMatrix" and
"Dtmc.Transition.Matrix.HMatrix". The package still uses @hmatrix@ internally
and therefore requires a BLAS/LAPACK implementation when it is built.
-}
module Dtmc () where
