# dtmc

A small Haskell library for type-safe discrete-time Markov chains.

The project is currently experimental. The first goal is to encode basic Markov-chain invariants clearly and test them carefully.

## Core idea

A transition matrix for a finite discrete-time Markov chain is a row-stochastic matrix.

That means:

1. the matrix is square;
2. every entry is non-negative;
3. every row sums to one.

In mathematics, if `A` and `B` are row-stochastic matrices, then `AB` is also row-stochastic.

This package reflects that idea in two layers.

## Two layers of safety

### 1. Runtime-verified, then type-carried

A raw `Matrix Double` is not trusted directly.

Instead, it must pass through:

    mkStochastic :: Matrix Double -> Maybe Stochastic

If the matrix satisfies the stochastic invariant, it becomes a value of type:

    Stochastic

After that, the type carries the fact that the matrix has already been verified.

### 2. Compiler-preserved operations

Some operations preserve the stochastic invariant by theorem.

For example:

    mulStochastic :: Stochastic -> Stochastic -> Stochastic

The product of two stochastic matrices is stochastic, so the result can safely be given type `Stochastic`.

The proof is documented in the Haddock comment above the implementation.

## Build

    cabal build all

## Test

    cabal test all --test-show-details=direct

## Numerical backend

This package uses `hmatrix`, which links against system BLAS/LAPACK.

See `NUMERICS.md` for details.
