# Numerics

This package uses `hmatrix`, which links against a system BLAS/LAPACK backend.

That means numerical results are reproducible up to the behaviour of the BLAS/LAPACK implementation available on the platform. Tiny floating-point differences may occur across machines or operating systems.

## Current numerical policy

For now, `dtmc` uses a fixed tolerance when checking whether rows sum to one:

    epsilon = 1e-9

A matrix is accepted as row-stochastic when:

1. it is square;
2. all entries are non-negative;
3. every row sums to `1` up to the tolerance.

This is a pragmatic numerical check, not an exact symbolic proof.

## Future work

Possible future improvements include:

- configurable tolerances;
- clearer error messages instead of `Maybe`;
- exact rational stochastic matrices;
- separation between exact and floating-point backends.
