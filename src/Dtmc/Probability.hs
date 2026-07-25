{- |
Module      : Dtmc.Probability
Description : Probabilities of concrete trajectories of a DTMC.

Probabilities of explicit finite trajectories of a finite DTMC. For an initial
law @lambda@ and one-step matrix @P@, the probability of following the path
@(i_0, ..., i_m)@ is @lambda(i_0) * prod_{r=1}^{m} P(i_{r-1}, i_r)@.

This module is intentionally small: later milestones extend it with timed
joint and conditional probability queries.
-}
module Dtmc.Probability (
    pathProbability,
) where

import Data.Finite (
    Finite,
 )
import Data.List.NonEmpty (
    NonEmpty ((:|)),
 )
import Dtmc.Distribution (
    Distribution,
    probabilityAt,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    transitionProbability,
 )
import GHC.TypeNats (
    KnownNat,
 )

{- | The probability that the chain follows the exact trajectory
@(i_0, ..., i_m)@:

@
pathProbability lambda p (i_0 :| [i_1, ..., i_m])
    == 'probabilityAt' lambda i_0 * product [ 'transitionProbability' p i_(r-1) i_r | r <- [1 .. m] ]
@

which is @lambda(i_0) * prod_{r=1}^{m} P(i_{r-1}, i_r)@, the joint law
@P(X_0 = i_0, X_1 = i_1, ..., X_m = i_m)@ under the Markov property.

Edge behaviour:

* A one-state path @i_0 :| []@ returns @'probabilityAt' lambda i_0@, the
  initial probability of that state (the empty product is one).
* An impossible step, where some @P(i_{r-1}, i_r)@ is zero, makes the whole
  product zero without any special casing.
* 'NonEmpty' rules out the empty path at the type level, so there is no
  partial or undefined case.

The result is ordinary 'Double' arithmetic: factors are multiplied as stored,
with no clamping to @[0, 1]@, renormalisation, or revalidation.

Time: @O(m * n)@ for a path of @m@ steps over @n@ states, dominated by the
@m@ row reads.
-}
pathProbability ::
    (KnownNat n) =>
    Distribution n ->
    TransitionMatrix n ->
    NonEmpty (Finite n) ->
    Double
pathProbability initial p (i0 :| rest) =
    probabilityAt initial i0 * go i0 rest
  where
    go _ [] = 1
    go prev (next : more) =
        transitionProbability p prev next * go next more
