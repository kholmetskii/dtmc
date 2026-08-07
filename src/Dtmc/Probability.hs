{- |
Module      : Dtmc.Probability
Description : Scalar, trajectory, event, and conditional probabilities of a DTMC.

Probability queries for a finite DTMC. The scalar queries read a single number
from the chain: 'transitionProbability' and 'transitionProbabilityN' give the
one- and @k@-step transition probabilities @P(i, j)@ and @P^k(i, j)@, and
'probabilityAtTime' gives the time-@k@ marginal @P(X_k = j)@ from an initial
law. The event queries build on them: 'pathProbability' scores an explicit
consecutive trajectory, while 'probability' and 'conditionalProbability'
score conjunctions of timed observations @(X_t = i)@ at arbitrary times, using
the Markov property and time homogeneity to reduce each query to the scalar
functions.

All results are ordinary 'Double' arithmetic: no clamping to @[0, 1]@,
renormalisation, or revalidation is performed.
-}
module Dtmc.Probability (
    transitionProbability,
    transitionProbabilityN,
    probabilityAtTime,
    Observation (..),
    FiniteObservation,
    ProbabilityError (..),
    pathProbability,
    probability,
    conditionalProbability,
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
import Dtmc.Dynamics (
    evolveN,
 )
import Dtmc.Probability.Internal (
    NormalisedObservations (..),
    normalise,
 )
import Dtmc.TransitionMatrix (
    TransitionMatrix,
    matrixPower,
    rowAt,
 )
import GHC.TypeNats (
    KnownNat,
 )
import Numeric.Natural (
    Natural,
 )

{- | A single timed state observation. @At t i@ is the event @X_t = i@: the
chain occupies state @i@ at time @t@. The state type is deliberately
polymorphic so the same event vocabulary can be reused by future countable
chain representations; the finite-chain queries in this module specialise it
to 'FiniteObservation'.

A list of observations is read as their /conjunction/, e.g. @[At 2 c, At 5 d]@
is the event @X_2 = c@ and @X_5 = d@. The list order carries no meaning;
'probability' and 'conditionalProbability' sort by time internally.
-}
data Observation state
    = At Natural state
    deriving (Eq, Show)

-- | A timed observation whose state is bounded by a finite chain's dimension.
type FiniteObservation n = Observation (Finite n)

{- | Why a probability query has no defined value. The only cause is
conditioning on an event of probability exactly zero, for which @P(E | C)@ is
undefined.
-}
data ProbabilityError
    = -- | The condition has probability exactly zero.
      ZeroProbabilityCondition
    deriving (Eq, Show)

{- | The one-step transition probability @P(i, j) = P(X_1 = j | X_0 = i)@ from
source @i@ to destination @j@. By construction it reads the destination
coordinate of the source row:

@
transitionProbability p i j == 'probabilityAt' ('rowAt' p i) j
@

Both 'Finite' indices are bounded, so the query is total. The stored entry is
returned exactly, with no clamping, renormalisation, or revalidation.

Time: @O(n)@, dominated by extracting the source row.
-}
transitionProbability ::
    (KnownNat n) =>
    TransitionMatrix n ->
    -- | source @i@
    Finite n ->
    -- | destination @j@
    Finite n ->
    Double
transitionProbability p i = probabilityAt (rowAt p i)

{- | The @k@-step transition probability
@P^k(i, j) = P(X_k = j | X_0 = i)@, read from the @k@-th matrix power:

@
transitionProbabilityN k p i j == transitionProbability ('matrixPower' k p) i j
@

At @k = 0@ the power is the identity, so the result is the Kronecker delta:
one when @i == j@ and zero otherwise. Chapman-Kolmogorov holds in exact
arithmetic; the computed 'Double' is neither revalidated nor clamped, so
rounding from the repeated-squaring power is preserved.

Time: @O(n^2 + n^3 log(k + 1))@, dominated by 'matrixPower'.
-}
transitionProbabilityN ::
    (KnownNat n) =>
    Natural ->
    TransitionMatrix n ->
    -- | source @i@
    Finite n ->
    -- | destination @j@
    Finite n ->
    Double
transitionProbabilityN k p = transitionProbability (matrixPower k p)

{- | The marginal probability @P(X_k = j)@ that the chain occupies state @j@
after @k@ steps, started from the initial law @initial@ and driven by @p@. It
reads the target coordinate of the evolved distribution:

@
probabilityAtTime k initial p j == 'probabilityAt' ('evolveN' k initial p) j
@

At @k = 0@ no step has been taken, so the result is
@'probabilityAt' initial j@, the initial probability of @j@.

This is distinct from 'probabilityAt': 'probabilityAt' reads a coordinate from
an already-computed distribution, whereas 'probabilityAtTime' first evolves
@initial@ through @k@ steps of @p@ and then reads coordinate @j@. The 'Double'
result is not clamped, renormalised, or revalidated.

Time: @O(n^2 + n^3 log(k + 1))@, dominated by the matrix power inside
'evolveN'.
-}
probabilityAtTime ::
    (KnownNat n) =>
    Natural ->
    Distribution n ->
    TransitionMatrix n ->
    Finite n ->
    Double
probabilityAtTime k initial p = probabilityAt (evolveN k initial p)

{- | The probability that the chain follows the exact trajectory
@(i_0, ..., i_m)@ over consecutive times @0, 1, ..., m@:

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

{- | The probability of the conjunction of timed observations. The list is read
as an /and/ of the events @X_t = i@ and may be given in any order.

Writing the normalised, time-sorted observations as
@[At t_0 i_0, At t_1 i_1, ..., At t_k i_k]@ with @t_0 < t_1 < ... < t_k@, the
Markov property and time homogeneity give

@
probability lambda p obs
    == 'probabilityAtTime' t_0 lambda p i_0
         * product [ 'transitionProbabilityN' (t_r - t_(r-1)) p i_(r-1) i_r | r <- [1 .. k] ]
@

that is @P(X_{t_0} = i_0) * prod_{r=1}^{k} (P^{t_r - t_{r-1}})_{i_{r-1}, i_r}@.

Normalisation before scoring:

* observations are sorted by ascending time;
* exact duplicates such as @[At 2 c, At 2 c]@ collapse to @[At 2 c]@ and do not
  change the result;
* observations demanding different states at one time, such as
  @[At 2 c, At 2 d]@ with @c /= d@, describe an impossible (but valid) event
  and return exactly @0@.

Boundary cases:

* the empty conjunction returns exactly @1@ (the sure event);
* a single observation @[At t i]@ equals @'probabilityAtTime' t lambda p i@;
* observations at consecutive times @0, 1, ..., k@ agree with
  'pathProbability'.

The result is ordinary 'Double' arithmetic, with no clamping to @[0, 1]@,
renormalisation, or revalidation.

Time: @O(m log m)@ to sort @m@ observations, then one initial distribution
evolution (via 'probabilityAtTime') and one transition-matrix power per
distinct time gap (via 'transitionProbabilityN'). Each such power costs
@O(n^3 log t)@ for its exponent @t@, and none are shared across gaps.
-}
probability ::
    (KnownNat n) =>
    Distribution n ->
    TransitionMatrix n ->
    [FiniteObservation n] ->
    Double
probability initial p observations =
    case normalise [(t, i) | At t i <- observations] of
        Impossible -> 0
        Consistent [] -> 1
        Consistent ((t0, i0) : rest) ->
            probabilityAtTime t0 initial p i0 * gaps (t0, i0) rest
  where
    gaps _ [] = 1
    gaps (tPrev, iPrev) ((t, i) : more) =
        transitionProbabilityN (t - tPrev) p iPrev i * gaps (t, i) more

{- | The conditional probability @P(E | C)@ of an event @E@ given a condition
@C@, both conjunctions of timed observations. It is defined by

@
P(E | C) = P(E and C) / P(C)
@

and computed by composing 'probability':

@
denominator = probability initial p condition
numerator   = probability initial p (event <> condition)
@

Observations shared by @event@ and @condition@ are collapsed while normalising
the combined list for the numerator, so conditioning an observation on itself
yields @1@.

When the denominator is /exactly/ zero the ratio is undefined and the result is
@'Left' 'ZeroProbabilityCondition'@. The zero test is exact: no epsilon,
clamping, or snapping is applied, and tolerated input and floating-point error
are treated exactly as elsewhere in the library. Otherwise the result is
@'Right' (numerator / denominator)@ using ordinary 'Double' division.

Boundary cases:

* @conditionalProbability initial p event []@ is
  @'Right' (probability initial p event)@, since the empty condition is
  the sure event and has probability @1@;
* with a condition of positive probability, @conditionalProbability initial p
  [] condition@ is @'Right' 1@;
* if the event conflicts with a possible condition (different states at one
  time), the numerator is @0@ and the result is @'Right' 0@;
* if the condition is internally contradictory, or otherwise has probability
  exactly zero, the result is @'Left' 'ZeroProbabilityCondition'@.

Time: two 'probability' calls, so @O(m log m)@ sorting plus the matrix
powers described there, for @m@ the combined observation count.
-}
conditionalProbability ::
    (KnownNat n) =>
    Distribution n ->
    TransitionMatrix n ->
    -- | event @E@
    [FiniteObservation n] ->
    -- | condition @C@
    [FiniteObservation n] ->
    Either ProbabilityError Double
conditionalProbability initial p event condition =
    if denominator == 0
        then Left ZeroProbabilityCondition
        else Right (numerator / denominator)
  where
    denominator = probability initial p condition
    numerator = probability initial p (event <> condition)
