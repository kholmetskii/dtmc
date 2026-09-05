{- |
Module      : Dtmc.Analysis.FiniteTime
Description : Transition, event, and conditional probabilities.

Finite-time probability queries shared by dense finite matrices and locally
finite kernels. Kernels implement 'Transition'; initial laws may be either a
dense finite @DistributionVector@ or a @DistributionMap@ through the
'Distribution' abstraction. All calculations use finite reachable support
and perform no truncation, clamping, or renormalisation.
-}
module Dtmc.Analysis.FiniteTime (
    stepProbability,
    nStepProbability,
    probability,
    probabilityGiven,
    Observation (..),
    ConditionalProbabilityError (..),
) where

import Dtmc.Analysis.FiniteTime.Internal (
    NormalisedObservations (..),
    normalise,
 )
import Dtmc.Distribution (
    Distribution (..),
 )
import Dtmc.Distribution.Map (
    pointMass,
 )
import Dtmc.Dynamics (
    evolveN,
 )
import Dtmc.Transition (
    Transition (..),
 )
import Numeric.Natural (
    Natural,
 )

{- | A timed state observation. @At t i@ is the event @X_t = i@. A list of
observations denotes their conjunction; list order has no meaning.
-}
data Observation state
    = At Natural state -- ^ Require the supplied state at the specified time.
    deriving (Eq, Show)

-- | Why a conditional probability query has no defined value.
data ConditionalProbabilityError
    = -- | The condition has probability exactly zero.
      ZeroProbabilityCondition
    deriving (Eq, Show)

{- | Return the one-step transition probability
@P(X_1 = j | X_0 = i)@ through any locally finite 'Transition'.

Complexity: excluding 'transitionLaw', @O(log(s + 1))@ time and @O(1)@
temporary and result space for returned law support size @s@.
-}
stepProbability ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    kernel ->
    TransitionState kernel ->
    TransitionState kernel ->
    Double
stepProbability kernel source =
    probabilityAt (transitionLaw kernel source)

{- | Return the @k@-step transition probability @P(X_k = j | X_0 = i)@. At
@k = 0@ this is the Kronecker delta.

For the complexity bounds, @w@, @e@, and @u@ are per-step upper bounds on
stored source states, traversed transition edges, and accumulated destination
states, and @r@ bounds the final stored support.

Complexity: excluding 'transitionLaw',
@O(k (w + e log(u + 1) + u) + log(r + 1) + 1)@ time, @O(w + u)@ temporary
space, and @O(1)@ result space.
-}
nStepProbability ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    Natural ->
    kernel ->
    TransitionState kernel ->
    TransitionState kernel ->
    Double
nStepProbability steps kernel source =
    probabilityAt (evolveN steps (pointMass source) kernel)

{- | Compute the probability of a conjunction of timed observations.
Observation order has no meaning, duplicates collapse, and an empty
conjunction is exactly one. Contradictory observations at one time give
exactly zero without inspecting the initial distribution or kernel.

A state probability is represented by a singleton observation. A consecutive
path is represented by observations at times zero, one, and so on.

For the complexity bounds, @m@ is the supplied observation count, @q@ the
normalised count, @k@ the greatest observed time, and @s_0@ the initial stored
support size. Across all propagation steps and observation gaps, @w@, @e@,
and @u@ bound stored source states, traversed transition edges, and
accumulated destination states, while @r@ bounds support at each lookup. Let
@C = w + e log(u + 1) + u@.

Complexity: excluding the initial 'distributionWeights' call and all
'transitionLaw' evaluations,
@O(m log(m + 1) + s_0 + k C + q log(r + 1) + 1)@ time,
@O(m + w + u)@ temporary space, and @O(1)@ result space.
-}
probability ::
    ( Distribution distribution
    , Transition kernel
    , DistributionState distribution ~ TransitionState kernel
    , Ord (TransitionState kernel)
    ) =>
    distribution ->
    kernel ->
    [Observation (TransitionState kernel)] ->
    Double
probability initial kernel observations =
    case normalise [(time, state) | At time state <- observations] of
        Impossible -> 0
        Consistent [] -> 1
        Consistent ((firstTime, firstState) : rest) ->
            probabilityAt (evolveN firstTime initial kernel) firstState
                * gaps (firstTime, firstState) rest
  where
    gaps _ [] = 1
    gaps (previousTime, previousState) ((time, state) : more) =
        nStepProbability
            (time - previousTime)
            kernel
            previousState
            state
            * gaps (time, state) more

{- | Compute conditional probability @P(E | C)@ for two conjunctions of timed
observations. An exactly zero-probability condition returns
'ZeroProbabilityCondition' without evaluating the numerator; otherwise the
result is the ordinary 'Double' quotient of the joint and condition
probabilities.

For the complexity bounds, use the parameters from 'probability' across both
probability evaluations: @m@ is the total number of pairs processed, @q@ the
total normalised count, and @k@ the total number of propagation steps.
Let @C = w + e log(u + 1) + u@.

Complexity: excluding up to two initial 'distributionWeights' calls and all
'transitionLaw' evaluations,
@O(m log(m + 1) + s_0 + k C + q log(r + 1) + 1)@ time,
@O(m + w + u)@ temporary space, and @O(1)@ result space.
-}
probabilityGiven ::
    ( Distribution distribution
    , Transition kernel
    , DistributionState distribution ~ TransitionState kernel
    , Ord (TransitionState kernel)
    ) =>
    distribution ->
    kernel ->
    [Observation (TransitionState kernel)] ->
    [Observation (TransitionState kernel)] ->
    Either ConditionalProbabilityError Double
probabilityGiven initial kernel event condition =
    if denominator == 0
        then Left ZeroProbabilityCondition
        else Right (numerator / denominator)
  where
    denominator = probability initial kernel condition
    numerator = probability initial kernel (event <> condition)
