{- |
Module      : Dtmc.Analysis.FiniteTime
Description : Scalar, trajectory, event, and conditional probabilities.

Finite-time probability queries shared by dense finite matrices and locally
finite kernels. Kernels implement 'Transition'; initial laws may be either a
dense finite @DistributionVector@ or a @DistributionMap@ through the
'Distribution' abstraction. All calculations use finite reachable support
and perform no truncation, clamping, or renormalisation.
-}
module Dtmc.Analysis.FiniteTime (
    transitionProbability,
    transitionProbabilityN,
    probabilityAtTime,
    Observation (..),
    ConditionalProbabilityError (..),
    pathProbability,
    jointProbability,
    conditionalProbability,
) where

import Data.List.NonEmpty (
    NonEmpty ((:|)),
 )
import Dtmc.Analysis.FiniteTime.Internal (
    NormalisedObservations (..),
    normalise,
 )
import Dtmc.Distribution (
    Distribution (..),
 )
import Dtmc.Distribution.Map (
    pointMass,
    toDistributionMap,
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
    = At Natural state
    deriving (Eq, Show)

-- | Why a conditional probability query has no defined value.
data ConditionalProbabilityError
    = -- | The condition has probability exactly zero.
      ZeroProbabilityCondition
    deriving (Eq, Show)

{- | One-step transition probability @P(i,j)@ through any 'Transition'. An
absent destination has probability exactly zero.
-}
transitionProbability ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    kernel ->
    TransitionState kernel ->
    TransitionState kernel ->
    Double
transitionProbability kernel source =
    probabilityAt (transitionLaw kernel source)

{- | The @k@-step transition probability @P^k(i,j)@, computed by evolving a
point mass for exactly @k@ sparse steps. At @k = 0@ this is the Kronecker
delta. No state-space enumeration or truncation is performed.
-}
transitionProbabilityN ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    Natural ->
    kernel ->
    TransitionState kernel ->
    TransitionState kernel ->
    Double
transitionProbabilityN steps kernel source =
    probabilityAt (evolveN steps (pointMass source) kernel)

{- | Marginal probability @P(X_k = j)@ after @k@ transitions. The initial law
is converted to an equivalent map-backed finite-support representation once
for this query.
-}
probabilityAtTime ::
    ( Distribution distribution
    , Transition kernel
    , DistributionState distribution ~ TransitionState kernel
    , Ord (TransitionState kernel)
    ) =>
    Natural ->
    distribution ->
    kernel ->
    TransitionState kernel ->
    Double
probabilityAtTime steps initial kernel =
    probabilityAt (evolveN steps initial kernel)

{- | Probability of an explicit consecutive non-empty trajectory. A one-state
path returns its initial probability; longer paths multiply the initial mass
by every one-step transition probability.
-}
pathProbability ::
    ( Distribution distribution
    , Transition kernel
    , DistributionState distribution ~ TransitionState kernel
    , Ord (TransitionState kernel)
    ) =>
    distribution ->
    kernel ->
    NonEmpty (TransitionState kernel) ->
    Double
pathProbability initial kernel (initialState :| rest) =
    probabilityAt (toDistributionMap initial) initialState
        * go initialState rest
  where
    go _ [] = 1
    go previous (next : more) =
        transitionProbability kernel previous next * go next more

{- | Joint probability of a conjunction of timed state observations.
Observations are sorted by time, duplicates collapse, conflicting states at
one time return exactly zero, and the empty conjunction returns exactly one.
-}
jointProbability ::
    ( Distribution distribution
    , Transition kernel
    , DistributionState distribution ~ TransitionState kernel
    , Ord (TransitionState kernel)
    ) =>
    distribution ->
    kernel ->
    [Observation (TransitionState kernel)] ->
    Double
jointProbability initial kernel observations =
    case normalise [(time, state) | At time state <- observations] of
        Impossible -> 0
        Consistent [] -> 1
        Consistent ((firstTime, firstState) : rest) ->
            probabilityAtTime firstTime initial kernel firstState
                * gaps (firstTime, firstState) rest
  where
    gaps _ [] = 1
    gaps (previousTime, previousState) ((time, state) : more) =
        transitionProbabilityN
            (time - previousTime)
            kernel
            previousState
            state
            * gaps (time, state) more

{- | Conditional probability @P(E | C)@ for conjunctions of timed
observations. An exactly zero condition returns
@Left ZeroProbabilityCondition@; otherwise ordinary 'Double' division is
used without an epsilon or clamping.
-}
conditionalProbability ::
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
conditionalProbability initial kernel event condition =
    if denominator == 0
        then Left ZeroProbabilityCondition
        else Right (numerator / denominator)
  where
    denominator = jointProbability initial kernel condition
    numerator = jointProbability initial kernel (event <> condition)
