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
    stepProbability,
    nStepProbability,
    stateProbability,
    pathProbability,
    observationProbability,
    conditionalObservationProbability,
    Observation (..),
    ConditionalProbabilityError (..),
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
        stepProbability kernel previous next * go next more

{- | One-step transition probability @P(X_1 = j | X_0 = i)@ through any
locally finite 'Transition'.
-}
stepProbability ::
    (Transition kernel, Ord (TransitionState kernel)) =>
    kernel ->
    TransitionState kernel ->
    TransitionState kernel ->
    Double
stepProbability kernel source =
    probabilityAt (transitionLaw kernel source)

{- | The @k@-step transition probability @P(X_k = j | X_0 = i)@. At @k = 0@
this is the Kronecker delta.
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

{- | State probability @P(X_k = j)@ under an initial distribution.
-}
stateProbability ::
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
stateProbability steps initial kernel =
    probabilityAt (evolveN steps initial kernel)

{- | Probability of a conjunction of timed observations. Observation order
has no meaning, duplicates collapse, and an empty conjunction is exactly one.
-}
observationProbability ::
    ( Distribution distribution
    , Transition kernel
    , DistributionState distribution ~ TransitionState kernel
    , Ord (TransitionState kernel)
    ) =>
    distribution ->
    kernel ->
    [Observation (TransitionState kernel)] ->
    Double
observationProbability initial kernel observations =
    case normalise [(time, state) | At time state <- observations] of
        Impossible -> 0
        Consistent [] -> 1
        Consistent ((firstTime, firstState) : rest) ->
            stateProbability firstTime initial kernel firstState
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

{- | Conditional observation probability @P(E | C)@. An exactly
zero-probability condition returns 'ZeroProbabilityCondition'; otherwise the
result is the ordinary 'Double' quotient of joint probabilities.
-}
conditionalObservationProbability ::
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
conditionalObservationProbability initial kernel event condition =
    if denominator == 0
        then Left ZeroProbabilityCondition
        else Right (numerator / denominator)
  where
    denominator = observationProbability initial kernel condition
    numerator = observationProbability initial kernel (event <> condition)
