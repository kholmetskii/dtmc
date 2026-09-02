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
    = At Natural state
    deriving (Eq, Show)

-- | Why a conditional probability query has no defined value.
data ConditionalProbabilityError
    = -- | The condition has probability exactly zero.
      ZeroProbabilityCondition
    deriving (Eq, Show)

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

{- | Probability of a conjunction of timed observations. Observation order
has no meaning, duplicates collapse, and an empty conjunction is exactly one.

A state probability is represented by a singleton observation. A consecutive
path is represented by observations at times zero, one, and so on.
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

{- | Conditional probability @P(E | C)@ for two conjunctions of timed
observations. An exactly zero-probability condition returns
'ZeroProbabilityCondition'; otherwise the result is the ordinary 'Double'
quotient of joint probabilities.
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
