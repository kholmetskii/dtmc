module Dtmc.Analysis.ProbabilityOracle (
    TruncatedLaw,
    transitionWeight,
    trajectoryProbability,
    stateProbability,
    observationProbability,
    hittingLaw,
    returnLaw,
    visitLawBefore,
    raceProbabilityWithin,
    lawProbability,
    lawUnresolvedMass,
    lawFiniteExpectation,
) where

import Data.Finite (
    getFinite,
 )
import Data.List (
    findIndex,
 )
import Data.Map.Strict (
    Map,
 )
import Data.Map.Strict qualified as Map
import Dtmc.Analysis.Event (
    DiscreteEvent (..),
    includesInfiniteOutcome,
    matchesDiscreteEvent,
 )
import Dtmc.State (
    FiniteState,
    finiteStates,
    stateIndex,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    toRows,
 )
import Numeric.Natural (
    Natural,
 )

data WeightedPath state = WeightedPath [state] Double

{- | A finite prefix of a discrete law. 'lawUnresolvedMass' is the probability
whose event time is strictly beyond the stored horizon, including any atom at
infinity. This test-only type deliberately does not appear in the library API.
-}
data TruncatedLaw = TruncatedLaw
    { lawHorizon :: Natural
    , lawFiniteMasses :: Map Natural Double
    , lawUnresolvedMass :: Double
    }

toIndex :: (FiniteState state) => state -> Int
toIndex = fromIntegral . getFinite . stateIndex

transitionWeight ::
    (FiniteState state) =>
    TransitionMatrix state ->
    state ->
    state ->
    Double
transitionWeight matrix source destination =
    toRows matrix !! toIndex source !! toIndex destination

iterateNatural :: Natural -> (value -> value) -> value -> value
iterateNatural steps advance = go steps
  where
    go 0 value = value
    go remaining value = go (remaining - 1) (advance value)

weightedTrajectories ::
    (FiniteState state) =>
    Natural ->
    [(state, Double)] ->
    TransitionMatrix state ->
    [WeightedPath state]
weightedTrajectories steps initial matrix =
    iterateNatural steps advance initialPaths
  where
    initialPaths =
        [ WeightedPath [state] weight
        | (state, weight) <- initial
        , weight /= 0
        ]
    advance paths = paths >>= extend
    extend (WeightedPath path weight) =
        [ WeightedPath (path <> [destination]) (weight * probability)
        | destination <- finiteStates
        , let probability = transitionWeight matrix (last path) destination
        , probability /= 0
        ]

trajectoryProbability ::
    (FiniteState state) =>
    [(state, Double)] ->
    TransitionMatrix state ->
    [state] ->
    Double
trajectoryProbability _ _ [] = 0
trajectoryProbability initial matrix (first : rest) =
    initialMass first * go first rest
  where
    initialMass state =
        sum [weight | (candidate, weight) <- initial, candidate == state]
    go _ [] = 1
    go previous (next : more) =
        transitionWeight matrix previous next * go next more

stateProbability ::
    (FiniteState state) =>
    Natural ->
    [(state, Double)] ->
    TransitionMatrix state ->
    state ->
    Double
stateProbability time initial matrix destination =
    sum
        [ weight
        | WeightedPath path weight <- weightedTrajectories time initial matrix
        , last path == destination
        ]

observationProbability ::
    (FiniteState state) =>
    Natural ->
    [(state, Double)] ->
    TransitionMatrix state ->
    [(Natural, state)] ->
    Double
observationProbability horizon initial matrix observations =
    sum
        [ weight
        | WeightedPath path weight <- weightedTrajectories horizon initial matrix
        , all (matches path) observations
        ]
  where
    matches path (time, expected) =
        path !! fromIntegral time == expected

lawFromFirstOccurrence ::
    Natural ->
    [WeightedPath state] ->
    ([state] -> Maybe Natural) ->
    TruncatedLaw
lawFromFirstOccurrence horizon paths occurrence =
    TruncatedLaw horizon masses unresolved
  where
    (masses, unresolved) = foldr addPath (Map.empty, 0) paths
    addPath (WeightedPath path weight) (known, unknown) =
        case occurrence path of
            Nothing -> (known, unknown + weight)
            Just time -> (Map.insertWith (+) time weight known, unknown)

hittingLaw ::
    (FiniteState state) =>
    Natural ->
    TransitionMatrix state ->
    (state -> Bool) ->
    state ->
    TruncatedLaw
hittingLaw horizon matrix isTarget initial =
    lawFromFirstOccurrence horizon paths firstHit
  where
    paths = weightedTrajectories horizon [(initial, 1)] matrix
    firstHit path = fromIntegral <$> findIndex isTarget path

returnLaw ::
    (FiniteState state) =>
    Natural ->
    TransitionMatrix state ->
    state ->
    TruncatedLaw
returnLaw horizon matrix initial =
    lawFromFirstOccurrence horizon paths firstReturn
  where
    paths = weightedTrajectories horizon [(initial, 1)] matrix
    firstReturn path =
        fromIntegral . (+ 1) <$> findIndex (== initial) (drop 1 path)

visitLawBefore ::
    (FiniteState state) =>
    Natural ->
    [(state, Double)] ->
    TransitionMatrix state ->
    (state -> Bool) ->
    TruncatedLaw
visitLawBefore bound initial matrix isVisited =
    TruncatedLaw bound masses 0
  where
    steps
        | bound == 0 = 0
        | otherwise = bound - 1
    paths = weightedTrajectories steps initial matrix
    count path =
        fromIntegral (length (filter isVisited (take (fromIntegral bound) path)))
    masses =
        Map.fromListWith
            (+)
            [(count path, weight) | WeightedPath path weight <- paths]

raceProbabilityWithin ::
    (FiniteState state) =>
    Natural ->
    TransitionMatrix state ->
    (state -> Bool) ->
    (state -> Bool) ->
    state ->
    Double
raceProbabilityWithin horizon matrix isSuccessful isCompeting initial =
    sum
        [ weight
        | WeightedPath path weight <-
            weightedTrajectories horizon [(initial, 1)] matrix
        , wins path
        ]
  where
    wins path =
        case (findIndex isSuccessful path, findIndex isCompeting path) of
            (Just successfulTime, Just competingTime) ->
                successfulTime < competingTime
            (Just _, Nothing) -> True
            _ -> False

lawProbability :: DiscreteEvent -> TruncatedLaw -> Maybe Double
lawProbability event law
    | eventKnown event (lawHorizon law) =
        Just (finiteMass + unresolvedContribution)
    | otherwise = Nothing
  where
    finiteMass =
        sum
            [ mass
            | (value, mass) <- Map.toList (lawFiniteMasses law)
            , matchesDiscreteEvent event value
            ]
    unresolvedContribution
        | includesInfiniteOutcome event = lawUnresolvedMass law
        | otherwise = 0

eventKnown :: DiscreteEvent -> Natural -> Bool
eventKnown event horizon =
    case event of
        EqualTo threshold -> threshold <= horizon
        LessThan threshold -> threshold <= horizon + 1
        AtMost threshold -> threshold <= horizon
        GreaterThan threshold -> threshold <= horizon
        AtLeast threshold -> threshold <= horizon + 1

lawFiniteExpectation :: TruncatedLaw -> Maybe Double
lawFiniteExpectation law
    | lawUnresolvedMass law == 0 =
        Just
            ( sum
                [ fromIntegral value * mass
                | (value, mass) <- Map.toList (lawFiniteMasses law)
                ]
            )
    | otherwise = Nothing
