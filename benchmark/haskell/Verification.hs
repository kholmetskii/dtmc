{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Verification (writeVerification) where

import Data.Aeson (Value (Null), encode, object, toJSON, (.=))
import Data.ByteString.Lazy qualified as BL
import Data.Finite (Finite, getFinite)
import Data.Proxy (Proxy)
import Dataset
import Dtmc.Analysis.Absorption qualified as Absorption
import Dtmc.Analysis.Classification qualified as Classification
import Dtmc.Analysis.Event (DiscreteEvent (AtMost, EqualTo, GreaterThan))
import Dtmc.Analysis.Expectation (Expectation (..))
import Dtmc.Analysis.FiniteTime qualified as FiniteTime
import Dtmc.Analysis.HittingTime qualified as Hitting
import Dtmc.Analysis.ReturnTime qualified as Return
import Dtmc.Analysis.Stationary qualified as Stationary
import Dtmc.Analysis.VisitCount qualified as VisitCount
import Dtmc.Distribution.Vector qualified as Vector
import Dtmc.Dynamics qualified as Dynamics
import Dtmc.State (finiteStates)
import Dtmc.Transition.Matrix qualified as Matrix
import Force (eitherOrFail)
import GHC.TypeNats (KnownNat)

writeVerification :: FilePath -> FilePath -> Int -> [SomeEntry] -> IO ()
writeVerification outputPath dataRoot maxSize entries = do
    results <- traverse (verifySome dataRoot maxSize) entries
    BL.writeFile
        outputPath
        ( encode
            ( object
                [ "schema_version" .= (1 :: Int)
                , "maximum_size" .= maxSize
                , "datasets" .= results
                ]
            )
        )

verifySome :: FilePath -> Int -> SomeEntry -> IO Value
verifySome dataRoot maxSize (SomeEntry (_ :: Proxy n) entry)
    | entrySize entry > maxSize =
        pure
            ( object
                [ "id" .= entryId entry
                , "size" .= entrySize entry
                , "skipped" .= True
                ]
            )
    | otherwise = verifyDataset <$> loadDataset @n dataRoot entry

verifyDataset :: forall n. (KnownNat n) => Dataset n -> Value
verifyDataset dataset =
    object
        [ "id" .= entryId entry
        , "family" .= family
        , "size" .= entrySize entry
        , "seed" .= entrySeed entry
        , "sha256" .= entrySha256 entry
        , "evolve_1" .= Vector.toList (Dynamics.evolveVector initial matrix)
        , "evolve_10" .= Vector.toList (Dynamics.evolveVectorN 10 initial matrix)
        , "power_10" .= Matrix.toRows (Matrix.power 10 matrix)
        , "finite_time_step" .= FiniteTime.stepProbability matrix initialState visitTarget
        , "finite_time_n_step_10" .= finiteTimeNStep 10
        , "finite_time_n_step_100" .= finiteTimeNStep 100
        , "finite_time_observation_10" .= finiteTimeObservation 10
        , "finite_time_observation_100" .= finiteTimeObservation 100
        , "classes" .= classes
        , "irreducible" .= Classification.irreducible matrix
        , "stationary" .= stationary
        , "hitting_probability" .= hittingProbabilities
        , "hitting_time" .= map expectationValue hittingTimes
        , "hitting_bounded_exact_10" .= hittingBounded EqualTo 10
        , "hitting_bounded_exact_100" .= hittingBounded EqualTo 100
        , "hitting_bounded_at_most_10" .= hittingBounded AtMost 10
        , "hitting_bounded_at_most_100" .= hittingBounded AtMost 100
        , "hitting_bounded_greater_than_10" .= hittingBounded GreaterThan 10
        , "hitting_bounded_greater_than_100" .= hittingBounded GreaterThan 100
        , "race_probability" .= raceProbabilities
        , "return_mean" .= map expectationValue returnMeans
        , "return_bounded_10" .= returnBounded 10
        , "return_bounded_100" .= returnBounded 100
        , "visits_bounded_expectation_10" .= visitsBounded 10
        , "visits_bounded_expectation_100" .= visitsBounded 100
        , "chain_period" .= chainPeriod
        , "cyclic_classes" .= cyclicClasses
        , "fundamental" .= fundamental
        , "absorption_time" .= absorptionTimes
        , "absorption_probability" .= absorptionProbabilities
        , "occupation" .= occupation
        ]
  where
    entry = datasetEntry dataset
    family = entryFamily entry
    matrix :: Matrix.TransitionMatrix (Finite n)
    matrix = eitherOrFail (Matrix.fromRows (datasetRows dataset))
    initial :: Vector.DistributionVector (Finite n)
    initial = eitherOrFail (Vector.fromList (datasetInitialWeights dataset))
    targets = datasetTargets dataset
    competing = datasetCompeting dataset
    states = finiteStates :: [Finite n]
    initialState =
        case states of
            [] -> error "verification requires a non-empty state space"
            first : _ -> first
    visitTarget =
        case targets of
            [] -> error "verification requires at least one target"
            target : _ -> target
    finiteTimeNStep steps =
        FiniteTime.nStepProbability steps matrix initialState visitTarget
    finiteTimeObservation time =
        FiniteTime.probability
            initial
            matrix
            [FiniteTime.At time visitTarget]
    hittingBounded event horizon =
        Hitting.probabilityGivenInitialState
            (event horizon)
            matrix
            isTarget
            initialState
    isTarget state = state `elem` targets
    index :: Finite n -> Int
    index = fromIntegral . getFinite
    classes =
        [ map index (Classification.classMembers communicatingClass)
        | communicatingClass <- Classification.communicatingClasses matrix
        ]
    stationary =
        [ object
            [ "members" .= map index members
            , "weights" .= Vector.toList distribution
            ]
        | (members, distribution) <-
            eitherOrFail (Stationary.stationaryDistributions matrix)
        ]
    hittingProbabilityAt = Hitting.eventualProbabilityGivenInitialState matrix targets
    hittingProbabilities =
        eitherOrFail (sequence [hittingProbabilityAt state | state <- states])
    hittingTimeAt = Hitting.expectationGivenInitialState matrix targets
    hittingTimes = eitherOrFail (sequence [hittingTimeAt state | state <- states])
    raceProbabilities
        | family `notElem` ["dense", "low-outdegree"] = Null
        | otherwise =
            let atState = Hitting.raceProbabilityGivenInitialState matrix targets competing
             in toJSON (eitherOrFail (sequence [atState state | state <- states]))
    returnMeanAt = Return.expectationGivenInitialState matrix
    returnMeans = eitherOrFail (sequence [returnMeanAt state | state <- states])
    returnBounded bound =
        Return.probabilityGivenInitialState (AtMost bound) matrix initialState
    visitsBounded bound =
        VisitCount.boundedExpectationGivenInitialState
            bound
            initialState
            matrix
            (== visitTarget)
    chainPeriod
        | family /= "periodic" = Null
        | otherwise = toJSON (Classification.chainPeriod matrix)
    cyclicClasses
        | family /= "periodic" = Null
        | otherwise =
            case Classification.cyclicClasses matrix of
                Nothing -> Null
                Just groups -> toJSON (map (map index) groups)
    fundamental
        | family /= "absorbing" = Null
        | otherwise =
            let (transient, values) = eitherOrFail (Absorption.fundamentalMatrix matrix)
             in object ["states" .= map index transient, "values" .= values]
    absorptionTimes
        | family /= "absorbing" = Null
        | otherwise =
            let atState = Absorption.expectationGivenInitialState matrix
                values = eitherOrFail (sequence [atState state | state <- states])
             in toJSONExpectations values
    absorptionProbabilities
        | family /= "absorbing" = Null
        | otherwise =
            let absorbing = datasetAbsorbing dataset
                transient = Classification.transientStates matrix
                valuesFor target =
                    eitherOrFail
                        ( sequence
                            [ Absorption.probabilityGivenInitialState matrix target state
                            | state <- transient
                            ]
                        )
             in object
                    [ "absorbing" .= map index absorbing
                    , "transient" .= map index transient
                    , "values" .= map valuesFor absorbing
                    ]
    occupation
        | family `notElem` ["absorbing", "reducible"] = Null
        | otherwise =
            let values = eitherOrFail (VisitCount.occupationMatrix matrix)
             in toJSONExpectationMatrix values

expectationValue :: Expectation -> Value
expectationValue expectation =
    case expectation of
        FiniteExpectation value -> object ["finite" .= value]
        InfiniteExpectation -> Null

toJSONExpectations :: [Expectation] -> Value
toJSONExpectations =
    toJSON . map expectationValue

toJSONExpectationMatrix :: [[Expectation]] -> Value
toJSONExpectationMatrix =
    toJSON . map (map expectationValue)
