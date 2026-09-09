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
import Dtmc.Analysis.Expectation (Expectation (..))
import Dtmc.Analysis.HittingTime qualified as Hitting
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
        , "classes" .= classes
        , "irreducible" .= Classification.irreducible matrix
        , "stationary" .= stationary
        , "hitting_probability" .= hittingProbabilities
        , "hitting_time" .= map expectationValue hittingTimes
        , "fundamental" .= fundamental
        , "absorption_time" .= absorptionTimes
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
    states = finiteStates :: [Finite n]
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
