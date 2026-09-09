{-# LANGUAGE ScopedTypeVariables #-}

module Force (
    checksumClasses,
    checksumExpectations,
    checksumFundamental,
    checksumMatrix,
    checksumOccupation,
    checksumStationary,
    checksumStates,
    checksumVector,
    eitherOrFail,
) where

import Data.Finite (Finite, getFinite)
import Dtmc.Analysis.Classification (CommClass (..))
import Dtmc.Analysis.Expectation (Expectation (..))
import Dtmc.Distribution.Vector qualified as Vector
import Dtmc.State (FiniteState, stateIndex)
import Dtmc.Transition.Matrix qualified as Matrix

checksumValues :: [Double] -> Double
checksumValues = foldl' (+) 0

checksumVector :: Vector.DistributionVector state -> Double
checksumVector = checksumValues . Vector.toList

checksumMatrix :: Matrix.TransitionMatrix state -> Double
checksumMatrix = checksumValues . concat . Matrix.toRows

checksumClasses :: (FiniteState state) => [CommClass state] -> Double
checksumClasses classes =
    checksumValues
        [ fromIntegral (getFinite (stateIndex member))
        | communicatingClass <- classes
        , member <- classMembers communicatingClass
        ]

checksumStates :: [Finite n] -> Double
checksumStates =
    checksumValues . map (fromIntegral . getFinite)

checksumExpectations :: [Expectation] -> Double
checksumExpectations = foldl' step 0
  where
    step total expectation =
        let value = case expectation of
                FiniteExpectation finite -> finite
                InfiniteExpectation -> 1
            next = total + value
         in next

checksumFundamental :: ([state], [[Double]]) -> Double
checksumFundamental (_, rows) = checksumValues (concat rows)

checksumOccupation :: [[Expectation]] -> Double
checksumOccupation = checksumExpectations . concat

checksumStationary ::
    (FiniteState state) =>
    [([state], Vector.DistributionVector state)] ->
    Double
checksumStationary = foldl' step 0
  where
    step total (members, distribution) =
        let memberTotal =
                checksumValues
                    (map (fromIntegral . getFinite . stateIndex) members)
            next = total + memberTotal + checksumVector distribution
         in next

eitherOrFail :: (Show problem) => Either problem value -> value
eitherOrFail = either (error . ("benchmark operation failed: " ++) . show) id
