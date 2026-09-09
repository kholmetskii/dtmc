{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Cases (benchmarksFor) where

import Control.DeepSeq (NFData (rnf))
import Control.Exception (evaluate)
import Criterion.Main
import Data.Finite (Finite)
import Data.Proxy (Proxy (Proxy))
import Data.Vector.Unboxed qualified as U
import Dataset
import Dtmc.Analysis.Absorption qualified as Absorption
import Dtmc.Analysis.Classification qualified as Classification
import Dtmc.Analysis.HittingTime qualified as Hitting
import Dtmc.Analysis.Stationary qualified as Stationary
import Dtmc.Analysis.VisitCount qualified as VisitCount
import Dtmc.Distribution.Vector qualified as Vector
import Dtmc.Dynamics qualified as Dynamics
import Dtmc.Simulation qualified as Simulation
import Dtmc.State (finiteStates)
import Dtmc.Transition.Matrix qualified as Matrix
import Force
import GHC.Exts (RealWorld)
import GHC.TypeNats (KnownNat)
import Numeric.Natural (Natural)
import System.Random.MWC qualified as MWC

newtype Prepared n = Prepared {unPrepared :: Matrix.TransitionMatrix (Finite n)}

instance NFData (Prepared n) where
    rnf (Prepared matrix) = checksumMatrix matrix `seq` ()

data PreparedDynamics n = PreparedDynamics
    { dynamicsMatrix :: Matrix.TransitionMatrix (Finite n)
    , dynamicsInitial :: Vector.DistributionVector (Finite n)
    }

instance NFData (PreparedDynamics n) where
    rnf prepared =
        checksumMatrix (dynamicsMatrix prepared)
            `seq` checksumVector (dynamicsInitial prepared)
            `seq` ()

newtype PreparedGen = PreparedGen {unPreparedGen :: MWC.Gen RealWorld}

instance NFData PreparedGen where
    rnf (PreparedGen generator) = generator `seq` ()

benchmarksFor ::
    forall n.
    (KnownNat n) =>
    Proxy n ->
    FilePath ->
    Entry ->
    Benchmark
benchmarksFor _ dataRoot entry =
    env (loadDataset @n dataRoot entry) $ \dataset ->
        bgroup
            (entryId entry)
            ( commonBenchmarks entry dataset
                ++ absorptionBenchmarks entry dataset
                ++ occupationBenchmarks entry dataset
                ++ simulationBenchmarks entry dataset
            )

commonBenchmarks :: forall n. (KnownNat n) => Entry -> Dataset n -> [Benchmark]
commonBenchmarks entry dataset =
    [ bgroup
        "construction"
        [ bench "public-consumed" $
            nf (checksumMatrix . matrixFromRows @n) (datasetRows dataset)
        ]
    , env (preparedDynamics dataset) $ \ ~(PreparedDynamics matrix initial) ->
        bgroup
            "evolution"
            [ bench "one-step" $
                nf (checksumVector . Dynamics.evolveVector initial) matrix
            , bench "k-10" $
                nf (checksumVector . Dynamics.evolveVectorN 10 initial) matrix
            , bench "k-100" $
                nf (checksumVector . Dynamics.evolveVectorN 100 initial) matrix
            ]
    ]
        ++ powerCases
        ++ structureCases
        ++ numericalCases
  where
    targets = datasetTargets dataset
    powerCases
        | entrySize entry > 500 = []
        | otherwise =
            [ env (preparedMatrix dataset) $ \ ~(Prepared matrix) ->
                bgroup
                    "power"
                    [ bench (show powerValue) $
                        nf (checksumMatrix . Matrix.power powerValue) matrix
                    | powerValue <- [2, 10, 100]
                    ]
            ]
    structureCases
        | entryFamily entry == "dense" && entrySize entry > 500 = []
        | otherwise =
            [ bgroup
                "structure"
                [ bench "classes-lifecycle" $
                    nf (classesLifecycle (Proxy @n)) (datasetRows dataset)
                , bench "classes-cold" $
                    perRunEnv (preparedGraph dataset) $ \prepared ->
                        evaluate $!
                            checksumClasses
                                (Classification.communicatingClasses (unPrepared prepared))
                , env (preparedWarmClasses dataset) $ \ ~(Prepared matrix) ->
                    bgroup
                        "classes-warm"
                        [ bench "access" $
                            whnf Classification.communicatingClasses matrix
                        , bench "consumed" $
                            nf
                                (checksumClasses . Classification.communicatingClasses)
                                matrix
                        ]
                , bench "irreducible-cold" $
                    perRunEnv (preparedGraph dataset) $ \prepared ->
                        evaluate $! Classification.irreducible (unPrepared prepared)
                , env (preparedWarmIrreducible dataset) $ \ ~(Prepared matrix) ->
                    bench "irreducible-warm" $
                        nf Classification.irreducible matrix
                ]
            ]
    numericalCases
        | entrySize entry > 500 = []
        | otherwise =
            [ bench "stationary" $
                perRunEnv (preparedGraph dataset) $ \prepared ->
                    evaluate $!
                        checksumStationary
                            ( eitherOrFail
                                (Stationary.stationaryDistributions (unPrepared prepared))
                            )
            , bench "hitting-probability" $
                perRunEnv (preparedGraph dataset) $ \prepared ->
                    evaluate $!
                        hittingProbabilityChecksum targets (unPrepared prepared)
            , bench "hitting-time" $
                perRunEnv (preparedGraph dataset) $ \prepared ->
                    evaluate $!
                        hittingExpectationChecksum targets (unPrepared prepared)
            ]

absorptionBenchmarks :: forall n. (KnownNat n) => Entry -> Dataset n -> [Benchmark]
absorptionBenchmarks entry dataset
    | entryFamily entry /= "absorbing" = []
    | entrySize entry > 500 = []
    | otherwise =
        [ bench "fundamental-matrix" $
            perRunEnv (preparedGraph dataset) $ \prepared ->
                evaluate $!
                    checksumFundamental
                        ( eitherOrFail
                            (Absorption.fundamentalMatrix (unPrepared prepared))
                        )
        , bench "absorption-time" $
            perRunEnv (preparedGraph dataset) $ \prepared ->
                evaluate $! absorptionExpectationChecksum (unPrepared prepared)
        ]

occupationBenchmarks :: forall n. (KnownNat n) => Entry -> Dataset n -> [Benchmark]
occupationBenchmarks entry dataset
    | entryFamily entry `notElem` ["absorbing", "reducible"] = []
    | entrySize entry > 250 = []
    | otherwise =
        [ bench "occupation-matrix" $
            perRunEnv (preparedGraph dataset) $ \prepared ->
                evaluate $!
                    checksumOccupation
                        ( eitherOrFail
                            (VisitCount.occupationMatrix (unPrepared prepared))
                        )
        ]

simulationBenchmarks :: forall n. (KnownNat n) => Entry -> Dataset n -> [Benchmark]
simulationBenchmarks entry dataset
    | entryFamily entry `notElem` ["dense", "low-outdegree"] = []
    | otherwise =
        [ env (preparedMatrix dataset) $ \ ~(Prepared matrix) ->
            bgroup
                "simulation"
                [ bench (show transitions) $
                    perRunEnv (preparedGenerator (entrySeed entry)) $ \generator ->
                        simulationChecksum
                            transitions
                            matrix
                            (unPreparedGen generator)
                | transitions <- [10000, 100000]
                ]
        ]

matrixFromRows :: forall n. (KnownNat n) => [[Double]] -> Matrix.TransitionMatrix (Finite n)
matrixFromRows = eitherOrFail . Matrix.fromRows

vectorFromWeights :: forall n. (KnownNat n) => [Double] -> Vector.DistributionVector (Finite n)
vectorFromWeights = eitherOrFail . Vector.fromList

preparedMatrix :: forall n. (KnownNat n) => Dataset n -> IO (Prepared n)
preparedMatrix dataset = do
    let matrix = matrixFromRows @n (datasetRows dataset)
    _ <- evaluate (checksumMatrix matrix)
    pure (Prepared matrix)

preparedDynamics :: forall n. (KnownNat n) => Dataset n -> IO (PreparedDynamics n)
preparedDynamics dataset = do
    Prepared matrix <- preparedMatrix dataset
    let initial = vectorFromWeights @n (datasetInitialWeights dataset)
    _ <- evaluate (checksumVector initial)
    pure (PreparedDynamics matrix initial)

preparedGraph :: forall n. (KnownNat n) => Dataset n -> IO (Prepared n)
preparedGraph dataset = do
    prepared <- preparedMatrix dataset
    case finiteStates :: [Finite n] of
        [] -> pure prepared
        first : _ -> do
            _ <- evaluate (Classification.supportEdge (unPrepared prepared) first first)
            pure prepared

preparedWarmClasses :: forall n. (KnownNat n) => Dataset n -> IO (Prepared n)
preparedWarmClasses dataset = do
    prepared <- preparedGraph dataset
    _ <- evaluate (checksumClasses (Classification.communicatingClasses (unPrepared prepared)))
    pure prepared

preparedWarmIrreducible :: forall n. (KnownNat n) => Dataset n -> IO (Prepared n)
preparedWarmIrreducible dataset = do
    prepared <- preparedGraph dataset
    _ <- evaluate (Classification.irreducible (unPrepared prepared))
    pure prepared

preparedGenerator :: Int -> IO PreparedGen
preparedGenerator seed =
    PreparedGen
        <$> MWC.initialize
            ( U.fromList
                [ fromIntegral seed
                , fromIntegral (seed * 1664525 + 1013904223)
                , 0x9e3779b9
                , 0x243f6a88
                ]
            )

classesLifecycle :: forall n. (KnownNat n) => Proxy n -> [[Double]] -> Double
classesLifecycle _ = checksumClasses . Classification.communicatingClasses . matrixFromRows @n

hittingProbabilityChecksum ::
    forall n.
    (KnownNat n) =>
    [Finite n] ->
    Matrix.TransitionMatrix (Finite n) ->
    Double
hittingProbabilityChecksum targets matrix =
    weightedEither
        [lookupProbability state | state <- finiteStates]
  where
    lookupProbability = Hitting.eventualProbabilityGivenInitialState matrix targets

hittingExpectationChecksum ::
    forall n.
    (KnownNat n) =>
    [Finite n] ->
    Matrix.TransitionMatrix (Finite n) ->
    Double
hittingExpectationChecksum targets matrix =
    checksumExpectations
        (eitherOrFail (sequence [lookupExpectation state | state <- finiteStates]))
  where
    lookupExpectation = Hitting.expectationGivenInitialState matrix targets

absorptionExpectationChecksum ::
    forall n.
    (KnownNat n) =>
    Matrix.TransitionMatrix (Finite n) ->
    Double
absorptionExpectationChecksum matrix =
    checksumExpectations
        (eitherOrFail (sequence [lookupExpectation state | state <- finiteStates]))
  where
    lookupExpectation = Absorption.expectationGivenInitialState matrix

weightedEither :: (Show problem) => [Either problem Double] -> Double
weightedEither =
    foldl'
        (\total value -> total + eitherOrFail value)
        0

simulationChecksum ::
    forall n.
    (KnownNat n) =>
    Natural ->
    Matrix.TransitionMatrix (Finite n) ->
    MWC.Gen RealWorld ->
    IO Double
simulationChecksum transitions matrix generator = do
    case finiteStates :: [Finite n] of
        [] -> error "benchmark simulation requires a non-empty state space"
        initial : _ -> do
            result <- Simulation.simulateMatrix transitions matrix initial generator
            pure (checksumStates (eitherOrFail result))
