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
import Dtmc.Analysis.Event (DiscreteEvent (AtMost, EqualTo, GreaterThan))
import Dtmc.Analysis.FiniteTime qualified as FiniteTime
import Dtmc.Analysis.HittingTime qualified as Hitting
import Dtmc.Analysis.ReturnTime qualified as Return
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

data PreparedScalarLookup n = PreparedScalarLookup
    { scalarLookup :: Finite n -> Double
    , scalarLookupState :: Finite n
    }

instance NFData (PreparedScalarLookup n) where
    rnf prepared =
        scalarLookup prepared
            `seq` scalarLookupState prepared
            `seq` ()

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
                ++ finiteTimeBenchmarks entry dataset
                ++ boundedHittingBenchmarks entry dataset
                ++ absorptionBenchmarks entry dataset
                ++ occupationBenchmarks entry dataset
                ++ returnBenchmarks entry dataset
                ++ visitBenchmarks entry dataset
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
    competing = datasetCompeting dataset
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
                ( [ bench "classes-lifecycle" $
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
                    ++ periodicStructureCases
                )
            ]
      where
        periodicStructureCases
            | entryFamily entry /= "periodic" = []
            | otherwise =
                [ bench "period-cold" $
                    perRunEnv (preparedGraph dataset) $ \prepared ->
                        evaluate $!
                            checksumPeriod
                                (Classification.chainPeriod (unPrepared prepared))
                , env (preparedWarmPeriod dataset) $ \ ~(Prepared matrix) ->
                    bench "period-warm" $
                        nf (checksumPeriod . Classification.chainPeriod) matrix
                , bench "cyclic-classes-cold" $
                    perRunEnv (preparedGraph dataset) $ \prepared ->
                        evaluate $!
                            checksumCyclicClasses
                                (Classification.cyclicClasses (unPrepared prepared))
                , env (preparedWarmCyclicClasses dataset) $ \ ~(Prepared matrix) ->
                    bench "cyclic-classes-warm" $
                        nf
                            (checksumCyclicClasses . Classification.cyclicClasses)
                            matrix
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
            , bgroup
                "hitting-probability"
                [ bench "cold-all-states" $
                    perRunEnv (preparedGraph dataset) $ \prepared ->
                        evaluate $!
                            hittingProbabilityChecksum targets (unPrepared prepared)
                , env (preparedWarmHittingProbability dataset) $ \prepared ->
                    bench "warm-lookup" $
                        whnf
                            (\fixture -> scalarLookup fixture (scalarLookupState fixture))
                            prepared
                ]
            , bgroup
                "hitting-time"
                [ bench "cold-all-states" $
                    perRunEnv (preparedGraph dataset) $ \prepared ->
                        evaluate $!
                            hittingExpectationChecksum targets (unPrepared prepared)
                , env (preparedWarmHittingExpectation dataset) $ \prepared ->
                    bench "warm-lookup" $
                        whnf
                            (\fixture -> scalarLookup fixture (scalarLookupState fixture))
                            prepared
                ]
            ]
                ++ raceCases
      where
        raceCases
            | entryFamily entry `notElem` ["dense", "low-outdegree"] = []
            | otherwise =
                [ bench "race/forward-committor" $
                    perRunEnv (preparedGraph dataset) $ \prepared ->
                        evaluate $!
                            raceProbabilityChecksum
                                targets
                                competing
                                (unPrepared prepared)
                ]

finiteTimeBenchmarks :: forall n. (KnownNat n) => Entry -> Dataset n -> [Benchmark]
finiteTimeBenchmarks entry dataset
    | entrySize entry > 100 = []
    | otherwise =
        [ bench "finite-time/step" $
            perRunEnv (preparedMatrix dataset) $ \prepared ->
                evaluate $!
                    FiniteTime.stepProbability
                        (unPrepared prepared)
                        source
                        target
        ]
            ++ [ bench ("finite-time/n-step/k-" ++ show steps) $
                    perRunEnv (preparedMatrix dataset) $ \prepared ->
                        evaluate $!
                            FiniteTime.nStepProbability
                                steps
                                (unPrepared prepared)
                                source
                                target
               | steps <- [10, 100]
               ]
            ++ [ bench ("finite-time/observation/k-" ++ show time) $
                    perRunEnv (preparedDynamics dataset) $ \prepared ->
                        evaluate $!
                            FiniteTime.probability
                                (dynamicsInitial prepared)
                                (dynamicsMatrix prepared)
                                [FiniteTime.At time target]
               | time <- [10, 100]
               ]
  where
    source = benchmarkState @n
    target = benchmarkTarget dataset

boundedHittingBenchmarks :: forall n. (KnownNat n) => Entry -> Dataset n -> [Benchmark]
boundedHittingBenchmarks entry dataset
    | entrySize entry > 100 = []
    | otherwise =
        [ bench ("hitting/bounded/" ++ label ++ "/k-" ++ show horizon) $
            perRunEnv (preparedMatrix dataset) $ \prepared ->
                evaluate $!
                    Hitting.probabilityGivenInitialState
                        (event horizon)
                        (unPrepared prepared)
                        isTarget
                        initial
        | (label, event) <-
            [ ("exact", EqualTo)
            , ("at-most", AtMost)
            , ("greater-than", GreaterThan)
            ]
        , horizon <- [10, 100]
        ]
  where
    initial = benchmarkState @n
    targets = datasetTargets dataset
    isTarget state = state `elem` targets

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
        , bench "absorption/probabilities" $
            perRunEnv (preparedGraph dataset) $ \prepared ->
                evaluate $!
                    absorptionProbabilityChecksum (unPrepared prepared)
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

returnBenchmarks :: forall n. (KnownNat n) => Entry -> Dataset n -> [Benchmark]
returnBenchmarks entry dataset = meanRecurrenceCase ++ boundedCases
  where
    meanRecurrenceCase
        | entryFamily entry `notElem` ["dense", "low-outdegree"] = []
        | entrySize entry > 500 = []
        | otherwise =
            [ bench "return/mean-recurrence" $
                perRunEnv (preparedGraph dataset) $ \prepared ->
                    evaluate $! returnExpectationChecksum (unPrepared prepared)
            ]
    boundedCases
        | entrySize entry > 100 = []
        | otherwise =
            [ bench ("return/bounded/k-" ++ show bound) $
                perRunEnv (preparedMatrix dataset) $ \prepared ->
                    evaluate $!
                        Return.probabilityGivenInitialState
                            (AtMost bound)
                            (unPrepared prepared)
                            (benchmarkState @n)
            | bound <- [10, 100]
            ]

visitBenchmarks :: forall n. (KnownNat n) => Entry -> Dataset n -> [Benchmark]
visitBenchmarks entry dataset
    | entrySize entry > 100 = []
    | otherwise =
        [ bench ("visits/bounded-expectation/k-" ++ show bound) $
            perRunEnv (preparedMatrix dataset) $ \prepared ->
                evaluate $!
                    VisitCount.boundedExpectationGivenInitialState
                        bound
                        (benchmarkState @n)
                        (unPrepared prepared)
                        (== target)
        | bound <- [10, 100]
        ]
  where
    target = benchmarkTarget dataset

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

preparedWarmPeriod :: forall n. (KnownNat n) => Dataset n -> IO (Prepared n)
preparedWarmPeriod dataset = do
    prepared <- preparedGraph dataset
    _ <- evaluate (checksumPeriod (Classification.chainPeriod (unPrepared prepared)))
    pure prepared

preparedWarmCyclicClasses :: forall n. (KnownNat n) => Dataset n -> IO (Prepared n)
preparedWarmCyclicClasses dataset = do
    prepared <- preparedGraph dataset
    _ <- evaluate (checksumCyclicClasses (Classification.cyclicClasses (unPrepared prepared)))
    pure prepared

preparedWarmHittingProbability ::
    forall n.
    (KnownNat n) =>
    Dataset n ->
    IO (PreparedScalarLookup n)
preparedWarmHittingProbability dataset = do
    Prepared matrix <- preparedGraph dataset
    let rawLookup = Hitting.eventualProbabilityGivenInitialState matrix (datasetTargets dataset)
        lookupValue = eitherOrFail . rawLookup
        state = benchmarkState @n
    _ <- evaluate (lookupValue state)
    pure (PreparedScalarLookup lookupValue state)

preparedWarmHittingExpectation ::
    forall n.
    (KnownNat n) =>
    Dataset n ->
    IO (PreparedScalarLookup n)
preparedWarmHittingExpectation dataset = do
    Prepared matrix <- preparedGraph dataset
    let rawLookup = Hitting.expectationGivenInitialState matrix (datasetTargets dataset)
        lookupValue queryState = checksumExpectations [eitherOrFail (rawLookup queryState)]
        state = benchmarkState @n
    _ <- evaluate (lookupValue state)
    pure (PreparedScalarLookup lookupValue state)

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

raceProbabilityChecksum ::
    forall n.
    (KnownNat n) =>
    [Finite n] ->
    [Finite n] ->
    Matrix.TransitionMatrix (Finite n) ->
    Double
raceProbabilityChecksum successful competing matrix =
    weightedEither
        [lookupProbability state | state <- finiteStates]
  where
    lookupProbability =
        Hitting.raceProbabilityGivenInitialState matrix successful competing

returnExpectationChecksum ::
    forall n.
    (KnownNat n) =>
    Matrix.TransitionMatrix (Finite n) ->
    Double
returnExpectationChecksum matrix =
    checksumExpectations
        (eitherOrFail (sequence [lookupExpectation state | state <- finiteStates]))
  where
    lookupExpectation = Return.expectationGivenInitialState matrix

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

absorptionProbabilityChecksum ::
    forall n.
    (KnownNat n) =>
    Matrix.TransitionMatrix (Finite n) ->
    Double
absorptionProbabilityChecksum =
    checksumAbsorption . eitherOrFail . Absorption.probabilityMatrix

benchmarkState :: forall n. (KnownNat n) => Finite n
benchmarkState =
    case finiteStates of
        [] -> error "benchmark requires a non-empty state space"
        first : _ -> first

benchmarkTarget :: forall n. (KnownNat n) => Dataset n -> Finite n
benchmarkTarget dataset =
    case datasetTargets dataset of
        [] -> error "benchmark dataset requires at least one target"
        target : _ -> target

checksumPeriod :: Maybe Natural -> Double
checksumPeriod = maybe 0 fromIntegral

checksumCyclicClasses :: forall n. Maybe [[Finite n]] -> Double
checksumCyclicClasses = maybe 0 (checksumStates . concat)

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
