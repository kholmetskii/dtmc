{-# LANGUAGE GADTs #-}

module Main (main) where

import Cases (benchmarksFor)
import Control.Monad (when)
import Criterion.Main (defaultMain)
import Criterion.Types (Benchmark)
import Data.Maybe (fromMaybe)
import Dataset
import System.Environment (getArgs, lookupEnv)
import Text.Read (readMaybe)
import Verification (writeVerification)

main :: IO ()
main = do
    arguments <- getArgs
    dataRoot <- fromMaybe "benchmark/data/generated" <$> lookupEnv "DTMC_BENCH_DATA"
    manifest <- loadManifest (dataRoot ++ "/manifest.json")
    family <- lookupEnv "DTMC_BENCH_FAMILY"
    size <- readEnvironment "DTMC_BENCH_SIZE"
    seed <- readEnvironment "DTMC_BENCH_SEED"
    let selected = selectEntries family size seed manifest
    when (null selected) $
        fail "no datasets matched DTMC_BENCH_FAMILY/DTMC_BENCH_SIZE/DTMC_BENCH_SEED"
    typed <- traverse (either fail pure . toSomeEntry) selected
    case arguments of
        ["--verify-json", outputPath] -> do
            maximumSize <- fromMaybe 100 <$> readEnvironment "DTMC_VERIFY_MAX_SIZE"
            writeVerification outputPath dataRoot maximumSize typed
        _ -> defaultMain (map (benchmarkFor dataRoot) typed)

benchmarkFor :: FilePath -> SomeEntry -> Benchmark
benchmarkFor dataRoot (SomeEntry proxy entry) =
    benchmarksFor proxy dataRoot entry

readEnvironment :: (Read value) => String -> IO (Maybe value)
readEnvironment name = do
    raw <- lookupEnv name
    case raw of
        Nothing -> pure Nothing
        Just value ->
            case readMaybe value of
                Nothing -> fail ("environment variable " ++ name ++ " is invalid: " ++ value)
                Just parsed -> pure (Just parsed)
