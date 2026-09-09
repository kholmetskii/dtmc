{-# LANGUAGE GADTs #-}

module Main (main) where

import Cases (benchmarksFor)
import Criterion.Main (defaultMain)
import Criterion.Types (Benchmark)
import Dataset
import System.Environment (getArgs, lookupEnv)
import Text.Read (readMaybe)
import Verification (writeVerification)

main :: IO ()
main = do
    arguments <- getArgs
    dataRoot <- maybe "benchmark/data/generated" id <$> lookupEnv "DTMC_BENCH_DATA"
    manifest <- loadManifest (dataRoot ++ "/manifest.json")
    family <- lookupEnv "DTMC_BENCH_FAMILY"
    size <- readEnvironment "DTMC_BENCH_SIZE"
    seed <- readEnvironment "DTMC_BENCH_SEED"
    let selected = selectEntries family size seed manifest
    if null selected
        then fail "no datasets matched DTMC_BENCH_FAMILY/DTMC_BENCH_SIZE/DTMC_BENCH_SEED"
        else pure ()
    typed <- traverse (either fail pure . toSomeEntry) selected
    case arguments of
        ["--verify-json", outputPath] -> do
            maximumSize <- maybe 100 id <$> readEnvironment "DTMC_VERIFY_MAX_SIZE"
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
