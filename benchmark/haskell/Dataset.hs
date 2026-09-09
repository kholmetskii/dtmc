{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

module Dataset (
    Dataset (..),
    Entry (..),
    Manifest (..),
    SomeEntry (..),
    datasetTargets,
    loadDataset,
    loadManifest,
    selectEntries,
    toSomeEntry,
) where

import Control.DeepSeq (NFData)
import Control.Monad (replicateM, unless)
import Data.Aeson (FromJSON (parseJSON), eitherDecodeFileStrict', withObject, (.:))
import Data.Binary.Get (Get, getDoublele, isEmpty, runGetOrFail)
import Data.ByteString.Lazy qualified as BL
import Data.Finite (Finite)
import Data.List (sortOn)
import Data.Proxy (Proxy (Proxy))
import Dtmc.State (finiteStates)
import GHC.Generics (Generic)
import GHC.TypeNats (KnownNat)
import System.FilePath ((</>))

data Entry = Entry
    { entryId :: String
    , entryFamily :: String
    , entrySize :: Int
    , entrySeed :: Int
    , entryFile :: FilePath
    , entrySha256 :: String
    , entryTargetIndices :: [Int]
    , entryAbsorbingIndices :: [Int]
    }
    deriving stock (Eq, Show, Generic)
    deriving anyclass (NFData)

instance FromJSON Entry where
    parseJSON = withObject "benchmark dataset" $ \value ->
        Entry
            <$> value .: "id"
            <*> value .: "family"
            <*> value .: "size"
            <*> value .: "seed"
            <*> value .: "file"
            <*> value .: "sha256"
            <*> value .: "targets"
            <*> value .: "absorbing"

newtype Manifest = Manifest {manifestEntries :: [Entry]}
    deriving stock (Eq, Show)

instance FromJSON Manifest where
    parseJSON = withObject "benchmark manifest" $ \value ->
        Manifest <$> value .: "datasets"

data Dataset n = Dataset
    { datasetEntry :: Entry
    , datasetRows :: [[Double]]
    , datasetInitialWeights :: [Double]
    }
    deriving stock (Generic)
    deriving anyclass (NFData)

data SomeEntry where
    SomeEntry :: (KnownNat n) => Proxy n -> Entry -> SomeEntry

loadManifest :: FilePath -> IO Manifest
loadManifest path = do
    decoded <- eitherDecodeFileStrict' path
    case decoded of
        Left problem -> fail ("cannot decode benchmark manifest: " ++ problem)
        Right manifest -> pure manifest

selectEntries ::
    Maybe String ->
    Maybe Int ->
    Maybe Int ->
    Manifest ->
    [Entry]
selectEntries family size seed =
    sortOn (\entry -> (entryFamily entry, entrySize entry, entrySeed entry))
        . filter matches
        . manifestEntries
  where
    matches entry =
        maybe True (== entryFamily entry) family
            && maybe True (== entrySize entry) size
            && maybe True (== entrySeed entry) seed

toSomeEntry :: Entry -> Either String SomeEntry
toSomeEntry entry =
    case entrySize entry of
        10 -> Right (SomeEntry (Proxy @10) entry)
        25 -> Right (SomeEntry (Proxy @25) entry)
        50 -> Right (SomeEntry (Proxy @50) entry)
        100 -> Right (SomeEntry (Proxy @100) entry)
        250 -> Right (SomeEntry (Proxy @250) entry)
        500 -> Right (SomeEntry (Proxy @500) entry)
        1000 -> Right (SomeEntry (Proxy @1000) entry)
        unsupported -> Left ("unsupported type-level benchmark size: " ++ show unsupported)

loadDataset :: forall n. FilePath -> Entry -> IO (Dataset n)
loadDataset dataRoot entry = do
    let size = entrySize entry
        expected = size * size + size
        path = dataRoot </> entryFile entry
    bytes <- BL.readFile path
    values <-
        case runGetOrFail (payload expected) bytes of
            Left (_, offset, problem) ->
                fail ("cannot decode " ++ path ++ " at byte " ++ show offset ++ ": " ++ problem)
            Right (_, _, decoded) -> pure decoded
    let (matrixValues, initial) = splitAt (size * size) values
        rows = chunksOf size matrixValues
    unless (length rows == size && all ((== size) . length) rows) $
        fail ("invalid matrix dimensions in " ++ path)
    pure
        Dataset
            { datasetEntry = entry
            , datasetRows = rows
            , datasetInitialWeights = initial
            }

datasetTargets :: forall n. (KnownNat n) => Dataset n -> [Finite n]
datasetTargets dataset = map stateAtIndex (entryTargetIndices (datasetEntry dataset))
  where
    states = finiteStates :: [Finite n]
    stateAtIndex index
        | index < 0 || index >= length states =
            error ("target index outside state space: " ++ show index)
        | otherwise = states !! index

payload :: Int -> Get [Double]
payload count = do
    values <- replicateM count getDoublele
    empty <- isEmpty
    unless empty (fail "trailing bytes")
    pure values

chunksOf :: Int -> [value] -> [[value]]
chunksOf _ [] = []
chunksOf width values =
    let (prefix, suffix) = splitAt width values
     in prefix : chunksOf width suffix
