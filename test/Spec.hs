module Main where

import Test.Hspec

import qualified Dtmc.StochasticSpec

main :: IO ()
main =
    hspec $ do
        Dtmc.StochasticSpec.spec