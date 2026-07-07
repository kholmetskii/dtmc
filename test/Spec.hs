module Main where

import Test.Hspec ( hspec )

import qualified Dtmc.StochasticSpec

main :: IO ()
main =
  hspec $ do
    Dtmc.StochasticSpec.spec