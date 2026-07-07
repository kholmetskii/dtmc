module Main where

import Test.Hspec ( hspec )

import qualified Dtmc.StochasticMatrixSpec

main :: IO ()
main =
  hspec $ do
    Dtmc.StochasticMatrixSpec.spec