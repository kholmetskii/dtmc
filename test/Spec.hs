module Main where

import Test.Hspec ( hspec )

import qualified Dtmc.StochasticMatrixSpec
import qualified Dtmc.SimulationSpec

main :: IO ()
main =
  hspec $ do
    Dtmc.StochasticMatrixSpec.spec
    Dtmc.SimulationSpec.spec