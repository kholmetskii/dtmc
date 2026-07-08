{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.SimulationSpec where

import Control.Monad (replicateM)
import Data.Either (fromRight)
import Data.Finite (Finite, finite)
import Dtmc.Simulation (step)
import Dtmc.StochasticMatrix (StochasticMatrix, mkStochasticMatrix)
import Numeric.LinearAlgebra (fromLists)
import Test.Hspec
  ( Spec
  , describe
  , it
  , shouldBe
  , shouldSatisfy
  )
import qualified System.Random.MWC as MWC

spec :: Spec
spec = do
  describe "step" $ do
    it "always stays in state 0 for an absorbing state 0" $ do
      gen <- MWC.create

      let matrix = absorbingMatrix
          state0 = finite @2 0

      nextStates <- replicateM 20 (step matrix state0 gen)

      nextStates `shouldSatisfy` all (== state0)

    it "always stays in state 1 for an absorbing state 1" $ do
      gen <- MWC.create

      let matrix = absorbingMatrix
          state1 = finite @2 1

      nextStates <- replicateM 20 (step matrix state1 gen)

      nextStates `shouldSatisfy` all (== state1)

absorbingMatrix :: StochasticMatrix 2
absorbingMatrix =
  fromRight
    (error "absorbingMatrix should be stochastic")
    ( mkStochasticMatrix @2
        ( fromLists
            [ [1.0, 0.0]
            , [0.0, 1.0]
            ]
        )
    )