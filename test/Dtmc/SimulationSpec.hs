module Dtmc.SimulationSpec (spec) where

import Control.Monad (replicateM)
import Control.Monad.ST (runST)
import Data.Finite (Finite)
import Dtmc.Simulation (step)
import Dtmc.TransitionMatrix (TransitionMatrix, mkTransitionMatrix)
import qualified Numeric.LinearAlgebra.Static as S
import qualified System.Random.MWC as MWC
import Test.Hspec ( Spec, it, shouldBe )

cyclic3 :: TransitionMatrix 3
cyclic3 =
  either (error . show) id $
    mkTransitionMatrix (S.matrix [0, 1, 0, 0, 0, 1, 1, 0, 0])

-- Deliberately asymmetric: state 0 is absorbing, state 1 is not.
-- (The identity matrix would be symmetric and useless as an oracle.)
absorbing2 :: TransitionMatrix 2
absorbing2 =
  either (error . show) id $
    mkTransitionMatrix (S.matrix [1, 0, 0.3, 0.7])

-- Everything runs in ST with a fixed seed (MWC.create), so the tests are
-- DETERMINISTIC (docs/TESTING.md T5). This is also the only place where step's
-- PrimMonad polymorphism is actually exercised: production path IO, test path ST.
orbit3 :: [Finite 3]
orbit3 = runST $ do
  g <- MWC.create
  s1 <- step cyclic3 0 g
  s2 <- step cyclic3 s1 g
  s3 <- step cyclic3 s2 g
  pure [s1, s2, s3]

absorbed :: [Finite 2]
absorbed = runST $ do
  g <- MWC.create
  replicateM 50 (step absorbing2 0 g)

spec :: Spec
spec = do
  it "follows the 3-cycle 0 -> 1 -> 2 -> 0 (transposition would give 0 -> 2 -> 1)" $
    orbit3 `shouldBe` [1, 2, 0]

  it "never leaves an absorbing state" $
    absorbed `shouldBe` replicate 50 0

  it "runs in ST, exercising the PrimMonad polymorphism of step" $
    length orbit3 `shouldBe` 3