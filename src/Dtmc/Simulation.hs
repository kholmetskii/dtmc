{- |
Module      : Dtmc.Simulation
Description : Sampling states and running the chain forward.

Random sampling from dense or sparse state distributions, plus shared
simulation through any locally finite 'Transition'. Failures are returned as
'SimulationError' values. A validation failure leaves the supplied MWC
generator unchanged; successfully validated sampling passes it to the
categorical backend in any 'PrimMonad'.
-}
module Dtmc.Simulation (
    SimulationError (..),
    sample,
    step,
    simulateN,
) where

import Control.Monad.Primitive (
    PrimMonad,
    PrimState,
 )
import Dtmc.Distribution (
    Distribution (..),
 )
import Dtmc.Simplex.Internal (
    simplexTolerance,
 )
import Dtmc.Transition (
    Transition (..),
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.Natural (
    Natural,
 )
import System.Random.MWC qualified as MWC
import System.Random.MWC.Distributions qualified as MWCD

{- | Why sampling could not produce a state. Weight indices refer to the order
returned by 'distributionWeights'. Input errors are detected before the random
generator is used.
-}
data SimulationError
    = -- | The distribution stores no states.
      EmptySupport
    | -- | Zero-based index of a weight that is @NaN@ or infinite.
      NonFiniteWeight Int
    | -- | Zero-based index and value of a weight below @-1e-9@.
      NegativeWeight Int Double
    | -- | Finite individual weights overflowed while being summed.
      NonFiniteTotal
    | -- | The repaired weights have a zero or negative total.
      NonPositiveTotal Double
    | -- | Impossible backend index and the stored support size.
      SampleIndexOutOfBounds Int Int
    deriving (Eq, Show)

{- | Draw a state from any finite-support 'Distribution'. Before sampling,
stored weights in @[-1e-9, 0)@ are replaced by zero; the categorical sampler
scales by the resulting total, so no explicit renormalisation is stored.

Returns 'Left' for empty support, non-finite weights or totals, weights below
@-1e-9@, or a non-positive repaired total. Validation happens before the
generator is advanced.

Time and temporary space: @O(s)@ for stored support size @s@.
-}
sample ::
    (Distribution distribution, PrimMonad m) =>
    distribution ->
    MWC.Gen (PrimState m) ->
    m (Either SimulationError (DistributionState distribution))
sample distribution generator =
    case prepareEntries (distributionWeights distribution) of
        Left problem -> pure (Left problem)
        Right (states, weights) -> do
            index <- MWCD.categorical weights generator
            pure
                ( case atMay states index of
                    Nothing -> Left (SampleIndexOutOfBounds index (length states))
                    Just state -> Right state
                )

prepareEntries :: [(state, Double)] -> Either SimulationError ([state], LA.Vector Double)
prepareEntries [] = Left EmptySupport
prepareEntries entries = do
    repaired <- traverse repairWeight (zip [0 ..] (map snd entries))
    let total = foldl' (+) 0 repaired
    validateTotal total
    pure (map fst entries, LA.fromList repaired)

validateTotal :: Double -> Either SimulationError ()
validateTotal total
    | isNaN total || isInfinite total = Left NonFiniteTotal
    | total <= 0 = Left (NonPositiveTotal total)
    | otherwise = Right ()

repairWeight :: (Int, Double) -> Either SimulationError Double
repairWeight (index, weight)
    | isNaN weight || isInfinite weight = Left (NonFiniteWeight index)
    | weight < negate simplexTolerance = Left (NegativeWeight index weight)
    | weight < 0 = Right 0
    | otherwise = Right weight

atMay :: [value] -> Int -> Maybe value
atMay _ index | index < 0 = Nothing
atMay values index =
    case drop index values of
        [] -> Nothing
        value : _ -> Just value

{- | Sample one transition from a state through any 'Transition'. Passing
each result back with the same generator advances one trajectory. The finite
transition law inherits the checked repair behavior of 'sample'.
-}
step ::
    (PrimMonad m, Transition kernel) =>
    kernel ->
    TransitionState kernel ->
    MWC.Gen (PrimState m) ->
    m (Either SimulationError (TransitionState kernel))
step kernel state =
    sample (transitionLaw kernel state)

{- | Simulate exactly @k@ transitions through any 'Transition'. On success,
return the trajectory including its initial state, with length @k + 1@. Stop
at the first invalid transition law and return its 'SimulationError'. At
@k = 0@, return the initial state without inspecting the kernel or advancing
the generator.
-}
simulateN ::
    (PrimMonad m, Transition kernel) =>
    Natural ->
    kernel ->
    TransitionState kernel ->
    MWC.Gen (PrimState m) ->
    m (Either SimulationError [TransitionState kernel])
simulateN transitions kernel initial generator =
    go transitions initial [initial]
  where
    go 0 _ reversed = pure (Right (reverse reversed))
    go remaining current reversed = do
        result <- step kernel current generator
        case result of
            Left problem -> pure (Left problem)
            Right next -> go (remaining - 1) next (next : reversed)
