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
    simulate,
    simulateMatrix,
) where

import Control.Monad.Primitive (
    PrimMonad,
    PrimState,
 )
import Data.Array qualified as Array
import Data.Array.Unboxed qualified as Unboxed
import Data.List qualified as List
import Dtmc.Distribution (
    Distribution (..),
 )
import Dtmc.Simplex.Internal (
    simplexTolerance,
 )
import Dtmc.State (
    FiniteState,
 )
import Dtmc.State.Internal (
    stateFromInt,
    stateIndexInt,
 )
import Dtmc.Transition (
    Transition (..),
 )
import Dtmc.Transition.Matrix.Internal (
    TransitionMatrix,
    unTransitionMatrix,
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

Complexity: excluding 'distributionWeights', @O(s + 1)@ time and @O(s)@
temporary space for stored support size @s@; result space is @O(1)@.
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
    let total = List.foldl' (+) 0 repaired
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

{- | Sample one transition from a state through any 'Transition'. Passing each
result back with the same generator advances one trajectory. The returned
finite-support law inherits the checked repair behaviour of 'sample'.

Complexity: excluding 'transitionLaw' and 'distributionWeights', @O(s + 1)@
time and @O(s)@ temporary space for stored support size @s@; result space is
@O(1)@.
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

Let @s@ bound the stored support size of every transition law encountered.

Complexity: excluding 'transitionLaw' and 'distributionWeights',
@O(k (s + 1) + 1)@ time, @O(k + s + 1)@ temporary space, and @O(k + 1)@
result space.
-}
simulate ::
    (PrimMonad m, Transition kernel) =>
    Natural ->
    kernel ->
    TransitionState kernel ->
    MWC.Gen (PrimState m) ->
    m (Either SimulationError [TransitionState kernel])
simulate transitions kernel initial generator =
    go transitions initial [initial]
  where
    go 0 _ reversed = pure (Right (reverse reversed))
    go remaining current reversed = do
        result <- step kernel current generator
        case result of
            Left problem -> pure (Left problem)
            Right next -> go (remaining - 1) next (next : reversed)

-- | A validated cumulative row used by 'simulateMatrix'. The final positive
-- index is retained as a defensive fallback if the random generator returns
-- the upper endpoint of its requested floating-point interval.
data PreparedMatrixRow = PreparedMatrixRow
    !(Unboxed.UArray Int Int)
    !(LA.Vector Double)
    !Double
    !Int

{- | Simulate a finite transition matrix without converting every visited row
to a map-backed distribution. Each distinct row is validated and converted to
a shared cumulative vector on first use; sampling that row thereafter uses a
binary search.

This has the same validation and generator-advancement contract as 'simulate':
an invalid visited row is reported before drawing a random number, unvisited
rows are not inspected, and a zero-step simulation inspects neither the matrix
nor the generator.

For matrix dimension @n@, @k@ transitions, and @r@ distinct visited rows,
complexity is @O(r n + k log n)@ time, @O(r n + n + k)@ temporary space, and
@O(k)@ result space. The row cache lives only for this simulation call.
-}
simulateMatrix ::
    forall state m.
    (FiniteState state, PrimMonad m) =>
    Natural ->
    TransitionMatrix state ->
    state ->
    MWC.Gen (PrimState m) ->
    m (Either SimulationError [state])
simulateMatrix 0 _ initial _ = pure (Right [initial])
simulateMatrix transitions matrix initial generator =
    go transitions initial [initial]
  where
    stored = unTransitionMatrix matrix
    dimension = LA.rows stored

    -- 'Array' is lazy in its elements, so only rows reached by the trajectory
    -- are converted and validated. Forced entries are then shared.
    preparedRows :: Array.Array Int (Either SimulationError PreparedMatrixRow)
    preparedRows =
        Array.listArray
            (0, dimension - 1)
            [prepareMatrixRow (LA.toList row) | row <- LA.toRows stored]

    go 0 _ reversed = pure (Right (reverse reversed))
    go remaining current reversed =
        case preparedRows Array.! stateIndexInt current of
            Left problem -> pure (Left problem)
            Right prepared -> do
                index <- samplePreparedMatrixRow prepared generator
                case stateFromInt index of
                    Nothing ->
                        pure
                            (Left (SampleIndexOutOfBounds index dimension))
                    Just next ->
                        go (remaining - 1) next (next : reversed)

prepareMatrixRow :: [Double] -> Either SimulationError PreparedMatrixRow
prepareMatrixRow weights = do
    let entries = filter ((/= 0) . snd) (zip [0 ..] weights)
    if null entries then Left EmptySupport else pure ()
    repaired <- traverse repairWeight (zip [0 ..] (map snd entries))
    let cumulative = drop 1 (List.scanl' (+) 0 repaired)
        total = List.foldl' (+) 0 repaired
    validateTotal total
    let lastPositive =
            List.foldl'
                (\latest (index, weight) -> if weight > 0 then index else latest)
                0
                (zip [0 ..] repaired)
        outcomes =
            Unboxed.listArray
                (0, length entries - 1)
                (map fst entries)
    pure
        ( PreparedMatrixRow
            outcomes
            (LA.fromList cumulative)
            total
            lastPositive
        )

samplePreparedMatrixRow ::
    (PrimMonad m) =>
    PreparedMatrixRow ->
    MWC.Gen (PrimState m) ->
    m Int
samplePreparedMatrixRow
    (PreparedMatrixRow outcomes cumulative total lastPositive)
    generator = do
    target <- MWC.uniformR (0, total) generator
    pure (outcomes Unboxed.! firstGreater target cumulative lastPositive)

firstGreater :: Double -> LA.Vector Double -> Int -> Int
firstGreater target cumulative fallback = search 0 (LA.size cumulative)
  where
    search lower upper
        | lower >= upper =
            if lower < LA.size cumulative then lower else fallback
        | cumulative `LA.atIndex` middle > target = search lower middle
        | otherwise = search (middle + 1) upper
      where
        middle = lower + (upper - lower) `quot` 2
