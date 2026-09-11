{-# LANGUAGE TypeFamilies #-}

{- |
Module      : Dtmc.Analysis.FiniteHorizon.Internal
Description : Dense finite-state primitives for finite-horizon analyses.

Private BLAS-backed recurrences shared by finite-time, hitting-time,
return-time, and bounded-visit analyses. Public analysis modules retain their
locally finite map recurrences and select these implementations only when a
transition representation supplies a suitable dense backend.
-}
module Dtmc.Analysis.FiniteHorizon.Internal (
    availableDenseBackend,
    branchingDenseBackend,
    denseStepProbability,
    denseProbabilityAfter,
    denseExactHittingProbability,
    denseLowerHittingProbability,
    denseUpperHittingProbability,
    denseExactReturnProbability,
    denseLowerReturnProbability,
    denseUpperReturnProbability,
    denseBoundedVisitExpectation,
) where

import Data.Map.Strict qualified as Map
import Data.Vector.Storable qualified as Storable
import Data.Vector.Storable.Mutable qualified as Mutable
import Dtmc.Transition.Internal (
    DenseTransitionBackend (..),
    Transition (..),
    denseRowBranches,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.Natural (
    Natural,
 )

availableDenseBackend ::
    (Transition transition) =>
    transition ->
    Maybe (DenseTransitionBackend (TransitionState transition))
availableDenseBackend = transitionDenseBackend

branchingDenseBackend ::
    (Transition transition) =>
    transition ->
    [TransitionState transition] ->
    Maybe (DenseTransitionBackend (TransitionState transition))
branchingDenseBackend transition sources = do
    backend <- availableDenseBackend transition
    if any (denseRowBranches backend) sources then Just backend else Nothing

denseStepProbability ::
    DenseTransitionBackend state ->
    state ->
    state ->
    Double
denseStepProbability backend source destination =
    denseTransitionMatrix backend
        `LA.atIndex` ( denseTransitionIndex backend source
                     , denseTransitionIndex backend destination
                     )

denseProbabilityAfter ::
    Natural ->
    [(state, Double)] ->
    DenseTransitionBackend state ->
    state ->
    Double
denseProbabilityAfter steps weights backend destination =
    evolved Storable.! denseTransitionIndex backend destination
  where
    evolved = iterateDense steps backend (weightsVector backend weights)

denseExactHittingProbability ::
    Natural ->
    DenseTransitionBackend state ->
    (state -> Bool) ->
    state ->
    Double
denseExactHittingProbability time backend isTarget initial =
    go time (pointVector backend initial)
  where
    (targetMask, survivorMask) = hittingMasks backend initial isTarget

    go 0 _ = 0
    go remaining survivors =
        let advanced = advanceDense backend survivors
            hitMass = targetMask LA.<.> advanced
         in if remaining == 1
                then hitMass
                else
                    let next = Storable.zipWith (*) survivorMask advanced
                     in next `seq`
                            if Storable.all (== 0) next
                                then 0
                                else go (remaining - 1) next

denseLowerHittingProbability ::
    Natural ->
    DenseTransitionBackend state ->
    (state -> Bool) ->
    state ->
    Double
denseLowerHittingProbability steps backend isTarget initial =
    go steps (pointVector backend initial) 0
  where
    (targetMask, survivorMask) = hittingMasks backend initial isTarget

    go 0 _ total = total
    go remaining survivors total =
        let advanced = advanceDense backend survivors
            hitMass = targetMask LA.<.> advanced
            cumulative = total + hitMass
            next = Storable.zipWith (*) survivorMask advanced
         in cumulative `seq`
                next `seq`
                    if Storable.all (== 0) next
                        then cumulative
                        else go (remaining - 1) next cumulative

denseUpperHittingProbability ::
    Natural ->
    DenseTransitionBackend state ->
    (state -> Bool) ->
    state ->
    Double
denseUpperHittingProbability time backend isTarget initial =
    go time (pointVector backend initial)
  where
    (_, survivorMask) = hittingMasks backend initial isTarget

    go 0 survivors = Storable.sum survivors
    go remaining survivors =
        let advanced = advanceDense backend survivors
            next = Storable.zipWith (*) survivorMask advanced
         in next `seq`
                if Storable.all (== 0) next
                    then 0
                    else go (remaining - 1) next

denseExactReturnProbability ::
    Natural ->
    DenseTransitionBackend state ->
    state ->
    Double
denseExactReturnProbability time backend initial =
    go time (pointVector backend initial)
  where
    initialIndex = denseTransitionIndex backend initial

    go 0 _ = 0
    go remaining survivors =
        let advanced = advanceDense backend survivors
            returnMass = advanced Storable.! initialIndex
         in if remaining == 1
                then returnMass
                else go (remaining - 1) (clearCoordinate initialIndex advanced)

denseLowerReturnProbability ::
    Natural ->
    DenseTransitionBackend state ->
    state ->
    Double
denseLowerReturnProbability bound backend initial =
    go bound (pointVector backend initial) 0
  where
    initialIndex = denseTransitionIndex backend initial

    go remaining _ total | remaining <= 1 = total
    go remaining survivors total =
        let advanced = advanceDense backend survivors
            returnMass = advanced Storable.! initialIndex
            cumulative = total + returnMass
            next = clearCoordinate initialIndex advanced
         in cumulative `seq` next `seq` go (remaining - 1) next cumulative

denseUpperReturnProbability ::
    Natural ->
    DenseTransitionBackend state ->
    state ->
    Double
denseUpperReturnProbability time backend initial =
    go time (pointVector backend initial)
  where
    initialIndex = denseTransitionIndex backend initial

    go 0 survivors = Storable.sum survivors
    go remaining survivors =
        let advanced = advanceDense backend survivors
            next = clearCoordinate initialIndex advanced
         in next `seq` go (remaining - 1) next

denseBoundedVisitExpectation ::
    Natural ->
    Map.Map state Double ->
    DenseTransitionBackend state ->
    (state -> Bool) ->
    Double
denseBoundedVisitExpectation 0 _ _ _ = 0
denseBoundedVisitExpectation bound weights backend isVisited =
    go bound (weightsVector backend (Map.toList weights)) 0
  where
    visitedMask =
        Storable.fromList
            [if isVisited state then 1 else 0 | state <- denseTransitionStates backend]

    go remaining current expectation =
        let visitProbability = visitedMask LA.<.> current
            cumulative = expectation + visitProbability
         in if remaining == 1
                then cumulative
                else
                    let next = advanceDense backend current
                     in cumulative `seq` next `seq` go (remaining - 1) next cumulative

advanceDense ::
    DenseTransitionBackend state ->
    LA.Vector Double ->
    LA.Vector Double
advanceDense backend = (LA.tr (denseTransitionMatrix backend) LA.#>)

iterateDense ::
    Natural ->
    DenseTransitionBackend state ->
    LA.Vector Double ->
    LA.Vector Double
iterateDense 0 _ weights = weights
iterateDense remaining backend weights =
    let next = advanceDense backend weights
     in next `seq` iterateDense (remaining - 1) backend next

pointVector :: DenseTransitionBackend state -> state -> LA.Vector Double
pointVector backend selected =
    Storable.generate dimension $ \index ->
        if index == selectedIndex then 1 else 0
  where
    dimension = LA.rows (denseTransitionMatrix backend)
    selectedIndex = denseTransitionIndex backend selected

weightsVector ::
    DenseTransitionBackend state ->
    [(state, Double)] ->
    LA.Vector Double
weightsVector backend weights =
    Storable.accum
        (+)
        (Storable.replicate dimension 0)
        [ (denseTransitionIndex backend state, weight)
        | (state, weight) <- weights
        , weight /= 0
        ]
  where
    dimension = LA.rows (denseTransitionMatrix backend)

hittingMasks ::
    DenseTransitionBackend state ->
    state ->
    (state -> Bool) ->
    (LA.Vector Double, LA.Vector Double)
hittingMasks backend initial isTarget =
    (targetMask, Storable.map (1 -) targetMask)
  where
    initialIndex = denseTransitionIndex backend initial
    targetMask =
        Storable.fromList
            [ if index /= initialIndex && isTarget state then 1 else 0
            | (index, state) <- zip [0 ..] (denseTransitionStates backend)
            ]

clearCoordinate :: Int -> LA.Vector Double -> LA.Vector Double
clearCoordinate index =
    Storable.modify $ \mutable -> Mutable.unsafeWrite mutable index 0
