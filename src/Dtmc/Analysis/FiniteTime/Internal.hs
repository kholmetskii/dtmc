{- |
Module      : Dtmc.Analysis.FiniteTime.Internal
Description : Normalised timed observations (unsafe underbelly).

The private normal form behind the event and conditional probability queries in
"Dtmc.Analysis.FiniteTime". 'normalise' is the intended way to build a
'NormalisedObservations': it establishes the invariant that a 'Consistent' list
holds exactly one @(time, state)@ entry per distinct time, in ascending time
order. Building 'Consistent' directly can break that invariant and give the
scoring in "Dtmc.Analysis.FiniteTime" a wrong answer.
-}
module Dtmc.Analysis.FiniteTime.Internal (
    NormalisedObservations (..),
    normalise,
) where

import Data.List (
    sortBy,
 )
import Data.Ord (
    comparing,
 )
import Numeric.Natural (
    Natural,
 )

{- | A conjunction of timed state observations after sorting, de-duplication,
and consistency checking.
-}
data NormalisedObservations state
    = -- | Two observations demand different states at one time.
      Impossible
    | -- | Distinct times in ascending order, each with one required state.
      Consistent [(Natural, state)]

{- | Sort @(time, state)@ pairs by ascending time, collapse exact duplicates,
and detect contradictions. Pairs requiring different states at the same time
yield 'Impossible'; otherwise the result is 'Consistent' with one entry per
distinct time in ascending order, so consecutive entries always have strictly
increasing times.

Sorting is @O(m log m)@ for @m@ pairs; the collapsing fold is @O(m)@.
-}
normalise :: (Eq state) => [(Natural, state)] -> NormalisedObservations state
normalise pairs =
    foldr insert (Consistent []) (sortBy (comparing fst) pairs)
  where
    insert _ Impossible = Impossible
    insert step (Consistent []) = Consistent [step]
    insert (t, i) (Consistent ((t', i') : rest))
        | t == t' && i == i' = Consistent ((t', i') : rest)
        | t == t' = Impossible
        | otherwise = Consistent ((t, i) : (t', i') : rest)
