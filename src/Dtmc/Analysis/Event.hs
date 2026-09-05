{- |
Module      : Dtmc.Analysis.Event
Description : Closed comparisons for discrete analysis quantities.

A small, closed vocabulary for comparing a discrete random quantity with a
finite threshold. Analysis modules use this selector for exact masses, lower
tails, and upper tails without accepting arbitrary predicates or constructing
a general event algebra.
-}
module Dtmc.Analysis.Event (
    DiscreteEvent (..),
    matches,
    includesInfiniteOutcome,
) where

import Numeric.Natural (
    Natural,
 )

{- | A comparison between a non-negative integer-valued quantity and a finite
threshold.

If the quantity can equal infinity, that atom belongs to 'GreaterThan' and
'AtLeast' and to none of the other events. Eventual hitting, eventual return,
and infinitely many visits remain separate queries because this type contains
only finite thresholds.
-}
data DiscreteEvent
    = EqualTo Natural -- ^ The quantity equals the threshold.
    | LessThan Natural -- ^ The quantity is strictly below the threshold.
    | AtMost Natural -- ^ The quantity does not exceed the threshold.
    | GreaterThan Natural -- ^ The quantity is strictly above the threshold.
    | AtLeast Natural -- ^ The quantity is at least the threshold.
    deriving (Eq, Ord, Show)

{- | Test whether a finite value satisfies a 'DiscreteEvent'. This is the
literal closed comparison semantics and performs no probability calculation.
For the separate infinite outcome, use 'includesInfiniteOutcome'.

Complexity: @O(1)@ time and @O(1)@ space.
-}
matches :: DiscreteEvent -> Natural -> Bool
matches event value =
    case event of
        EqualTo threshold -> value == threshold
        LessThan threshold -> value < threshold
        AtMost threshold -> value <= threshold
        GreaterThan threshold -> value > threshold
        AtLeast threshold -> value >= threshold

{- | Test whether the event contains the infinity atom of an extended-natural
quantity. Every finite threshold is below infinity, so precisely the two
upper-tail comparisons contain it.

Complexity: @O(1)@ time and @O(1)@ space.
-}
includesInfiniteOutcome :: DiscreteEvent -> Bool
includesInfiniteOutcome event =
    case event of
        EqualTo _ -> False
        LessThan _ -> False
        AtMost _ -> False
        GreaterThan _ -> True
        AtLeast _ -> True
