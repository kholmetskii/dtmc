{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}

{- |
Module      : Dtmc.State.Internal
Description : Checked conversions for finite-state indices.

Shared conversions between named finite states and the integer indices used by
dynamic graph and linear-algebra code. Keeping the reverse conversion here
prevents an arbitrary 'Int' from being passed directly to 'Data.Finite.finite',
which wraps out-of-range values modulo the state-space cardinality.
-}
module Dtmc.State.Internal (
    stateCardinalityInt,
    stateIndexInt,
    stateFromInt,
) where

import Data.Finite (
    finite,
    getFinite,
 )
import Data.Proxy (
    Proxy (Proxy),
 )
import Dtmc.State (
    Cardinality,
    FiniteState,
    stateAt,
    stateIndex,
 )
import GHC.TypeNats (
    natVal,
 )

{- | Return the number of states as a runtime 'Int'.

Complexity: @O(1)@ time and @O(1)@ space.
-}
stateCardinalityInt :: forall state. (FiniteState state) => Int
stateCardinalityInt =
    fromIntegral (natVal (Proxy @(Cardinality state)))

{- | Return the canonical zero-based integer index of a state.

Complexity: the cost of 'stateIndex' plus @O(1)@ time and @O(1)@ space.
-}
stateIndexInt :: (FiniteState state) => state -> Int
stateIndexInt = fromIntegral . getFinite . stateIndex

{- | Recover the state at a runtime integer index. Returns 'Nothing' rather
than wrapping a negative or out-of-range integer modulo the state count.

Complexity: the cost of 'stateAt' plus @O(1)@ time and @O(1)@ space.
-}
stateFromInt :: forall state. (FiniteState state) => Int -> Maybe state
stateFromInt index
    | index < 0 || index >= stateCardinalityInt @state = Nothing
    | otherwise = Just (stateAt (finite (fromIntegral index)))
