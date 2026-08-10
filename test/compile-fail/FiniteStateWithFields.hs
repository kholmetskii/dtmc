{-# LANGUAGE DeriveGeneric #-}

-- This fixture must fail to compile with:
--
-- FiniteState: constructors with fields are unsupported
module FiniteStateWithFields where

import Dtmc.State (
    FiniteState,
 )
import GHC.Generics (
    Generic,
 )

data InvalidState = InvalidState Int
    deriving (Generic)

instance FiniteState InvalidState
