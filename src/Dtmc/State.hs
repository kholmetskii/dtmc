{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE DefaultSignatures #-}
{-# LANGUAGE EmptyCase #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE UndecidableSuperClasses #-}

{- |
Module      : Dtmc.State
Description : Canonical indexing for finite named state types.

'FiniteState' identifies a globally finite state type with a canonical total
indexing by @'Finite' ('Cardinality' state)@. The generic implementation
supports enumeration types whose constructors have no fields. Enable
@DeriveAnyClass@ and @DeriveGeneric@, then derive both 'Generic' and
'FiniteState':

@
data Weather = Dry | Wet | Storm
    deriving (Eq, Ord, Show, Generic, FiniteState)
@

After deriving 'Generic', declaring an empty @instance FiniteState Weather@
is a more explicit way to select the same generic defaults.

Constructor declaration order determines vector and matrix order. A stock
derived 'Ord' instance has the same order and is the intended companion;
handwritten 'FiniteState' and 'Ord' instances are trusted to preserve the
documented ordering and bijection laws.
-}
module Dtmc.State (
    type Cardinality,
    type GenericCardinality,
    FiniteState,
    finiteStates,
    stateIndex,
    stateAt,
) where

import Data.Finite (
    Finite,
    finite,
    finites,
    getFinite,
 )
import Data.Kind (
    Type,
 )
import Data.Proxy (
    Proxy (Proxy),
 )
import GHC.Generics (
    C,
    D,
    Generic (Rep, from, to),
    M1 (M1),
    U1 (U1),
    V1,
    type (:+:) (L1, R1),
 )
import GHC.TypeLits (
    ErrorMessage (Text),
    TypeError,
 )
import GHC.TypeNats (
    KnownNat,
    Nat,
    natVal,
    type (+),
 )

{- | The number of inhabitants of a finite state type. For 'Finite', this is
its existing type-level bound; for every other type, it is derived from its
'Generic' representation. Users cannot override this closed family
independently of that representation.

A generic constructor carrying any fields reduces to a custom 'TypeError'.
-}
type family Cardinality (state :: Type) :: Nat where
    Cardinality (Finite n) = n
    Cardinality state = GenericCardinality (Rep state)

{- | A finite state type with a canonical bijection to
@'Finite' ('Cardinality' state)@.

Instances must satisfy:

* @stateAt (stateIndex state) == state@;
* @stateIndex (stateAt index) == index@;
* @finiteStates == map stateAt finites@;
* @finiteStates@ is strictly ascending according to 'Ord'.

The generic defaults satisfy these laws for fieldless enumeration types with
a stock derived 'Ord' instance. Handwritten method implementations are trusted
to satisfy them, but their 'Cardinality' still comes from the supported
'Generic' representation. Empty state types are supported: their state list is
empty and 'stateAt' has an uninhabited 'Finite 0' domain.
-}
class (Ord state, KnownNat (Cardinality state)) => FiniteState state where
    {- | Return every state exactly once, in canonical index order.

    Complexity: implementation-dependent.
    -}
    finiteStates :: [state]

    {- | Convert a state to its total, statically bounded index.

    Complexity: implementation-dependent.
    -}
    stateIndex :: state -> Finite (Cardinality state)

    {- | Recover the state at a statically bounded index.

    Complexity: implementation-dependent.
    -}
    stateAt :: Finite (Cardinality state) -> state

    default finiteStates ::
        ( Generic state
        , GenericFiniteState (Rep state)
        ) =>
        [state]
    finiteStates = map to genericStates

    default stateIndex ::
        ( Generic state
        , GenericFiniteState (Rep state)
        ) =>
        state ->
        Finite (Cardinality state)
    stateIndex = finite . genericIndex . from

    default stateAt ::
        ( Generic state
        , GenericFiniteState (Rep state)
        ) =>
        Finite (Cardinality state) ->
        state
    stateAt = to . genericAt . fromIntegral . getFinite

instance (KnownNat n) => FiniteState (Finite n) where
    finiteStates = finites
    stateIndex = id
    stateAt = id

instance FiniteState ()

instance FiniteState Bool

instance FiniteState Ordering

{- | The type-level cardinality of a 'Generic' representation. This advanced
helper underlies 'Cardinality'; ordinary users should use 'FiniteState'
instead.
-}
type family GenericCardinality (representation :: Type -> Type) :: Nat where
    GenericCardinality (M1 D metadata representation) =
        GenericCardinality representation
    GenericCardinality (left :+: right) =
        GenericCardinality left + GenericCardinality right
    GenericCardinality (M1 C metadata U1) = 1
    GenericCardinality (M1 C metadata fields) =
        TypeError
            ( 'Text
                "FiniteState: constructors with fields are unsupported"
            )
    GenericCardinality V1 = 0

-- Generic machinery implementing the canonical state/index bijection for
-- fieldless enumeration representations.
class GenericFiniteState representation where
    genericStates :: [representation value]
    genericIndex :: representation value -> Integer
    genericAt :: Integer -> representation value

instance
    (GenericFiniteState representation) =>
    GenericFiniteState (M1 D metadata representation)
    where
    genericStates = map M1 genericStates
    genericIndex (M1 value) = genericIndex value
    genericAt = M1 . genericAt

instance
    ( GenericFiniteState left
    , GenericFiniteState right
    , KnownNat (GenericCardinality left)
    ) =>
    GenericFiniteState (left :+: right)
    where
    genericStates =
        map L1 (genericStates @left)
            ++ map R1 (genericStates @right)

    genericIndex (L1 value) = genericIndex value
    genericIndex (R1 value) = genericCardinality @left + genericIndex value

    genericAt index
        | index < genericCardinality @left = L1 (genericAt index)
        | otherwise =
            R1 (genericAt (index - genericCardinality @left))

genericCardinality ::
    forall representation.
    (KnownNat (GenericCardinality representation)) =>
    Integer
genericCardinality =
    fromIntegral (natVal (Proxy @(GenericCardinality representation)))

instance {-# OVERLAPPING #-} GenericFiniteState (M1 C metadata U1) where
    genericStates = [M1 U1]
    genericIndex (M1 U1) = 0
    genericAt _ = M1 U1

instance
    {-# OVERLAPPABLE #-}
    ( TypeError
        ( 'Text
            "FiniteState: constructors with fields are unsupported"
        )
    ) =>
    GenericFiniteState (M1 C metadata fields)
    where
    genericStates = unsupportedConstructorFields
    genericIndex _ = unsupportedConstructorFields
    genericAt _ = unsupportedConstructorFields

-- Required only to complete an instance made unusable by its 'TypeError'.
unsupportedConstructorFields :: value
unsupportedConstructorFields =
    error "Dtmc.State: constructors with fields are unsupported"

instance GenericFiniteState V1 where
    genericStates = []
    genericIndex value = case value of {}
    genericAt _ =
        error "Dtmc.State.stateAt: unreachable Finite 0 index"
