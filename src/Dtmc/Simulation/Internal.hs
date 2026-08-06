{- |
Module      : Dtmc.Simulation.Internal
Description : Reserved internal companion for simulation (currently empty).

"Dtmc.Simulation" exposes no raw carrier type or unsafe primitive of its own --
it samples directly from the carriers in "Dtmc.Distribution.Internal" and
"Dtmc.Simplex.Internal" -- so this companion module is intentionally empty. It
exists to keep the public/internal module pairing uniform across the library
and to give future unchecked sampling helpers a home without a later
public-facing move.
-}
module Dtmc.Simulation.Internal () where
