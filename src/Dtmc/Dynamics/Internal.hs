{- |
Module      : Dtmc.Dynamics.Internal
Description : Reserved internal companion for forward dynamics (currently empty).

"Dtmc.Dynamics" exposes no raw carrier type or unsafe primitive of its own -- it
operates directly on the carriers in "Dtmc.Distribution.Internal" and
"Dtmc.TransitionMatrix.Internal" -- so this companion module is intentionally
empty. It exists to keep the public/internal module pairing uniform across the
library and to give future unchecked dynamics helpers a home without a later
public-facing move.
-}
module Dtmc.Dynamics.Internal () where
