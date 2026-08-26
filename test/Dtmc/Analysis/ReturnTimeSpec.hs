module Dtmc.Analysis.ReturnTimeSpec (
    spec,
) where

import Dtmc.Analysis.TimeSpecSupport (
    returnTimeSpec,
 )
import Test.Hspec (
    Spec,
 )

spec :: Spec
spec = returnTimeSpec
