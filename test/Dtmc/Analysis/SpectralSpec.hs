{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.SpectralSpec (
    spec,
) where

import Data.Complex (
    magnitude,
    realPart,
 )
import Data.Finite (
    Finite,
 )
import Dtmc.Analysis.Spectral (
    secondLargestModulus,
    spectralGap,
    spectrum,
 )
import Dtmc.TestSupport (
    approxEq,
    genTransitionMatrix,
 )
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    TransitionMatrixError,
    matrixPower,
    mkTransitionMatrix,
    unTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    it,
    shouldBe,
    shouldSatisfy,
 )
import Test.Hspec.QuickCheck (
    prop,
 )
import Test.QuickCheck (
    conjoin,
    counterexample,
    forAll,
    property,
 )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

-- Section 2.4 exercise: eigenvalues 1, 1/2 and 1/5, all real and distinct.
distinctSpectrum :: TransitionMatrix (Finite 3)
distinctSpectrum =
    checked
        ( mkTransitionMatrix
            (S.matrix [0.3, 0.5, 0.2, 0.1, 0.6, 0.3, 0.1, 0.1, 0.8] :: S.Sq 3)
        )

-- Section 4.2 example: a conjugate pair (-5 +- i sqrt 5) / 10 of modulus
-- sqrt (3 / 10).
complexSpectrum :: TransitionMatrix (Finite 3)
complexSpectrum =
    checked
        ( mkTransitionMatrix
            (S.matrix [0, 1, 0, 0.4, 0, 0.6, 0.5, 0.5, 0] :: S.Sq 3)
        )

-- Irreducible and aperiodic: eigenvalues 1 and 1/2.
twoState :: TransitionMatrix (Finite 2)
twoState =
    checked
        ( mkTransitionMatrix
            (S.matrix [0.9, 0.1, 0.4, 0.6] :: S.Sq 2)
        )

-- Period 3: the three cube roots of unity sit on the unit circle.
threeCycle :: TransitionMatrix (Finite 3)
threeCycle =
    checked
        ( mkTransitionMatrix
            (S.matrix [0, 1, 0, 0, 0, 1, 1, 0, 0] :: S.Sq 3)
        )

singleton :: TransitionMatrix (Finite 1)
singleton =
    checked (mkTransitionMatrix (S.matrix [1] :: S.Sq 1))

closeList :: [Double] -> [Double] -> Bool
closeList expected actual =
    length expected == length actual
        && and (zipWith (approxEq 1e-8) expected actual)

-- The diagonalisation formula of section 2.4, used only as an oracle:
-- P^n = Q D^n Q^-1 for a diagonalisable P.
diagonalisedPower :: Int -> TransitionMatrix (Finite 3) -> [[Double]]
diagonalisedPower steps p =
    LA.toLists (LA.cmap realPart (vectors LA.<> scaled LA.<> LA.inv vectors))
  where
    (values, vectors) = LA.eig (S.extract (unTransitionMatrix p))
    scaled = LA.diag (LA.cmap (^ steps) values)

squaredPower :: Int -> TransitionMatrix (Finite 3) -> [[Double]]
squaredPower steps p =
    LA.toLists (S.extract (unTransitionMatrix (matrixPower (fromIntegral steps) p)))

spec :: Spec
spec = do
    describe "spectrum" $ do
        it "reproduces the real eigenvalues of the notes" $
            map magnitude (spectrum distinctSpectrum)
                `shouldSatisfy` closeList [1, 0.5, 0.2]

        it "reproduces the complex pair of the notes" $
            map magnitude (spectrum complexSpectrum)
                `shouldSatisfy` closeList [1, sqrt 0.3, sqrt 0.3]

        it "puts every cube root of unity on the unit circle" $
            map magnitude (spectrum threeCycle)
                `shouldSatisfy` closeList [1, 1, 1]

        prop "starts at one and never exceeds it" $
            forAll (genTransitionMatrix @3) $ \m ->
                case mkTransitionMatrix m ::
                        Either TransitionMatrixError (TransitionMatrix (Finite 3)) of
                    Left err -> counterexample (show err) False
                    Right p ->
                        let values = map magnitude (spectrum p)
                         in conjoin
                                [ counterexample
                                    ("no Perron root: " <> show values)
                                    (property (any (approxEq 1e-8 1) values))
                                , counterexample
                                    ("modulus above one: " <> show values)
                                    (property (all (<= 1 + 1e-8) values))
                                ]

    describe "secondLargestModulus" $ do
        it "is the convergence rate of an ergodic chain" $
            secondLargestModulus twoState `shouldSatisfy` maybe False (approxEq 1e-8 0.5)

        it "is one for a periodic chain" $
            secondLargestModulus threeCycle `shouldSatisfy` maybe False (approxEq 1e-8 1)

        it "is absent for a one-state chain" $
            secondLargestModulus singleton `shouldBe` Nothing

    describe "spectralGap" $ do
        it "complements the second largest modulus" $
            spectralGap twoState `shouldSatisfy` maybe False (approxEq 1e-8 0.5)

        it "vanishes for a periodic chain" $
            spectralGap threeCycle `shouldSatisfy` maybe False (approxEq 1e-8 0)

    describe "diagonalisation" $
        it "reproduces matrixPower on a chain with distinct eigenvalues" $
            -- The proposition of section 2.4, checked against repeated
            -- squaring. Not exported: matrixPower computes it directly.
            concat (diagonalisedPower 5 distinctSpectrum)
                `shouldSatisfy` closeList (concat (squaredPower 5 distinctSpectrum))
