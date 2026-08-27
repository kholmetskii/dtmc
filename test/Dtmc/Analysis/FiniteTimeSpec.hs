{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.FiniteTimeSpec (
    spec,
) where

import Data.Finite (
    Finite,
    finites,
    getFinite,
 )
import Data.List.NonEmpty (
    NonEmpty ((:|)),
 )
import Dtmc.Analysis.FiniteTime (
    ConditionalProbabilityError (..),
    Observation (..),
    conditionalProbability,
    jointProbability,
    pathProbability,
    probabilityAtTime,
    transitionProbability,
    transitionProbabilityN,
 )
import Dtmc.Distribution (
    probabilityAt,
 )
import Dtmc.Distribution.Map qualified as DistributionMap
import Dtmc.Distribution.Vector (
    DistributionVector,
    mkDistributionVector,
 )
import Dtmc.Dynamics (
    evolveVector,
    evolveVectorN,
 )
import Dtmc.State qualified
import Dtmc.TestSupport (
    approxEq,
    genSimplexPoint,
    genTransitionMatrix,
    testTolerance,
 )
import Dtmc.Transition.Kernel qualified as Kernel
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    matrixPower,
    mkTransitionMatrix,
    rowAt,
    unTransitionMatrix,
 )
import GHC.Generics (
    Generic,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Numeric.Natural (
    Natural,
 )
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
    choose,
    conjoin,
    counterexample,
    forAll,
    property,
    (===),
 )

-- A three-state chain with several impossible one-step transitions.
chain :: TransitionMatrix (Finite 3)
chain =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 0.5
                , 0.5
                , 0.0
                , 0.0
                , 0.2
                , 0.8
                , 1.0
                , 0.0
                , 0.0
                ] ::
                S.Sq 3
            )

initial :: DistributionVector (Finite 3)
initial =
    either (error . show) id $
        mkDistributionVector (S.vector [0.6, 0.3, 0.1] :: S.R 3)

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

asTransitionKernel ::
    (Dtmc.State.FiniteState state) =>
    TransitionMatrix state ->
    Kernel.TransitionKernel state
asTransitionKernel matrix =
    Kernel.transitionKernel $ \source ->
        checked $
            DistributionMap.mkDistributionMap
                [ (destination, transitionProbability matrix source destination)
                | destination <- Dtmc.State.finiteStates
                ]

kernelChain :: Kernel.TransitionKernel (Finite 3)
kernelChain = asTransitionKernel chain

mapInitial :: DistributionMap.DistributionMap (Finite 3)
mapInitial =
    checked $
        DistributionMap.mkDistributionMap
            [ (state, probabilityAt initial state)
            | state <- Dtmc.State.finiteStates
            ]

simpleRandomWalk :: Kernel.TransitionKernel Integer
simpleRandomWalk =
    Kernel.transitionKernel $ \state ->
        checked
            (DistributionMap.mkDistributionMap [(state - 1, 0.5), (state + 1, 0.5)])

closeTo :: Double -> Double -> Bool
closeTo = approxEq testTolerance

{- | Hold for a @Right@ whose 'Double' is within 'testTolerance' of the
expected value; fail for any @Left@ or out-of-tolerance value.
-}
rightCloseTo :: Double -> Either ConditionalProbabilityError Double -> Bool
rightCloseTo expected (Right actual) = approxEq testTolerance actual expected
rightCloseTo _ (Left _) = False

rightResultsClose :: Either error Double -> Either error Double -> Bool
rightResultsClose (Right left) (Right right) = closeTo left right
rightResultsClose (Left _) (Left _) = True
rightResultsClose _ _ = False

data NamedPhase = PhaseA | PhaseB | PhaseC
    deriving (Eq, Ord, Show, Generic)

instance Dtmc.State.FiniteState NamedPhase

namedCycle :: TransitionMatrix NamedPhase
namedCycle =
    either (error . show) id $
        mkTransitionMatrix @NamedPhase
            (S.matrix [0, 1, 0, 0, 0, 1, 1, 0, 0] :: S.Sq 3)

twoState :: TransitionMatrix (Finite 2)
twoState =
    either (error . show) id $
        mkTransitionMatrix
            (S.matrix [0.9, 0.1, 0.4, 0.6] :: S.Sq 2)

twoStateSquared :: TransitionMatrix (Finite 2)
twoStateSquared =
    either (error . show) id $
        mkTransitionMatrix
            (S.matrix [0.85, 0.15, 0.6, 0.4] :: S.Sq 2)

closedFormTransition :: TransitionMatrix (Finite 3)
closedFormTransition =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 0.1
                , 0.5
                , 0.4
                , 0.1
                , 0.8
                , 0.1
                , 0.0
                , 0.5
                , 0.5
                ] ::
                S.Sq 3
            )

closedFormProbability :: Int -> Double
closedFormProbability n =
    5 / 63 + 5 / 18 * (0.1 ^ n) - 5 / 14 * (0.3 ^ n)

{- | Five-state transition matrix over states @[A, B, C, D, E]@ used by the
probability examples.
-}
observationMatrix :: TransitionMatrix (Finite 5)
observationMatrix =
    either (error . show) id $
        mkTransitionMatrix
            ( S.matrix
                [ 0
                , 0
                , 0
                , 1
                , 0
                , 1 / 3
                , 0
                , 0
                , 0
                , 2 / 3
                , 0
                , 0
                , 0
                , 0
                , 1
                , 0
                , 0
                , 1 / 3
                , 2 / 3
                , 0
                , 1 / 4
                , 1 / 4
                , 0
                , 0
                , 1 / 2
                ] ::
                S.Sq 5
            )

-- | Initial law @lambda = [1/4, 1/2, 0, 1/4, 0]@ for the probability examples.
observationInitial :: DistributionVector (Finite 5)
observationInitial =
    either (error . show) id $
        mkDistributionVector (S.vector [1 / 4, 1 / 2, 0, 1 / 4, 0] :: S.R 5)

spec :: Spec
spec = do
    describe "Observation" $ do
        it "is polymorphic in the state type" $
            (At 2 "rain" :: Observation String) `shouldBe` At 2 "rain"

    describe "transitionProbability" $ do
        prop "agrees with rowAt then probabilityAt" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                case mkTransitionMatrix @(Finite 3) matrix of
                    Right p ->
                        conjoin
                            [ transitionProbability p i j
                                === probabilityAt (rowAt p i) j
                            | i <- finites
                            , j <- finites
                            ]
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

        it "uses named state constructors" $
            transitionProbability namedCycle PhaseB PhaseC
                `shouldBe` 1

    describe "transitionProbabilityN" $ do
        it "is the Kronecker delta at exponent zero" $
            let ijs =
                    [(i, j) | i <- finites, j <- finites] ::
                        [(Finite 2, Finite 2)]
             in map (uncurry (transitionProbabilityN 0 twoState)) ijs
                    `shouldBe` map (\(i, j) -> if i == j then 1 else 0) ijs

        prop "agrees with transitionProbability at exponent one" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                case mkTransitionMatrix @(Finite 3) matrix of
                    Right p ->
                        conjoin
                            [ property $
                                approxEq
                                    testTolerance
                                    (transitionProbabilityN 1 p i j)
                                    (transitionProbability p i j)
                            | i <- finites
                            , j <- finites
                            ]
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

        it "matches a hand-computed square at exponent two" $
            sequence_
                [ transitionProbabilityN 2 twoState i j
                    `shouldSatisfy` closeTo (probabilityAt (rowAt twoStateSquared i) j)
                | i <- finites :: [Finite 2]
                , j <- finites :: [Finite 2]
                ]

        prop "agrees with the corresponding matrixPower entry" $
            forAll (genTransitionMatrix @3) $ \matrix ->
                case mkTransitionMatrix @(Finite 3) matrix of
                    Right p ->
                        let power = S.extract (unTransitionMatrix (matrixPower 4 p))
                         in conjoin
                                [ property $
                                    approxEq
                                        testTolerance
                                        (transitionProbabilityN 4 p i j)
                                        ( power
                                            `LA.atIndex` ( fromIntegral (getFinite i)
                                                         , fromIntegral (getFinite j)
                                                         )
                                        )
                                | i <- finites
                                , j <- finites
                                ]
                    Left err ->
                        counterexample
                            ("generated matrix was rejected: " <> show err)
                            False

        it "preserves named state types" $
            transitionProbabilityN 2 namedCycle PhaseA PhaseC
                `shouldBe` 1

    describe "transitionProbabilityN hand-computed regressions" $ do
        it "gives P^3(E, D) = 3/8 for the five-state chain" $
            transitionProbabilityN 3 observationMatrix 4 3
                `shouldSatisfy` closeTo (3 / 8)

        it "matches the three-state P^n(2, 0) closed form" $
            mapM_
                ( \n ->
                    transitionProbabilityN n closedFormTransition 2 0
                        `shouldSatisfy` closeTo
                            (closedFormProbability (fromIntegral n))
                )
                ([0, 1, 2, 3, 5, 10, 20] :: [Natural])

    describe "probabilityAtTime" $ do
        it "returns the initial probability at time zero" $
            conjoin
                [ probabilityAtTime 0 initial chain state
                    === probabilityAt initial state
                | state <- finites
                ]

        prop "agrees with probabilityAt of evolveVectorN"
            $ forAll
                ( (,,)
                    <$> choose (0, 6 :: Int)
                    <*> genSimplexPoint 3
                    <*> genTransitionMatrix @3
                )
            $ \(k, entries, matrix) ->
                case ( mkDistributionVector @(Finite 3) (S.vector entries :: S.R 3)
                     , mkTransitionMatrix matrix
                     ) of
                    (Right mu, Right p) ->
                        conjoin
                            [ property $
                                approxEq
                                    testTolerance
                                    (probabilityAtTime (fromIntegral k) mu p state)
                                    (probabilityAt (evolveVectorN (fromIntegral k) mu p) state)
                            | state <- finites
                            ]
                    result ->
                        counterexample
                            ("generated input was rejected: " <> show result)
                            False

        prop "agrees with repeated evolveVector for small exponents"
            $ forAll
                ( (,,)
                    <$> choose (0, 6 :: Int)
                    <*> genSimplexPoint 3
                    <*> genTransitionMatrix @3
                )
            $ \(k, entries, matrix) ->
                case ( mkDistributionVector @(Finite 3) (S.vector entries :: S.R 3)
                     , mkTransitionMatrix matrix
                     ) of
                    (Right mu, Right p) ->
                        let iterated = iterate (`evolveVector` p) mu !! k
                         in conjoin
                                [ property $
                                    approxEq
                                        testTolerance
                                        (probabilityAtTime (fromIntegral k) mu p state)
                                        (probabilityAt iterated state)
                                | state <- finites
                                ]
                    result ->
                        counterexample
                            ("generated input was rejected: " <> show result)
                            False

    describe "Transition realization independence" $ do
        it "computes transition probabilities on an infinite state type" $ do
            transitionProbabilityN 2 simpleRandomWalk 0 0
                `shouldSatisfy` closeTo 0.5
            transitionProbabilityN 3 simpleRandomWalk 0 0
                `shouldBe` 0

        prop "gives matrices and equivalent kernels the same transition powers" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left problem -> counterexample (show problem) False
                    Right matrix ->
                        let kernel = asTransitionKernel matrix
                         in conjoin
                                [ property $
                                    closeTo
                                        (transitionProbabilityN time matrix source destination)
                                        (transitionProbabilityN time kernel source destination)
                                | source <- finites :: [Finite 3]
                                , destination <- finites :: [Finite 3]
                                , time <- [0 .. 4]
                                ]

        it "matches finite trajectory and observation queries" $ do
            pathProbability mapInitial kernelChain (0 :| [1, 2])
                `shouldSatisfy` closeTo
                    (pathProbability initial chain (0 :| [1, 2]))
            jointProbability
                mapInitial
                kernelChain
                [At 3 2, At 0 0, At 1 1]
                `shouldSatisfy` closeTo
                    (jointProbability initial chain [At 3 2, At 0 0, At 1 1])

        it "matches finite conditional probability queries" $
            rightResultsClose
                (conditionalProbability mapInitial kernelChain [At 2 2] [At 0 0])
                (conditionalProbability initial chain [At 2 2] [At 0 0])
                `shouldBe` True

    describe "pathProbability" $ do
        it "returns the initial probability for a one-state path" $
            approxEq
                testTolerance
                (pathProbability initial chain (0 :| []))
                (probabilityAt initial 0)
                `shouldBe` True

        it "is lambda_i * P(i, j) for a two-state path" $
            approxEq
                testTolerance
                (pathProbability initial chain (0 :| [1]))
                (0.6 * 0.5)
                `shouldBe` True

        it "is the product of initial and transition probabilities" $
            approxEq
                testTolerance
                (pathProbability initial chain (0 :| [1, 2]))
                (0.6 * 0.5 * 0.8)
                `shouldBe` True

        it "is zero for a path with an impossible transition" $ do
            approxEq
                testTolerance
                (pathProbability initial chain (0 :| [2]))
                0
                `shouldBe` True
            approxEq
                testTolerance
                (pathProbability initial chain (0 :| [1, 0]))
                0
                `shouldBe` True

        prop "a one-state path equals the initial probability" $
            forAll ((,) <$> genSimplexPoint 3 <*> genTransitionMatrix @3) $
                \(entries, matrix) ->
                    case ( mkDistributionVector @(Finite 3) (S.vector entries :: S.R 3)
                         , mkTransitionMatrix matrix
                         ) of
                        (Right mu, Right p) ->
                            conjoin
                                [ pathProbability mu p (i :| [])
                                    === probabilityAt mu i
                                | i <- [0, 1, 2]
                                ]
                        result ->
                            counterexample
                                ("generated input was rejected: " <> show result)
                                False

        prop "a two-state path equals lambda_i * P(i, j)" $
            forAll ((,) <$> genSimplexPoint 3 <*> genTransitionMatrix @3) $
                \(entries, matrix) ->
                    case ( mkDistributionVector @(Finite 3) (S.vector entries :: S.R 3)
                         , mkTransitionMatrix matrix
                         ) of
                        (Right mu, Right p) ->
                            conjoin
                                [ property $
                                    approxEq
                                        testTolerance
                                        (pathProbability mu p (i :| [j]))
                                        ( probabilityAt mu i
                                            * transitionProbability p i j
                                        )
                                | i <- [0, 1, 2]
                                , j <- [0, 1, 2]
                                ]
                        result ->
                            counterexample
                                ("generated input was rejected: " <> show result)
                                False

    describe "jointProbability" $ do
        it "returns exactly one for no observations" $
            jointProbability initial chain [] `shouldBe` 1

        it "equals probabilityAtTime for a single observation" $
            approxEq
                testTolerance
                (jointProbability initial chain [At 1 1])
                (probabilityAtTime 1 initial chain 1)
                `shouldBe` True

        it "is unchanged by observation order" $
            approxEq
                testTolerance
                (jointProbability initial chain [At 0 0, At 1 1])
                (jointProbability initial chain [At 1 1, At 0 0])
                `shouldBe` True

        it "is unchanged by duplicate observations" $
            approxEq
                testTolerance
                (jointProbability initial chain [At 0 0, At 0 0, At 1 1])
                (jointProbability initial chain [At 0 0, At 1 1])
                `shouldBe` True

        it "is exactly zero for conflicting states at one time" $
            jointProbability initial chain [At 0 0, At 0 1] `shouldBe` 0

        it "agrees with pathProbability over times 0, 1, 2" $
            approxEq
                testTolerance
                (jointProbability initial chain [At 0 0, At 1 1, At 2 2])
                (pathProbability initial chain (0 :| [1, 2]))
                `shouldBe` True

        it "is exactly zero through an impossible transition" $
            jointProbability initial chain [At 0 0, At 1 2] `shouldBe` 0

        it "matches a hand-computed multi-gap example" $
            approxEq
                testTolerance
                ( jointProbability
                    observationInitial
                    observationMatrix
                    [At 2 2, At 3 4, At 6 3]
                )
                (5 / 96)
                `shouldBe` True

        prop "a single observation equals probabilityAtTime" $
            forAll ((,) <$> genSimplexPoint 3 <*> genTransitionMatrix @3) $
                \(entries, matrix) ->
                    case ( mkDistributionVector @(Finite 3) (S.vector entries :: S.R 3)
                         , mkTransitionMatrix matrix
                         ) of
                        (Right mu, Right p) ->
                            conjoin
                                [ property $
                                    approxEq
                                        testTolerance
                                        (jointProbability mu p [At t i])
                                        (probabilityAtTime t mu p i)
                                | t <- [0, 1, 2]
                                , i <- [0, 1, 2]
                                ]
                        result ->
                            counterexample
                                ("generated input was rejected: " <> show result)
                                False

        prop "is invariant under observation order" $
            forAll ((,) <$> genSimplexPoint 3 <*> genTransitionMatrix @3) $
                \(entries, matrix) ->
                    case ( mkDistributionVector @(Finite 3) (S.vector entries :: S.R 3)
                         , mkTransitionMatrix matrix
                         ) of
                        (Right mu, Right p) ->
                            property $
                                approxEq
                                    testTolerance
                                    (jointProbability mu p [At 1 1, At 3 2])
                                    (jointProbability mu p [At 3 2, At 1 1])
                        result ->
                            counterexample
                                ("generated input was rejected: " <> show result)
                                False

    describe "conditionalProbability" $ do
        it "returns the event probability for an empty condition" $
            conditionalProbability initial chain [At 1 1] []
                `shouldSatisfy` rightCloseTo
                    (jointProbability initial chain [At 1 1])

        it "returns one for an empty event and a positive condition" $
            conditionalProbability initial chain [] [At 0 0]
                `shouldSatisfy` rightCloseTo 1

        it "returns one when conditioning an observation on itself" $
            conditionalProbability initial chain [At 1 1] [At 1 1]
                `shouldSatisfy` rightCloseTo 1

        it "ignores observations shared by event and condition" $
            conditionalProbability initial chain [At 1 1] [At 1 1, At 2 2]
                `shouldSatisfy` rightCloseTo 1

        it "returns zero for a conflict against a possible condition" $
            conditionalProbability initial chain [At 1 0] [At 1 1]
                `shouldSatisfy` rightCloseTo 0

        it "reports a zero-probability condition" $
            conditionalProbability initial chain [At 0 0] [At 0 0, At 1 2]
                `shouldBe` Left ZeroProbabilityCondition

        it "reports a contradictory condition" $
            conditionalProbability initial chain [At 0 0] [At 1 1, At 1 2]
                `shouldBe` Left ZeroProbabilityCondition

        it "is unaffected by event and condition ordering" $ do
            conditionalProbability initial chain [At 2 2, At 1 1] [At 0 0]
                `shouldSatisfy` rightCloseTo 0.4
            conditionalProbability initial chain [At 1 1, At 2 2] [At 0 0]
                `shouldSatisfy` rightCloseTo 0.4

    describe "conditionalProbability hand-computed regressions" $ do
        it "gives P(X10=D, X11=D | X3=A, X7=E) = 1/4" $
            conditionalProbability
                observationInitial
                observationMatrix
                [At 10 3, At 11 3]
                [At 3 0, At 7 4]
                `shouldSatisfy` rightCloseTo (1 / 4)

        it "accepts an out-of-order event and gives 15/92" $
            conditionalProbability
                observationInitial
                observationMatrix
                [At 6 3, At 2 2]
                [At 3 4]
                `shouldSatisfy` rightCloseTo (15 / 92)

        it "gives P(X2=C) = 5/36" $
            approxEq
                testTolerance
                (probabilityAtTime 2 observationInitial observationMatrix 2)
                (5 / 36)
                `shouldBe` True

        it "gives P(X3=E) = 23/72" $
            approxEq
                testTolerance
                (probabilityAtTime 3 observationInitial observationMatrix 4)
                (23 / 72)
                `shouldBe` True

        it "gives P(X2=C, X3=E, X6=D) = 5/96" $
            approxEq
                testTolerance
                ( jointProbability
                    observationInitial
                    observationMatrix
                    [At 2 2, At 3 4, At 6 3]
                )
                (5 / 96)
                `shouldBe` True
