{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeApplications #-}

module Dtmc.Analysis.CanonicalDifferentialSpec (
    spec,
) where

import Data.Finite (
    Finite,
    finites,
 )
import Data.List.NonEmpty (
    NonEmpty ((:|)),
 )
import Data.Maybe (
    fromMaybe,
 )
import Dtmc.Analysis.Event (
    DiscreteEvent (..),
 )
import Dtmc.Analysis.Expectation (
    Expectation (..),
 )
import Dtmc.Analysis.FiniteTime qualified as FT
import Dtmc.Analysis.HittingTime qualified as Hit
import Dtmc.Analysis.ProbabilityOracle qualified as Oracle
import Dtmc.Analysis.ReturnTime qualified as Return
import Dtmc.Analysis.VisitCount qualified as Visit
import Dtmc.Distribution (
    probabilityAt,
 )
import Dtmc.Distribution.Vector (
    DistributionVector,
    mkDistributionVector,
 )
import Dtmc.TestSupport
import Dtmc.Transition.Matrix (
    TransitionMatrix,
    mkTransitionMatrix,
 )
import Numeric.LinearAlgebra qualified as LA
import Numeric.LinearAlgebra.Static qualified as S
import Test.Hspec (
    Spec,
    describe,
    it,
 )
import Test.Hspec.QuickCheck (
    prop,
 )
import Test.QuickCheck (
    counterexample,
    forAll,
    property,
 )

initialWeights :: [(Finite 3, Double)]
initialWeights = zip finites [0.2, 0.3, 0.5]

initialDistribution :: DistributionVector (Finite 3)
initialDistribution =
    checked (mkDistributionVector (S.vector [0.2, 0.3, 0.5] :: S.R 3))

terminalChain :: TransitionMatrix (Finite 3)
terminalChain =
    checked
        ( mkTransitionMatrix
            ( S.matrix
                [ 0
                , 0.5
                , 0.5
                , 0
                , 0
                , 1
                , 0
                , 0
                , 1
                ] ::
                S.Sq 3
            )
        )

checked :: (Show error) => Either error value -> value
checked = either (error . show) id

entries :: S.R 3 -> [Double]
entries = LA.toList . S.extract

known :: Maybe Double -> Double
known = fromMaybe (error "oracle horizon does not determine this event")

close :: Double -> Double -> Bool
close = approxEq testTolerance

rightClose :: Double -> Either error Double -> Bool
rightClose expected = either (const False) (close expected)

expectationClose :: Maybe Double -> Expectation -> Bool
expectationClose Nothing InfiniteExpectation = True
expectationClose (Just expected) (FiniteExpectation actual) = close expected actual
expectationClose _ _ = False

finiteAndBoundedChecks :: TransitionMatrix (Finite 3) -> Bool
finiteAndBoundedChecks matrix =
    and
        [ transitionChecks
        , trajectoryChecks
        , observationChecks
        , hittingChecks
        , returnChecks
        , visitChecks
        ]
  where
    target state = state == (2 :: Finite 3)
    transitionChecks =
        and
            [ close
                (FT.stepProbability matrix source destination)
                (Oracle.transitionWeight matrix source destination)
            | source <- finites
            , destination <- finites
            ]
            && and
                [ close
                    (FT.nStepProbability time matrix source destination)
                    (Oracle.stateProbability time [(source, 1)] matrix destination)
                | time <- [0 .. 4]
                , source <- finites
                , destination <- finites
                ]
            && and
                [ close
                    (FT.stateProbability time initialDistribution matrix destination)
                    (Oracle.stateProbability time initialWeights matrix destination)
                | time <- [0 .. 4]
                , destination <- finites
                ]
    trajectoryChecks =
        close
            (FT.pathProbability initialDistribution matrix (0 :| [1, 2]))
            (Oracle.trajectoryProbability initialWeights matrix [0, 1, 2])
    observations = [(1, 1), (3, 2)]
    oracleJoint =
        Oracle.observationProbability 3 initialWeights matrix observations
    observationChecks =
        close
            (FT.observationProbability initialDistribution matrix [FT.At 1 1, FT.At 3 2])
            oracleJoint
            && conditionalChecks
    conditionalChecks =
        let denominator =
                Oracle.observationProbability 1 initialWeights matrix [(1, 1)]
            numerator = oracleJoint
            actual =
                FT.conditionalObservationProbability
                    initialDistribution
                    matrix
                    [FT.At 3 2]
                    [FT.At 1 1]
         in if denominator == 0
                then actual == Left FT.ZeroProbabilityCondition
                else either (const False) (close (numerator / denominator)) actual
    hittingChecks =
        and
            [ let law = Oracle.hittingLaw 4 matrix target source
                  exact = known (Oracle.lawProbability (EqualTo time) law)
                  dense = entries (hitProbabilityByState (EqualTo time) matrix [2])
               in close (Hit.probabilityGivenInitialState (EqualTo time) matrix target source) exact
                    && close (dense !! fromIntegral source) exact
            | source <- finites
            , time <- [0 .. 4]
            ]
            && and
                [ let law = Oracle.hittingLaw 4 matrix target source
                      bounded = known (Oracle.lawProbability (LessThan bound) law)
                      dense = entries (hitProbabilityByState (LessThan bound) matrix [2])
                   in close
                        (Hit.probabilityGivenInitialState (LessThan bound) matrix target source)
                        bounded
                        && close (dense !! fromIntegral source) bounded
                | source <- finites
                , bound <- [0 .. 5]
                ]
    returnChecks =
        and
            [ let law = Oracle.returnLaw 4 matrix source
                  exact = known (Oracle.lawProbability (EqualTo time) law)
                  dense = entries (returnProbabilityByState (EqualTo time) matrix)
               in close (Return.probabilityGivenInitialState (EqualTo time) matrix source) exact
                    && close (dense !! fromIntegral source) exact
            | source <- finites
            , time <- [0 .. 4]
            ]
            && and
                [ let law = Oracle.returnLaw 4 matrix source
                      bounded = known (Oracle.lawProbability (LessThan bound) law)
                      dense = entries (returnProbabilityByState (LessThan bound) matrix)
                   in close
                        (Return.probabilityGivenInitialState (LessThan bound) matrix source)
                        bounded
                        && close (dense !! fromIntegral source) bounded
                | source <- finites
                , bound <- [0 .. 5]
                ]
    visitChecks =
        and
            [ let law =
                    Oracle.visitLawBefore
                        bound
                        initialWeights
                        matrix
                        target
                  distribution =
                    Visit.boundedLaw
                        bound
                        initialDistribution
                        matrix
                        target
                  expected = known (Oracle.lawFiniteExpectation law)
               in and
                    [ close
                        (probabilityAt distribution count)
                        (known (Oracle.lawProbability (EqualTo count) law))
                        && close
                            ( Visit.boundedProbability
                                bound
                                (EqualTo count)
                                initialDistribution
                                matrix
                                target
                            )
                            (known (Oracle.lawProbability (EqualTo count) law))
                    | count <- [0 .. bound]
                    ]
                    && close
                        ( Visit.boundedExpectation
                            bound
                            initialDistribution
                            matrix
                            target
                        )
                        expected
            | bound <- [0 .. 4]
            ]

terminalChecks :: Bool
terminalChecks =
    and
        [ hittingEventualChecks
        , hittingRaceChecks
        , hittingExpectationChecks
        , returnEventualChecks
        , returnExpectationChecks
        , totalVisitChecks
        ]
  where
    states = finites :: [Finite 3]
    target state = state == (1 :: Finite 3)
    competing state = state == (2 :: Finite 3)
    hitLaws = [Oracle.hittingLaw 1 terminalChain target state | state <- states]
    returnLaws = [Oracle.returnLaw 1 terminalChain state | state <- states]
    eventual law = 1 - Oracle.lawUnresolvedMass law
    hitValues = map eventual hitLaws
    returnValues = map eventual returnLaws
    hittingEventualChecks =
        case hitEventualProbabilityByState terminalChain [1] of
            Left _ -> False
            Right dense ->
                and (zipWith close (entries dense) hitValues)
                    && and
                        [ rightClose expected (Hit.eventualProbabilityGivenInitialState terminalChain [1] state)
                        | (state, expected) <- zip states hitValues
                        ]
    raceValues =
        [ Oracle.raceProbabilityWithin
            1
            terminalChain
            target
            competing
            state
        | state <- states
        ]
    hittingRaceChecks =
        case hitRaceProbabilityByState terminalChain [1] [2] of
            Left _ -> False
            Right dense ->
                and (zipWith close (entries dense) raceValues)
                    && and
                        [ rightClose
                            expected
                            (Hit.raceProbabilityGivenInitialState terminalChain [1] [2] state)
                        | (state, expected) <- zip states raceValues
                        ]
    hitExpectations = map Oracle.lawFiniteExpectation hitLaws
    hittingExpectationChecks =
        case hitExpectationByState terminalChain [1] of
            Left _ -> False
            Right actual ->
                and (zipWith expectationClose hitExpectations actual)
                    && and
                        [ either
                            (const False)
                            (expectationClose expected)
                            (Hit.expectationGivenInitialState terminalChain [1] state)
                        | (state, expected) <- zip states hitExpectations
                        ]
    returnEventualChecks =
        case returnEventualProbabilityByState terminalChain of
            Left _ -> False
            Right dense ->
                and (zipWith close (entries dense) returnValues)
                    && and
                        [ rightClose expected (Return.eventualProbabilityGivenInitialState terminalChain state)
                        | (state, expected) <- zip states returnValues
                        ]
    returnExpectations = map Oracle.lawFiniteExpectation returnLaws
    returnExpectationChecks =
        case returnExpectationByState terminalChain of
            Left _ -> False
            Right actual ->
                and (zipWith expectationClose returnExpectations actual)
                    && and
                        [ either
                            (const False)
                            (expectationClose expected)
                            (Return.expectationGivenInitialState terminalChain state)
                        | (state, expected) <- zip states returnExpectations
                        ]
    visitLaws =
        [ Oracle.visitLawBefore 2 [(state, 1)] terminalChain target
        | state <- states
        ]
    visitExpectations = map Oracle.lawFiniteExpectation visitLaws
    totalVisitChecks =
        and
            [ case visitTotalProbabilityByState (EqualTo count) terminalChain 1 of
                Left _ -> False
                Right dense ->
                    and
                        [ let expected = known (Oracle.lawProbability (EqualTo count) law)
                           in close (entries dense !! fromIntegral state) expected
                                && rightClose
                                    expected
                                    (Visit.totalProbabilityGivenInitialState (EqualTo count) terminalChain 1 state)
                        | (state, law) <- zip states visitLaws
                        ]
            | count <- [0 .. 2]
            ]
            && case visitInfiniteProbabilityByState terminalChain 1 of
                Left _ -> False
                Right dense ->
                    entries dense == [0, 0, 0]
                        && all
                            (rightClose 0 . Visit.infiniteProbabilityGivenInitialState terminalChain 1)
                            states
            && case visitTotalExpectationByState terminalChain 1 of
                Left _ -> False
                Right actual ->
                    and (zipWith expectationClose visitExpectations actual)
                        && and
                            [ either
                                (const False)
                                (expectationClose expected)
                                (Visit.totalExpectationGivenInitialState terminalChain 1 state)
                            | (state, expected) <- zip states visitExpectations
                            ]

spec :: Spec
spec = do
    describe "canonical finite-horizon differential baseline" $ do
        prop "all finite and bounded queries match path enumeration (random @3)" $
            forAll (genTransitionMatrix @3) $ \rawMatrix ->
                case mkTransitionMatrix rawMatrix of
                    Left problem -> counterexample (show problem) False
                    Right matrix -> property (finiteAndBoundedChecks matrix)

    describe "canonical infinite-horizon differential baseline" $ do
        it "all eventual, race, expectation, and total-visit queries match a completed path law" $
            terminalChecks
