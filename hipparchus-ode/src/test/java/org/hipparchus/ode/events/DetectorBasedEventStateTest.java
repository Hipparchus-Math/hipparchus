/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.ode.events;


import java.lang.reflect.Field;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.solvers.BracketedUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.EquationsMapper;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.SecondaryODE;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.ode.nonstiff.LutherIntegrator;
import org.hipparchus.ode.nonstiff.interpolators.ClassicalRungeKuttaStateInterpolator;
import org.hipparchus.ode.sampling.DummyStepInterpolator;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.mockito.Mockito;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class DetectorBasedEventStateTest {

    // Unit test reproducing https://gitlab.orekit.org/orekit/orekit/-/issues/1808
    @Test
    void testEpochComparisonAtLeastSignificantBit() throws NoSuchFieldException, IllegalAccessException {
        // Epoch of event
        final double eventTime = 17016.237999999998;

        // Get the interpolated state at event time
        // It will return globalCurrent at 17016.238 sec since the difference between current and previous state times is smaller than the least bit
        final ODEStateAndDerivative globalCurrent = new ODEStateAndDerivative(17016.238, new double[2], new double[2]);
        final ODEStateAndDerivative globalPrevious = new ODEStateAndDerivative(17016.237999999998, new double[2], new double[2]);
        final ClassicalRungeKuttaStateInterpolator interpolator = new ClassicalRungeKuttaStateInterpolator(true, new double[2][2], globalPrevious, globalCurrent,
                                                                                                     globalPrevious, globalCurrent, null);
        final ODEStateAndDerivative interpolatedState = interpolator.getInterpolatedState(eventTime);
        Assertions.assertEquals(interpolatedState.getTime(), globalCurrent.getTime());
        Assertions.assertNotEquals(interpolatedState.getTime(), globalPrevious.getTime());

        // Configure the event state (failing before the fix)
        // Since detecting the event causing the numerical issue is tricky; we access the private field to simplify the workflow and directly set the necessary values causing the issue
        final DetectorBasedEventState es = new DetectorBasedEventState(new CloseEventsGenerator(16436.238, 17016.238, 720, 1e-10, 100));
        final Field pendingEvent = DetectorBasedEventState.class.getDeclaredField("pendingEvent");
        pendingEvent.setAccessible(true);
        pendingEvent.set(es, true);
        final Field pendingEventTime = DetectorBasedEventState.class.getDeclaredField("pendingEventTime");
        pendingEventTime.setAccessible(true);
        pendingEventTime.set(es, eventTime);
        final Field afterG = DetectorBasedEventState.class.getDeclaredField("afterG");
        afterG.setAccessible(true);
        afterG.set(es, 0.0); // Dummy value (this value is not interesting in that case)

        // Action & verify
        Assertions.assertNotNull(es.doEvent(interpolatedState));
    }

    @Test
    void testDoEventThrowsIfTimeMismatch() throws NoSuchFieldException, IllegalAccessException {
        // Initialization
        final ODEEventDetector detector = new DummyDetector();
        final DetectorBasedEventState eventState = new DetectorBasedEventState(detector);
        final Field pendingEvent = DetectorBasedEventState.class.getDeclaredField("pendingEvent");
        pendingEvent.setAccessible(true);
        pendingEvent.set(eventState, true);
        final Field pendingEventTime = DetectorBasedEventState.class.getDeclaredField("pendingEventTime");
        pendingEventTime.setAccessible(true);
        pendingEventTime.set(eventState, 1.0);
        final ODEStateAndDerivative state = new ODEStateAndDerivative(1.0001, new double[]{0.0}, new double[]{0.0});
        // Action & verify
        Assertions.assertThrows(org.hipparchus.exception.MathRuntimeException.class, () -> eventState.doEvent(state));
    }

    @Test
    void testInitSetsInitialValues() {
        // Initialize
        final ODEEventDetector detector = new DummyDetector();
        final DetectorBasedEventState eventState = new DetectorBasedEventState(detector);
        final ODEStateAndDerivative s0 = new ODEStateAndDerivative(0.0, new double[1], new double[1]);
        // Action
        eventState.init(s0, 10.0);
        // Verify
        Assertions.assertNotNull(eventState.getEventDetector());
        Assertions.assertEquals(detector, eventState.getEventDetector());
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, eventState.getEventTime());
    }

    @Test
    void testEvaluateStepNoEvent() {
        // Initialize
        final ODEEventDetector detector = new DummyDetector();
        final DetectorBasedEventState eventState = new DetectorBasedEventState(detector);
        final ODEStateAndDerivative s0 = new ODEStateAndDerivative(0.0, new double[1], new double[1]);
        final ODEStateAndDerivative s1 = new ODEStateAndDerivative(1.0, new double[1], new double[1]);
        final ODEStateInterpolator interpolator = Mockito.mock(ODEStateInterpolator.class);
        Mockito.when(interpolator.isForward()).thenReturn(true);
        Mockito.when(interpolator.getPreviousState()).thenReturn(s0);
        Mockito.when(interpolator.getCurrentState()).thenReturn(s1);
        // Action
        eventState.init(s0, 1.0);
        eventState.reinitializeBegin(interpolator);
        final boolean hasEvent = eventState.evaluateStep(interpolator);
        // Verify
        Assertions.assertFalse(hasEvent);
    }

    // JIRA: MATH-322
    @Test
    void closeEvents() throws MathIllegalArgumentException, MathIllegalStateException {

        final double r1  = 90.0;
        final double r2  = 135.0;
        final double gap = r2 - r1;

        final double tolerance = 0.1;
        DetectorBasedEventState es = new DetectorBasedEventState(new CloseEventsGenerator(r1, r2, 1.5 * gap, tolerance, 100));
        EquationsMapper mapper = new ExpandableODE(new OrdinaryDifferentialEquation() {
            @Override
            public int getDimension() {
                return 0;
            }
            @Override
            public double[] computeDerivatives(double t, double[] y) {
                return new double[0];
            }
        }).getMapper();

        double[] a = new double[0];
        ODEStateAndDerivative osdLongBefore = new ODEStateAndDerivative(r1 - 1.5 * gap, a, a);
        ODEStateAndDerivative osBefore      = new ODEStateAndDerivative(r1 - 0.5 * gap, a, a);
        ODEStateInterpolator interpolatorA  = new DummyStepInterpolator(true,
                                                                        osdLongBefore, osBefore,
                                                                        osdLongBefore, osBefore,
                                                                        mapper);
        es.reinitializeBegin(interpolatorA);
        assertFalse(es.evaluateStep(interpolatorA));

        ODEStateAndDerivative osdBetween    = new ODEStateAndDerivative(0.5 * (r1 + r2), a, a);
        ODEStateInterpolator interpolatorB  = new DummyStepInterpolator(true,
                                                                        osBefore, osdBetween,
                                                                        osBefore, osdBetween,
                                                                        mapper);
        assertTrue(es.evaluateStep(interpolatorB));
        assertEquals(r1, es.getEventTime(), tolerance);
        ODEStateAndDerivative osdAtEvent    = new ODEStateAndDerivative(es.getEventTime(), a, a);
        es.doEvent(osdAtEvent);

        ODEStateAndDerivative osdAfterSecond = new ODEStateAndDerivative(r2 + 0.4 * gap, a, a);
        ODEStateInterpolator interpolatorC  = new DummyStepInterpolator(true,
                                                                        osdAtEvent, osdAfterSecond,
                                                                        osdAtEvent, osdAfterSecond,
                                                                        mapper);
        assertTrue(es.evaluateStep(interpolatorC));
        assertEquals(r2, es.getEventTime(), tolerance);

    }

    // Jira: MATH-695
    @Test
    void testIssue695()
        throws MathIllegalArgumentException, MathIllegalStateException {

        OrdinaryDifferentialEquation equation = new OrdinaryDifferentialEquation() {
            @Override
            public int getDimension() {
                return 1;
            }
            @Override
            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { 1.0 };
            }
        };

        DormandPrince853Integrator integrator = new DormandPrince853Integrator(0.001, 1000, 1.0e-14, 1.0e-14);
        integrator.addEventDetector(new ResettingEvent(10.99, 0.1, 1.0e-9, 1000));
        integrator.addEventDetector(new ResettingEvent(11.01, 0.1, 1.0e-9, 1000));
        integrator.setInitialStepSize(3.0);

        double target = 30.0;
        ODEStateAndDerivative finalState =
                        integrator.integrate(equation, new ODEState(0.0, new double[1]), target);
        assertEquals(target, finalState.getTime(), 1.0e-10);
        assertEquals(32.0, finalState.getPrimaryState()[0], 1.0e-10);

    }

    // Jira: MATH-965
    @Test
    void testIssue965()
        throws MathIllegalArgumentException, MathIllegalStateException {

        ExpandableODE equation = new ExpandableODE(new OrdinaryDifferentialEquation() {
            @Override
            public int getDimension() {
                return 1;
            }
            @Override
            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { 2.0 };
            }
        });
        int index = equation.addSecondaryEquations(new SecondaryODE() {
            @Override
            public int getDimension() {
                return 1;
            }
            @Override
            public double[] computeDerivatives(double t, double[] primary,
                                           double[] primaryDot, double[] secondary) {
                return new double[] { -3.0 };
            }
        });
        assertEquals(1, index);

        DormandPrince853Integrator integrator = new DormandPrince853Integrator(0.001, 1000, 1.0e-14, 1.0e-14);
        integrator.addEventDetector(new SecondaryStateEvent(index, -3.0, 0.1, 1.0e-9, 1000));
        integrator.setInitialStepSize(3.0);

        ODEState initialState = new ODEState(0.0,
                                             new double[] { 0.0 },
                                             new double[][] { { 0.0 } });
        ODEStateAndDerivative finalState = integrator.integrate(equation, initialState, 30.0);
        assertEquals( 1.0, finalState.getTime(), 1.0e-10);
        assertEquals( 2.0, finalState.getPrimaryState()[0], 1.0e-10);
        assertEquals(-3.0, finalState.getSecondaryState(index)[0], 1.0e-10);

    }

    @ParameterizedTest
    @ValueSource(booleans = {true, false})
    void testAdaptableInterval(final boolean isForward) {
        // GIVEN
        final TestDetector detector = new TestDetector();
        final DetectorBasedEventState eventState = new DetectorBasedEventState(detector);
        final ODEStateInterpolator mockedInterpolator = Mockito.mock(ODEStateInterpolator.class);
        final ODEStateAndDerivative stateAndDerivative1 = getStateAndDerivative(1);
        final ODEStateAndDerivative stateAndDerivative2 = getStateAndDerivative(-1);
        if (isForward) {
            Mockito.when(mockedInterpolator.getCurrentState()).thenReturn(stateAndDerivative1);
            Mockito.when(mockedInterpolator.getPreviousState()).thenReturn(stateAndDerivative2);
        } else {
            Mockito.when(mockedInterpolator.getCurrentState()).thenReturn(stateAndDerivative2);
            Mockito.when(mockedInterpolator.getPreviousState()).thenReturn(stateAndDerivative1);
        }
        Mockito.when(mockedInterpolator.isForward()).thenReturn(isForward);
        // WHEN
        eventState.evaluateStep(mockedInterpolator);
        // THEN
        if (isForward) {
            Assertions.assertEquals(1, detector.triggeredForward);
            Assertions.assertEquals(0, detector.triggeredBackward);
        } else {
            Assertions.assertEquals(0, detector.triggeredForward);
            Assertions.assertEquals(1, detector.triggeredBackward);
        }
    }

    @Test
    void testEventsCloserThanThreshold()
            throws MathIllegalArgumentException, MathIllegalStateException {

        OrdinaryDifferentialEquation equation = new OrdinaryDifferentialEquation() {

            @Override
            public int getDimension() {
                return 1;
            }

            @Override
            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { 1.0 };
            }
        };

        LutherIntegrator integrator = new LutherIntegrator(20.0);
        CloseEventsGenerator eventsGenerator =
                new CloseEventsGenerator(9.0 - 1.0 / 128, 9.0 + 1.0 / 128, 1.0, 0.02, 1000);
        integrator.addEventDetector(eventsGenerator);
        double tEnd = integrator.integrate(equation, new ODEState(0.0, new double[1]), 100.0).getTime();
        assertEquals( 2, eventsGenerator.getCount());
        assertEquals( 9.0 + 1.0 / 128, tEnd, 1.0 / 32.0);

    }

    private static class SecondaryStateEvent implements ODEEventDetector {

        private final AdaptableInterval             maxCheck;
        private final int                           maxIter;
        private final BracketingNthOrderBrentSolver solver;
        private final int                           index;
        private final double                        target;

        public SecondaryStateEvent(final int index, final double target,
                                   final double maxCheck, final double threshold, final int maxIter) {
            this.maxCheck  = (s, isForward) -> maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new BracketingNthOrderBrentSolver(0, threshold, 0, 5);
            this.index     = index;
            this.target    = target;
        }

        @Override
        public AdaptableInterval getMaxCheckInterval() {
            return maxCheck;
        }

        @Override
        public int getMaxIterationCount() {
            return maxIter;
        }

        @Override
        public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
            return solver;
        }

        /** {@inheritDoc} */
        @Override
        public ODEEventHandler getHandler() {
            return (state, detector, increasing) -> Action.STOP;
        }

        @Override
        public double g(ODEStateAndDerivative s) {
            return s.getSecondaryState(index)[0] - target;
        }

    }

    private static class DummyDetector implements ODEEventDetector {
        @Override
        public AdaptableInterval getMaxCheckInterval() {
            return (state, isForward) -> 1.0;
        }
        @Override
        public int getMaxIterationCount() {
            return 10;
        }
        @Override
        public BracketingNthOrderBrentSolver getSolver() {
            return new BracketingNthOrderBrentSolver();
        }
        @Override
        public ODEEventHandler getHandler() {
            return (state, detector, increasing) -> Action.CONTINUE;
        }
        @Override
        public double g(ODEStateAndDerivative state) {
            return 1.0;
        }
    }

    private static class ResettingEvent implements ODEEventDetector {

        private static double lastTriggerTime = Double.NEGATIVE_INFINITY;
        private final AdaptableInterval             maxCheck;
        private final int                           maxIter;
        private final BracketingNthOrderBrentSolver solver;
        private final double                        tEvent;

        public ResettingEvent(final double tEvent,
                              final double maxCheck, final double threshold, final int maxIter) {
            this.maxCheck  = (s, isForward) -> maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new BracketingNthOrderBrentSolver(0, threshold, 0, 5);
            this.tEvent    = tEvent;
        }

        @Override
        public AdaptableInterval getMaxCheckInterval() {
            return maxCheck;
        }

        @Override
        public int getMaxIterationCount() {
            return maxIter;
        }

        @Override
        public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
            return solver;
        }

        @Override
        public double g(ODEStateAndDerivative s) {
            // the bug corresponding to issue 695 causes the g function
            // to be called at obsolete times t despite an event
            // occurring later has already been triggered.
            // When this occurs, the following assertion is violated
            assertTrue(s.getTime() >= lastTriggerTime,
                    "going backard in time! (" + s.getTime() + " < " + lastTriggerTime + ")");
            return s.getTime() - tEvent;
        }

        @Override
        public ODEEventHandler getHandler() {
            return new ODEEventHandler() {
                @Override
                public Action eventOccurred(ODEStateAndDerivative s, ODEEventDetector detector, boolean increasing) {
                    // remember in a class variable when the event was triggered
                    lastTriggerTime = s.getTime();
                    return Action.RESET_STATE;
                }

                @Override
                public ODEStateAndDerivative resetState(ODEEventDetector detector, ODEStateAndDerivative s) {
                    double[] y = s.getPrimaryState();
                    y[0] += 1.0;
                    return new ODEStateAndDerivative(s.getTime(), y, s.getPrimaryDerivative());
                }
            };
        }

    }

    private static class CloseEventsGenerator implements ODEEventDetector {

        private final AdaptableInterval             maxCheck;
        private final int                           maxIter;
        private final BracketingNthOrderBrentSolver solver;
        final double                                r1;
        final double                                r2;
        int                                         count;

        public CloseEventsGenerator(final double r1, final double r2,
                                    final double maxCheck, final double threshold, final int maxIter) {
            this.maxCheck  = (s, isForward) -> maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new BracketingNthOrderBrentSolver(0, threshold, 0, 5);
            this.r1        = r1;
            this.r2        = r2;
            this.count     = 0;
        }

        @Override
        public AdaptableInterval getMaxCheckInterval() {
            return maxCheck;
        }

        @Override
        public int getMaxIterationCount() {
            return maxIter;
        }

        @Override
        public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
            return solver;
        }

        @Override
        public double g(ODEStateAndDerivative s) {
            return (s.getTime() - r1) * (r2 - s.getTime());
        }

        @Override
        public ODEEventHandler getHandler() {
            return (state, detector, increasing) -> ++count < 2 ? Action.CONTINUE : Action.STOP;
        }

        public int getCount() {
            return count;
        }

    }

    private static ODEStateAndDerivative getStateAndDerivative(final double time) {
        return new ODEStateAndDerivative(time, new double[] {time}, new double[1]);
    }

    private static class TestDetector implements ODEEventDetector {

        int triggeredForward = 0;
        int triggeredBackward = 0;

        @Override
        public AdaptableInterval getMaxCheckInterval() {
            return (state, isForward) -> {
                if (isForward) {
                    triggeredForward++;
                } else {
                    triggeredBackward++;
                }
                return 1.;
            };
        }

        @Override
        public int getMaxIterationCount() {
            return 10;
        }

        @Override
        public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
            return new BracketingNthOrderBrentSolver();
        }

        @Override
        public ODEEventHandler getHandler() {
            return null;
        }

        @Override
        public double g(ODEStateAndDerivative state) {
            return 0.;
        }
    }

}
