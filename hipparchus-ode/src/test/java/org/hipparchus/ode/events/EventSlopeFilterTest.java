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

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.solvers.BracketedUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;

import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well19937a;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;
import org.mockito.Mockito;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class EventSlopeFilterTest {

    @Test
    void testInit() {
        // GIVEN
        final TestDetector detector = new TestDetector();
        final EventSlopeFilter<?> eventSlopeFilter = new EventSlopeFilter<>(detector,
                FilterType.TRIGGER_ONLY_DECREASING_EVENTS);
        // WHEN
        eventSlopeFilter.init(Mockito.mock(ODEStateAndDerivative.class), 1);
        // THEN
        assertTrue(detector.initialized);
        assertFalse(detector.resetted);
    }

    @ParameterizedTest
    @EnumSource(FilterType.class)
    void testReset(final FilterType type) {
        // GIVEN
        final TestDetector detector = new TestDetector();
        final EventSlopeFilter<?> eventSlopeFilter = new EventSlopeFilter<>(detector, type);
        // WHEN
        eventSlopeFilter.reset(Mockito.mock(ODEStateAndDerivative.class), 1);
        // THEN
        assertTrue(detector.resetted);
        assertFalse(detector.initialized);
    }

    private static class TestDetector extends AbstractODEDetector<TestDetector> {
        boolean initialized = false;
        boolean resetted = false;

        protected TestDetector() {
            super((state, isForward) -> 0, 1, new BracketingNthOrderBrentSolver(), (s, e, d) -> Action.CONTINUE);
        }

        @Override
        public void init(ODEStateAndDerivative s0, double t) {
            super.init(s0, t);
            initialized = true;
        }

        @Override
        public void reset(ODEStateAndDerivative intermediateState, double finalTime) {
            super.reset(intermediateState, finalTime);
            resetted = true;
        }

        @Override
        protected TestDetector create(AdaptableInterval newMaxCheck, int newmaxIter, BracketedUnivariateSolver<UnivariateFunction> newSolver, ODEEventHandler newHandler) {
            return null;
        }

        @Override
        public double g(ODEStateAndDerivative state) {
            return 1;
        }
    }

    @Test
    void testHistoryIncreasingForward() {

        // start point: g > 0
        testHistory(FilterType.TRIGGER_ONLY_INCREASING_EVENTS,
                0.5 * FastMath.PI, 30.5 * FastMath.PI, FastMath.PI, -1);

        // start point: g = 0
        testHistory(FilterType.TRIGGER_ONLY_INCREASING_EVENTS,
                0, 30.5 * FastMath.PI, FastMath.PI, -1);

        // start point: g < 0
        testHistory(FilterType.TRIGGER_ONLY_INCREASING_EVENTS,
                1.5 * FastMath.PI, 30.5 * FastMath.PI, FastMath.PI, +1);

    }

    @Test
    void testHistoryIncreasingBackward() {

        // start point: g > 0
        testHistory(FilterType.TRIGGER_ONLY_INCREASING_EVENTS,
                    0.5 * FastMath.PI, -30.5 * FastMath.PI, FastMath.PI, -1);

        // start point: g = 0
        testHistory(FilterType.TRIGGER_ONLY_INCREASING_EVENTS,
                    0, -30.5 * FastMath.PI, FastMath.PI, +1);

        // start point: g < 0
        testHistory(FilterType.TRIGGER_ONLY_INCREASING_EVENTS,
                    1.5 * FastMath.PI, -30.5 * FastMath.PI, FastMath.PI, -1);

    }

    @Test
    void testHistoryDecreasingForward() {

        // start point: g > 0
        testHistory(FilterType.TRIGGER_ONLY_DECREASING_EVENTS,
                    0.5 * FastMath.PI, 30.5 * FastMath.PI, 0, +1);

        // start point: g = 0
        testHistory(FilterType.TRIGGER_ONLY_DECREASING_EVENTS,
                    0, 30.5 * FastMath.PI, 0, +1);

        // start point: g < 0
        testHistory(FilterType.TRIGGER_ONLY_DECREASING_EVENTS,
                    1.5 * FastMath.PI, 30.5 * FastMath.PI, 0, +1);

    }

    @Test
    void testHistoryDecreasingBackward() {

        // start point: g > 0
        testHistory(FilterType.TRIGGER_ONLY_DECREASING_EVENTS,
                    0.5 * FastMath.PI, -30.5 * FastMath.PI, 0, -1);

        // start point: g = 0
        testHistory(FilterType.TRIGGER_ONLY_DECREASING_EVENTS,
                    0, -30.5 * FastMath.PI, 0, -1);

        // start point: g < 0
        testHistory(FilterType.TRIGGER_ONLY_DECREASING_EVENTS,
                    1.5 * FastMath.PI, -30.5 * FastMath.PI, 0, +1);

    }

    private void testHistory(FilterType type, double t0, double t1, double refSwitch, double signEven) {
        Event onlyIncreasing = new Event(Double.POSITIVE_INFINITY, 1.0e-10, 1000, false, true);
        EventSlopeFilter<Event> eventFilter = new EventSlopeFilter<>(onlyIncreasing, type);
        eventFilter.init(new ODEStateAndDerivative(t0,
                                                   new double[] { FastMath.sin(t0),  FastMath.cos(t0) },
                                                   new double[] { FastMath.cos(t0), -FastMath.sin(t0) }),
                         t1);

        // first pass to set up switches history for a long period
        double h = FastMath.copySign(0.05, t1 - t0);
        double n = (int) FastMath.floor((t1 - t0) / h);
        for (int i = 0; i < n; ++i) {
            double t = t0 + i * h;
            eventFilter.g(new ODEStateAndDerivative(t,
                                                    new double[] { FastMath.sin(t),  FastMath.cos(t) },
                                                    new double[] { FastMath.cos(t), -FastMath.sin(t) }));
        }

        // verify old events are preserved, even if randomly accessed
        RandomGenerator rng = new Well19937a(0xb0e7401265af8cd3l);
        for (int i = 0; i < 5000; i++) {
            double t = t0 + (t1 - t0) * rng.nextDouble();
            double g = eventFilter.g(new ODEStateAndDerivative(t,
                                                               new double[] { FastMath.sin(t),  FastMath.cos(t) },
                                                               new double[] { FastMath.cos(t), -FastMath.sin(t) }));
            int turn = (int) FastMath.floor((t - refSwitch) / (2 * FastMath.PI));
            if (turn % 2 == 0) {
                assertEquals( signEven * FastMath.sin(t), g, 1.0e-10);
            } else {
                assertEquals(-signEven * FastMath.sin(t), g, 1.0e-10);
            }
        }

    }

    @Test
    void testIncreasingOnly()
        throws MathIllegalArgumentException, MathIllegalStateException {
        ODEIntegrator integrator = new DormandPrince853Integrator(1.0e-3, 100.0, 1e-7, 1e-7);
        Event allEvents = new Event(0.1, 1.0e-7, 100, true, true);
        integrator.addEventDetector(allEvents);
        Event onlyIncreasing = new Event(0.1, 1.0e-7, 100, false, true);
        integrator.addEventDetector(new EventSlopeFilter<>(onlyIncreasing, FilterType.TRIGGER_ONLY_INCREASING_EVENTS));
        double t0 = 0.5 * FastMath.PI;
        double tEnd = 5.5 * FastMath.PI;
        double[] y = { 0.0, 1.0 };
        assertEquals(tEnd,
                            integrator.integrate(new SineCosine(), new ODEState(t0, y), tEnd).getTime(),
                            1.0e-7);

        assertEquals(5, allEvents.getEventCount());
        assertEquals(2, onlyIncreasing.getEventCount());

    }

    @Test
    void testDecreasingOnly()
        throws MathIllegalArgumentException, MathIllegalStateException {
        ODEIntegrator integrator = new DormandPrince853Integrator(1.0e-3, 100.0, 1e-7, 1e-7);
        Event allEvents = new Event(0.1, 1.0e-7, 1000, true, true);
        integrator.addEventDetector(allEvents);
        Event onlyDecreasing = new Event(0.1, 1.0e-7, 1000, true, false);
        integrator.addEventDetector(new EventSlopeFilter<>(onlyDecreasing, FilterType.TRIGGER_ONLY_DECREASING_EVENTS));
        double t0 = 0.5 * FastMath.PI;
        double tEnd = 5.5 * FastMath.PI;
        double[] y = { 0.0, 1.0 };
        assertEquals(tEnd,
                            integrator.integrate(new SineCosine(), new ODEState(t0, y), tEnd).getTime(),
                            1.0e-7);

        assertEquals(5, allEvents.getEventCount());
        assertEquals(3, onlyDecreasing.getEventCount());

    }

    @Test
    void testTwoOppositeFilters()
        throws MathIllegalArgumentException, MathIllegalStateException {
        ODEIntegrator integrator = new DormandPrince853Integrator(1.0e-3, 100.0, 1e-7, 1e-7);
        Event allEvents = new Event(0.1, 1.0e-7, 100, true, true);
        integrator.addEventDetector(allEvents);
        Event onlyIncreasing = new Event(0.1, 1.0e-7, 100, false, true);
        integrator.addEventDetector(new EventSlopeFilter<>(onlyIncreasing, FilterType.TRIGGER_ONLY_INCREASING_EVENTS));
        Event onlyDecreasing = new Event(0.1, 1.0e-7, 100, true, false);
        integrator.addEventDetector(new EventSlopeFilter<>(onlyDecreasing, FilterType.TRIGGER_ONLY_DECREASING_EVENTS));
        double t0 = 0.5 * FastMath.PI;
        double tEnd = 5.5 * FastMath.PI;
        double[] y = { 0.0, 1.0 };
        assertEquals(tEnd,
                            integrator.integrate(new SineCosine(), new ODEState(t0, y), tEnd).getTime(),
                            1.0e-7);

        assertEquals(5, allEvents.getEventCount());
        assertEquals(2, onlyIncreasing.getEventCount());
        assertEquals(3, onlyDecreasing.getEventCount());

    }

    private static class SineCosine implements OrdinaryDifferentialEquation {
        public int getDimension() {
            return 2;
        }

        public double[] computeDerivatives(double t, double[] y) {
            return new double[] { y[1], -y[0] };
        }
    }

    /** State events for this unit test. */
    protected static class Event implements ODEEventDetector {

        private final AdaptableInterval                             maxCheck;
        private final int                                           maxIter;
        private final BracketedUnivariateSolver<UnivariateFunction> solver;
        private final boolean                                       expectDecreasing;
        private final boolean                                       expectIncreasing;
        private int                                                 eventCount;

        public Event(final double maxCheck, final double threshold, final int maxIter,
                     boolean expectDecreasing, boolean expectIncreasing) {
            this.maxCheck         = (s, isForward) -> maxCheck;
            this.maxIter          = maxIter;
            this.solver           = new BracketingNthOrderBrentSolver(0, threshold, 0, 5);
            this.expectDecreasing = expectDecreasing;
            this.expectIncreasing = expectIncreasing;
        }

        public AdaptableInterval getMaxCheckInterval() {
            return maxCheck;
        }

        public int getMaxIterationCount() {
            return maxIter;
        }

        public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
            return solver;
        }

        public int getEventCount() {
            return eventCount;
        }

        @Override
        public void init(ODEStateAndDerivative s0, double t) {
            eventCount = 0;
        }

        public double g(ODEStateAndDerivative s) {
            return s.getPrimaryState()[0];
        }

        public ODEEventHandler getHandler() {
            return (ODEStateAndDerivative s, ODEEventDetector detector, boolean increasing) -> {
                if (increasing) {
                    assertTrue(expectIncreasing);
                } else {
                    assertTrue(expectDecreasing);
                }
                eventCount++;
                return Action.RESET_STATE;
            };
        }

    }

}
