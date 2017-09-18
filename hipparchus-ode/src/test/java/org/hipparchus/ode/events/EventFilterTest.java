/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
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
import org.junit.Assert;
import org.junit.Test;

public class EventFilterTest {

    @Test
    public void testHistoryIncreasingForward() {

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
    public void testHistoryIncreasingBackward() {

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
    public void testHistoryDecreasingForward() {

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
    public void testHistoryDecreasingBackward() {

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

    public void testHistory(FilterType type, double t0, double t1, double refSwitch, double signEven) {
        Event onlyIncreasing = new Event(false, true);
        EventFilter eventFilter =
                new EventFilter(onlyIncreasing, type);
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
                Assert.assertEquals( signEven * FastMath.sin(t), g, 1.0e-10);
            } else {
                Assert.assertEquals(-signEven * FastMath.sin(t), g, 1.0e-10);
            }
        }

    }

    @Test
    public void testIncreasingOnly()
        throws MathIllegalArgumentException, MathIllegalStateException {
        double e = 1e-15;
        ODEIntegrator integrator = new DormandPrince853Integrator(1.0e-3, 100.0, 1e-7, 1e-7);
        Event allEvents = new Event(true, true);
        integrator.addEventHandler(allEvents, 0.1, e, 1000,
                                   new BracketingNthOrderBrentSolver(1.0e-7, 5));
        Event onlyIncreasing = new Event(false, true);
        integrator.addEventHandler(new EventFilter(onlyIncreasing,
                                                   FilterType.TRIGGER_ONLY_INCREASING_EVENTS),
                                   0.1, e, 100,
                                   new BracketingNthOrderBrentSolver(1.0e-7, 5));
        double t0 = 0.5 * FastMath.PI;
        double tEnd = 5.5 * FastMath.PI;
        double[] y = { 0.0, 1.0 };
        Assert.assertEquals(tEnd,
                            integrator.integrate(new SineCosine(), new ODEState(t0, y), tEnd).getTime(),
                            1.0e-7);

        Assert.assertEquals(5, allEvents.getEventCount());
        Assert.assertEquals(2, onlyIncreasing.getEventCount());

    }

    @Test
    public void testDecreasingOnly()
        throws MathIllegalArgumentException, MathIllegalStateException {
        double e = 1e-15;
        ODEIntegrator integrator = new DormandPrince853Integrator(1.0e-3, 100.0, 1e-7, 1e-7);
        Event allEvents = new Event(true, true);
        integrator.addEventHandler(allEvents, 0.1, e, 1000,
                                   new BracketingNthOrderBrentSolver(1.0e-7, 5));
        Event onlyDecreasing = new Event(true, false);
        integrator.addEventHandler(new EventFilter(onlyDecreasing,
                                                   FilterType.TRIGGER_ONLY_DECREASING_EVENTS),
                                   0.1, e, 1000,
                                   new BracketingNthOrderBrentSolver(1.0e-7, 5));
        double t0 = 0.5 * FastMath.PI;
        double tEnd = 5.5 * FastMath.PI;
        double[] y = { 0.0, 1.0 };
        Assert.assertEquals(tEnd,
                            integrator.integrate(new SineCosine(), new ODEState(t0, y), tEnd).getTime(),
                            1.0e-7);

        Assert.assertEquals(5, allEvents.getEventCount());
        Assert.assertEquals(3, onlyDecreasing.getEventCount());

    }

    @Test
    public void testTwoOppositeFilters()
        throws MathIllegalArgumentException, MathIllegalStateException {
        double e = 1e-15;
        ODEIntegrator integrator = new DormandPrince853Integrator(1.0e-3, 100.0, 1e-7, 1e-7);
        Event allEvents = new Event(true, true);
        integrator.addEventHandler(allEvents, 0.1, e, 1000,
                                   new BracketingNthOrderBrentSolver(1.0e-7, 5));
        Event onlyIncreasing = new Event(false, true);
        integrator.addEventHandler(new EventFilter(onlyIncreasing,
                                                   FilterType.TRIGGER_ONLY_INCREASING_EVENTS),
                                   0.1, e, 1000,
                                   new BracketingNthOrderBrentSolver(1.0e-7, 5));
        Event onlyDecreasing = new Event(true, false);
        integrator.addEventHandler(new EventFilter(onlyDecreasing,
                                                   FilterType.TRIGGER_ONLY_DECREASING_EVENTS),
                                   0.1, e, 1000,
                                   new BracketingNthOrderBrentSolver(1.0e-7, 5));
        double t0 = 0.5 * FastMath.PI;
        double tEnd = 5.5 * FastMath.PI;
        double[] y = { 0.0, 1.0 };
        Assert.assertEquals(tEnd,
                            integrator.integrate(new SineCosine(), new ODEState(t0, y), tEnd).getTime(),
                            1.0e-7);

        Assert.assertEquals(5, allEvents.getEventCount());
        Assert.assertEquals(2, onlyIncreasing.getEventCount());
        Assert.assertEquals(3, onlyDecreasing.getEventCount());

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
    protected static class Event implements ODEEventHandler {

        private final boolean expectDecreasing;
        private final boolean expectIncreasing;
        private int eventCount;

        public Event(boolean expectDecreasing, boolean expectIncreasing) {
            this.expectDecreasing = expectDecreasing;
            this.expectIncreasing = expectIncreasing;
        }

        public int getEventCount() {
            return eventCount;
        }

        public void init(ODEStateAndDerivative s0, double t) {
            eventCount = 0;
        }

        public double g(ODEStateAndDerivative s) {
            return s.getPrimaryState()[0];
        }

        public Action eventOccurred(ODEStateAndDerivative s, boolean increasing) {
            if (increasing) {
                Assert.assertTrue(expectIncreasing);
            } else {
                Assert.assertTrue(expectDecreasing);
            }
            eventCount++;
            return Action.RESET_STATE;
        }

    }
}
