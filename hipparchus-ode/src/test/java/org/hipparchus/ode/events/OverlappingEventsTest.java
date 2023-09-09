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

import java.util.ArrayList;
import java.util.List;

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
import org.junit.Assert;
import org.junit.Test;

/** Tests for overlapping state events. Also tests an event function that does
 * not converge to zero, but does have values of opposite sign around its root.
 */
public class OverlappingEventsTest implements OrdinaryDifferentialEquation {

    /** Expected event times for first event. */
    private static final double[] EVENT_TIMES1 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
                                                  7.0, 8.0, 9.0};

    /** Expected event times for second event. */
    private static final double[] EVENT_TIMES2 = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0,
                                                  3.5, 4.0, 4.5, 5.0, 5.5, 6.0,
                                                  6.5, 7.0, 7.5, 8.0, 8.5, 9.0,
                                                  9.5};

    /** Test for events that occur at the exact same time, but due to numerical
     * calculations occur very close together instead. Uses event type 0. See
     * {@link org.hipparchus.ode.events.ODEEventHandler#g(double, double[])
     * ODEEventHandler.g(double, double[])}.
     */
    @Test
    public void testOverlappingEvents0()
        throws MathIllegalArgumentException, MathIllegalStateException {
        test(0);
    }

    /** Test for events that occur at the exact same time, but due to numerical
     * calculations occur very close together instead. Uses event type 1. See
     * {@link org.hipparchus.ode.events.ODEEventHandler#g(double, double[])
     * ODEEventHandler.g(double, double[])}.
     */
    @Test
    public void testOverlappingEvents1()
        throws MathIllegalArgumentException, MathIllegalStateException {
        test(1);
    }

    /** Test for events that occur at the exact same time, but due to numerical
     * calculations occur very close together instead.
     * @param eventType the type of events to use. See
     * {@link org.hipparchus.ode.events.ODEEventHandler#g(double, double[])
     * ODEEventHandler.g(double, double[])}.
     */
    public void test(int eventType)
        throws MathIllegalArgumentException, MathIllegalStateException {
        double e = 1e-15;
        ODEIntegrator integrator = new DormandPrince853Integrator(e, 100.0, 1e-7, 1e-7);
        ODEEventDetector evt1 = new Event(0.1, e, 999, 0, eventType);
        ODEEventDetector evt2 = new Event(0.1, e, 999, 1, eventType);
        integrator.addEventDetector(evt1);
        integrator.addEventDetector(evt2);
        double t = 0.0;
        double tEnd = 9.75;
        double[] y = {0.0, 0.0};
        List<Double> events1 = new ArrayList<Double>();
        List<Double> events2 = new ArrayList<Double>();
        while (t < tEnd) {
            final ODEStateAndDerivative finalState =
                            integrator.integrate(this, new ODEState(t, y), tEnd);
            t = finalState.getTime();
            y = finalState.getPrimaryState();
            //System.out.println("t=" + t + ",\t\ty=[" + y[0] + "," + y[1] + "]");
            if (y[0] >= 1.0) {
                y[0] = 0.0;
                events1.add(t);
                //System.out.println("Event 1 @ t=" + t);
            }
            if (y[1] >= 1.0) {
                y[1] = 0.0;
                events2.add(t);
                //System.out.println("Event 2 @ t=" + t);
            }
        }
        Assert.assertEquals(EVENT_TIMES1.length, events1.size());
        Assert.assertEquals(EVENT_TIMES2.length, events2.size());
        for(int i = 0; i < EVENT_TIMES1.length; i++) {
            Assert.assertEquals(EVENT_TIMES1[i], events1.get(i), 1e-7);
        }
        for(int i = 0; i < EVENT_TIMES2.length; i++) {
            Assert.assertEquals(EVENT_TIMES2[i], events2.get(i), 1e-7);
        }
        //System.out.println();
    }

    /** {@inheritDoc} */
    public int getDimension() {
        return 2;
    }

    /** {@inheritDoc} */
    public double[] computeDerivatives(double t, double[] y) {
        return new double[] { 1.0, 2.0 };
    }

    /** State events for this unit test. */
    private class Event implements ODEEventDetector {

        private final AdaptableInterval             maxCheck;
        private final int                           maxIter;
        private final BracketingNthOrderBrentSolver solver;
        private final int                           idx;
        private final int                           eventType;

        /** Constructor for the {@link Event} class.
         * @param maxCheck maximum checking interval, must be strictly positive (s)
         * @param threshold convergence threshold (s)
         * @param maxIter maximum number of iterations in the event time search
         * @param idx the index of the continuous variable to use
         * @param eventType the type of event to use. See {@link #g}
         */
        public Event(final double maxCheck, final double threshold, final int maxIter,
                     final int idx, final int eventType) {
            this.maxCheck  = s -> maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new BracketingNthOrderBrentSolver(0, threshold, 0, 5);
            this.idx       = idx;
            this.eventType = eventType;
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

        public ODEEventHandler getHandler() {
            return (state, detector, increasing) -> Action.STOP;
        }

        /** {@inheritDoc} */
        public double g(ODEStateAndDerivative s) {
            return (eventType == 0) ? s.getPrimaryState()[idx] >= 1.0 ? 1.0 : -1.0
                                    : s.getPrimaryState()[idx] - 1.0;
        }

    }

}
