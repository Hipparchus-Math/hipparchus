/*
 * Licensed to the Hipparchus project under one or more
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
package org.hipparchus.ode.events;

import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEIntegrator;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853FieldIntegrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.Decimal64;
import org.hipparchus.util.Decimal64Field;
import org.junit.Assert;
import org.junit.Test;

/**
 * Check events handlers and step handlers are called at consistent times.
 *
 * @author Luc Maisonobe
 */
public class EventsScheduling {

    @Test
    public void testForward() {
        doTest(0.0, 1.0, 32);
    }

    @Test
    public void testBackward() {
        doTest(1.0, 0.0, 32);
    }

    @Test
    public void testFieldForward() {
        doTestField(0.0, 1.0, 32);
    }

    @Test
    public void testFieldBackward() {
        doTestField(1.0, 0.0, 32);
    }

    private static void doTest(final double start, final double stop, final int expectedCalls) {

        final ODEIntegrator integrator =
                new DormandPrince853Integrator(10, 100.0, 1e-7, 1e-7);

        // checker that will be used in both step handler and events handlers
        // to check they are called in consistent order
        final ScheduleChecker checker = new ScheduleChecker(start, stop);
        integrator.addStepHandler((interpolator, islast) -> {
            checker.callTime(interpolator.getPreviousState().getTime());
            checker.callTime(interpolator.getCurrentState().getTime());
        });

        for (int i = 0; i < 10; ++i) {
            integrator.addEventHandler(new SimpleDetector(0.0625 * (i + 1), checker),
                                       1.0, 1.0e-9, 100);
        }

        final OrdinaryDifferentialEquation ode = new OrdinaryDifferentialEquation() {
            public int getDimension() {
                return 1;
            }
            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { 1 };
            }
        };

        final ODEState initialState = new ODEState(start, new double[] { 0.0 });

        integrator.integrate(new ExpandableODE(ode), initialState, stop);

        Assert.assertEquals(expectedCalls, checker.calls);

    }

    private static void doTestField(final double start, final double stop, final int expectedCalls) {

        final FieldODEIntegrator<Decimal64> integrator =
                new DormandPrince853FieldIntegrator<Decimal64>(Decimal64Field.getInstance(), 10, 100.0, 1e-7, 1e-7);

        // checker that will be used in both step handler and events handlers
        // to check they are called in consistent order
        final ScheduleChecker checker = new ScheduleChecker(start, stop);
        integrator.addStepHandler((interpolator, islast) -> {
            checker.callTime(interpolator.getPreviousState().getTime().getReal());
            checker.callTime(interpolator.getCurrentState().getTime().getReal());
        });

        for (int i = 0; i < 10; ++i) {
            integrator.addEventHandler(new SimpleFieldDetector(0.0625 * (i + 1), checker),
                                       1.0, 1.0e-9, 100);
        }

        final FieldOrdinaryDifferentialEquation<Decimal64> ode =
                        new FieldOrdinaryDifferentialEquation<Decimal64>() {
            public int getDimension() {
                return 1;
            }
            public Decimal64[] computeDerivatives(Decimal64 t, Decimal64[] y) {
                return new Decimal64[] { Decimal64.ONE };
            }
        };

        final FieldODEState<Decimal64> initialState =
                        new FieldODEState<>(new Decimal64(start), new Decimal64[] { Decimal64.ZERO });

        integrator.integrate(new FieldExpandableODE<>(ode), initialState, new Decimal64(stop));

        Assert.assertEquals(expectedCalls, checker.calls);

    }

    /** Checker for method calls scheduling. */
    private static class ScheduleChecker {

        private final double start;
        private final double stop;
        private double last;
        private int    calls;

        ScheduleChecker(final double start, final double stop) {
            this.start = start;
            this.stop  = stop;
            this.last  = Double.NaN;
            this.calls = 0;
        }

        void callTime(final double time) {
            if (!Double.isNaN(last)) {
                // check scheduling is always consistent with integration direction
                if (start < stop) {
                    // forward direction
                    Assert.assertTrue(time >= start);
                    Assert.assertTrue(time <= stop);
                    Assert.assertTrue(time >= last);
               } else {
                    // backward direction
                   Assert.assertTrue(time <= start);
                   Assert.assertTrue(time >= stop);
                   Assert.assertTrue(time <= last);
                }
            }
            last = time;
            ++calls;
        }

    }

    private static class SimpleDetector implements ODEEventHandler {

        private final double tEvent;
        private final ScheduleChecker checker;
        SimpleDetector(final double tEvent, final ScheduleChecker checker) {
            this.tEvent  = tEvent;
            this.checker = checker;
        }

        @Override
        public double g(final ODEStateAndDerivative state) {
            return state.getTime() - tEvent;
        }

        @Override
        public Action eventOccurred(final ODEStateAndDerivative state, final boolean increasing) {
            checker.callTime(state.getTime());
            return Action.CONTINUE;
        }
        
    }

    private static class SimpleFieldDetector implements FieldODEEventHandler<Decimal64> {

        private final double tEvent;
        private final ScheduleChecker checker;
        SimpleFieldDetector(final double tEvent, final ScheduleChecker checker) {
            this.tEvent  = tEvent;
            this.checker = checker;
        }

        @Override
        public Decimal64 g(final FieldODEStateAndDerivative<Decimal64> state) {
            return state.getTime().subtract(tEvent);
        }

        @Override
        public Action eventOccurred(final FieldODEStateAndDerivative<Decimal64> state, final boolean increasing) {
            checker.callTime(state.getTime().getReal());
            return Action.CONTINUE;
        }
        
    }

}
