/*
 * Licensed to the Hipparchus project under one or more
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
package org.hipparchus.ode.events;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketedUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
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
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
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
        integrator.addStepHandler((interpolator) -> {
            checker.callTime(interpolator.getPreviousState().getTime());
            checker.callTime(interpolator.getCurrentState().getTime());
        });

        for (int i = 0; i < 10; ++i) {
            integrator.addEventDetector(new SimpleDetector(0.0625 * (i + 1), checker,
                                                           1.0, 1.0e-9, 100));
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

        final FieldODEIntegrator<Binary64> integrator =
                new DormandPrince853FieldIntegrator<Binary64>(Binary64Field.getInstance(), 10, 100.0, 1e-7, 1e-7);

        // checker that will be used in both step handler and events handlers
        // to check they are called in consistent order
        final ScheduleChecker checker = new ScheduleChecker(start, stop);
        integrator.addStepHandler((interpolator) -> {
            checker.callTime(interpolator.getPreviousState().getTime().getReal());
            checker.callTime(interpolator.getCurrentState().getTime().getReal());
        });

        for (int i = 0; i < 10; ++i) {
            integrator.addEventDetector(new SimpleFieldDetector(0.0625 * (i + 1), checker, 1.0, 1.0e-9, 100));
        }

        final FieldOrdinaryDifferentialEquation<Binary64> ode =
                        new FieldOrdinaryDifferentialEquation<Binary64>() {
            public int getDimension() {
                return 1;
            }
            public Binary64[] computeDerivatives(Binary64 t, Binary64[] y) {
                return new Binary64[] { Binary64.ONE };
            }
        };

        final FieldODEState<Binary64> initialState =
                        new FieldODEState<>(new Binary64(start), new Binary64[] { Binary64.ZERO });

        integrator.integrate(new FieldExpandableODE<>(ode), initialState, new Binary64(stop));

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

    private static class SimpleDetector implements ODEEventDetector {

        private final AdaptableInterval             maxCheck;
        private final int                           maxIter;
        private final BracketingNthOrderBrentSolver solver;
        private final double                        tEvent;
        private final ScheduleChecker               checker;
        SimpleDetector(final double tEvent, final ScheduleChecker checker,
                       final double maxCheck, final double threshold, final int maxIter) {
            this.maxCheck  = s -> maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new BracketingNthOrderBrentSolver(0, threshold, 0, 5);
            this.tEvent    = tEvent;
            this.checker   = checker;
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
            return (state, detector, increasing) -> {
                checker.callTime(state.getTime());
                return Action.CONTINUE;
            };
        }

        @Override
        public double g(final ODEStateAndDerivative state) {
            return state.getTime() - tEvent;
        }

    }

    private static class SimpleFieldDetector implements FieldODEEventDetector<Binary64> {

        private final FieldAdaptableInterval<Binary64>             maxCheck;
        private final int                                          maxIter;
        private final BracketedRealFieldUnivariateSolver<Binary64> solver;
        private final double                                       tEvent;
        private final ScheduleChecker                              checker;

        SimpleFieldDetector(final double tEvent, final ScheduleChecker checker,
                            final double maxCheck, final double threshold, final int maxIter) {
            this.maxCheck  = s -> maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new FieldBracketingNthOrderBrentSolver<>(new Binary64(0),
                                                                      new Binary64(threshold),
                                                                      new Binary64(0),
                                                                      5);
            this.tEvent    = tEvent;
            this.checker   = checker;
        }

        public FieldAdaptableInterval<Binary64> getMaxCheckInterval() {
            return maxCheck;
        }

        public int getMaxIterationCount() {
            return maxIter;
        }

        public BracketedRealFieldUnivariateSolver<Binary64> getSolver() {
            return solver;
        }

        public FieldODEEventHandler<Binary64> getHandler() {
            return (state, detector, increasing) -> {
                checker.callTime(state.getTime().getReal());
                return Action.CONTINUE;
            };
        }
        
        @Override
        public Binary64 g(final FieldODEStateAndDerivative<Binary64> state) {
            return state.getTime().subtract(tEvent);
        }

    }

}
