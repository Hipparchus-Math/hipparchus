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
import org.hipparchus.ode.EquationsMapper;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.SecondaryODE;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.ode.nonstiff.LutherIntegrator;
import org.hipparchus.ode.sampling.DummyStepInterpolator;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.junit.Assert;
import org.junit.Test;

public class DetectorBasedEventStateTest {

    // JIRA: MATH-322
    @Test
    public void closeEvents() throws MathIllegalArgumentException, MathIllegalStateException {

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
        Assert.assertFalse(es.evaluateStep(interpolatorA));

        ODEStateAndDerivative osdBetween    = new ODEStateAndDerivative(0.5 * (r1 + r2), a, a);
        ODEStateInterpolator interpolatorB  = new DummyStepInterpolator(true,
                                                                        osBefore, osdBetween,
                                                                        osBefore, osdBetween,
                                                                        mapper);
        Assert.assertTrue(es.evaluateStep(interpolatorB));
        Assert.assertEquals(r1, es.getEventTime(), tolerance);
        ODEStateAndDerivative osdAtEvent    = new ODEStateAndDerivative(es.getEventTime(), a, a);
        es.doEvent(osdAtEvent);

        ODEStateAndDerivative osdAfterSecond = new ODEStateAndDerivative(r2 + 0.4 * gap, a, a);
        ODEStateInterpolator interpolatorC  = new DummyStepInterpolator(true,
                                                                        osdAtEvent, osdAfterSecond,
                                                                        osdAtEvent, osdAfterSecond,
                                                                        mapper);
        Assert.assertTrue(es.evaluateStep(interpolatorC));
        Assert.assertEquals(r2, es.getEventTime(), tolerance);

    }

    // Jira: MATH-695
    @Test
    public void testIssue695()
        throws MathIllegalArgumentException, MathIllegalStateException {

        OrdinaryDifferentialEquation equation = new OrdinaryDifferentialEquation() {
            public int getDimension() {
                return 1;
            }
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
        Assert.assertEquals(target, finalState.getTime(), 1.0e-10);
        Assert.assertEquals(32.0, finalState.getPrimaryState()[0], 1.0e-10);

    }

    private static class ResettingEvent implements ODEEventDetector {

        private static double lastTriggerTime = Double.NEGATIVE_INFINITY;
        private final AdaptableInterval             maxCheck;
        private final int                           maxIter;
        private final BracketingNthOrderBrentSolver solver;
        private final double                        tEvent;

        public ResettingEvent(final double tEvent,
                              final double maxCheck, final double threshold, final int maxIter) {
            this.maxCheck  = s -> maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new BracketingNthOrderBrentSolver(0, threshold, 0, 5);
            this.tEvent    = tEvent;
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

        public double g(ODEStateAndDerivative s) {
            // the bug corresponding to issue 695 causes the g function
            // to be called at obsolete times t despite an event
            // occurring later has already been triggered.
            // When this occurs, the following assertion is violated
            Assert.assertTrue("going backard in time! (" + s.getTime() + " < " + lastTriggerTime + ")",
                              s.getTime() >= lastTriggerTime);
            return s.getTime() - tEvent;
        }

        public ODEEventHandler getHandler() {
            return new ODEEventHandler() {
                public Action eventOccurred(ODEStateAndDerivative s, ODEEventDetector detector, boolean increasing) {
                    // remember in a class variable when the event was triggered
                    lastTriggerTime = s.getTime();
                    return Action.RESET_STATE;
                }

                public ODEStateAndDerivative resetState(ODEEventDetector detector, ODEStateAndDerivative s) {
                    double[] y = s.getPrimaryState();
                    y[0] += 1.0;
                    return new ODEStateAndDerivative(s.getTime(), y, s.getPrimaryDerivative());
                }
            };
        }

    }

    // Jira: MATH-965
    @Test
    public void testIssue965()
        throws MathIllegalArgumentException, MathIllegalStateException {

        ExpandableODE equation = new ExpandableODE(new OrdinaryDifferentialEquation() {
            public int getDimension() {
                return 1;
            }
            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { 2.0 };
            }
        });
        int index = equation.addSecondaryEquations(new SecondaryODE() {
            public int getDimension() {
                return 1;
            }
            public double[] computeDerivatives(double t, double[] primary,
                                           double[] primaryDot, double[] secondary) {
                return new double[] { -3.0 };
            }
        });
        Assert.assertEquals(1, index);

        DormandPrince853Integrator integrator = new DormandPrince853Integrator(0.001, 1000, 1.0e-14, 1.0e-14);
        integrator.addEventDetector(new SecondaryStateEvent(index, -3.0, 0.1, 1.0e-9, 1000));
        integrator.setInitialStepSize(3.0);

        ODEState initialState = new ODEState(0.0,
                                             new double[] { 0.0 },
                                             new double[][] { { 0.0 } });
        ODEStateAndDerivative finalState = integrator.integrate(equation, initialState, 30.0);
        Assert.assertEquals( 1.0, finalState.getTime(), 1.0e-10);
        Assert.assertEquals( 2.0, finalState.getPrimaryState()[0], 1.0e-10);
        Assert.assertEquals(-3.0, finalState.getSecondaryState(index)[0], 1.0e-10);

    }

    private static class SecondaryStateEvent implements ODEEventDetector {

        private final AdaptableInterval             maxCheck;
        private final int                           maxIter;
        private final BracketingNthOrderBrentSolver solver;
        private int                                 index;
        private final double                        target;

        public SecondaryStateEvent(final int index, final double target,
                                   final double maxCheck, final double threshold, final int maxIter) {
            this.maxCheck  = s -> maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new BracketingNthOrderBrentSolver(0, threshold, 0, 5);
            this.index     = index;
            this.target    = target;
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

        /** {@inheritDoc} */
        public ODEEventHandler getHandler() {
            return (state, detector, increasing) -> Action.STOP;
        }

        public double g(ODEStateAndDerivative s) {
            return s.getSecondaryState(index)[0] - target;
        }

    }

    @Test
    public void testEventsCloserThanThreshold()
        throws MathIllegalArgumentException, MathIllegalStateException {

        OrdinaryDifferentialEquation equation = new OrdinaryDifferentialEquation() {

            public int getDimension() {
                return 1;
            }

            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { 1.0 };
            }
        };

        LutherIntegrator integrator = new LutherIntegrator(20.0);
        CloseEventsGenerator eventsGenerator =
                        new CloseEventsGenerator(9.0 - 1.0 / 128, 9.0 + 1.0 / 128, 1.0, 0.02, 1000);
        integrator.addEventDetector(eventsGenerator);
        double tEnd = integrator.integrate(equation, new ODEState(0.0, new double[1]), 100.0).getTime();
        Assert.assertEquals( 2, eventsGenerator.getCount());
        Assert.assertEquals( 9.0 + 1.0 / 128, tEnd, 1.0 / 32.0);

    }

    private class CloseEventsGenerator implements ODEEventDetector {

        private final AdaptableInterval             maxCheck;
        private final int                           maxIter;
        private final BracketingNthOrderBrentSolver solver;
        final double                                r1;
        final double                                r2;
        int                                         count;

        public CloseEventsGenerator(final double r1, final double r2,
                                    final double maxCheck, final double threshold, final int maxIter) {
            this.maxCheck  = s -> maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new BracketingNthOrderBrentSolver(0, threshold, 0, 5);
            this.r1        = r1;
            this.r2        = r2;
            this.count     = 0;
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

        public double g(ODEStateAndDerivative s) {
            return (s.getTime() - r1) * (r2 - s.getTime());
        }

        public ODEEventHandler getHandler() {
            return (state, detector, increasing) -> ++count < 2 ? Action.CONTINUE : Action.STOP;
        }

        public int getCount() {
            return count;
        }

    }

}
