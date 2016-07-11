/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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

package org.hipparchus.ode.nonstiff;


import java.util.Arrays;
import java.util.List;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.LocalizedODEFormats;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEJacobiansProvider;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.SecondaryODE;
import org.hipparchus.ode.TestProblem1;
import org.hipparchus.ode.TestProblem3;
import org.hipparchus.ode.TestProblem4;
import org.hipparchus.ode.TestProblem5;
import org.hipparchus.ode.TestProblemHandler;
import org.hipparchus.ode.VariationalEquation;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.ODEEventHandler;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public abstract class EmbeddedRungeKuttaIntegratorAbstractTest {

    protected abstract EmbeddedRungeKuttaIntegrator
    createIntegrator(final double minStep, final double maxStep,
                     final double scalAbsoluteTolerance, final double scalRelativeTolerance);

    protected abstract EmbeddedRungeKuttaIntegrator
    createIntegrator(final double minStep, final double maxStep,
                     final double[] vecAbsoluteTolerance, final double[] vecRelativeTolerance);

    @Test
    public abstract void testForwardBackwardExceptions();

    protected void doTestForwardBackwardExceptions() {
        OrdinaryDifferentialEquation equations = new OrdinaryDifferentialEquation() {

            public int getDimension() {
                return 1;
            }

            public double[] computeDerivatives(double t, double[] y) {
                if (t < -0.5) {
                    throw new LocalException();
                } else {
                    throw new RuntimeException("oops");
                }
            }
        };

        EmbeddedRungeKuttaIntegrator integrator = createIntegrator(0.0, 1.0, 1.0e-10, 1.0e-10);

        try  {
            integrator.integrate(new ExpandableODE(equations),
                                 new ODEState(-1, new double[1]),
                                 0);
            Assert.fail("an exception should have been thrown");
          } catch(LocalException de) {
            // expected behavior
          }

          try  {
              integrator.integrate(new ExpandableODE(equations),
                                   new ODEState(0, new double[1]),
                                   1);
               Assert.fail("an exception should have been thrown");
          } catch(RuntimeException de) {
            // expected behavior
          }
    }

    protected static class LocalException extends RuntimeException {
        private static final long serialVersionUID = 20151208L;
    }

    @Test
    public void testMinStep() {
        try {
            TestProblem1 pb = new TestProblem1();
            double minStep = 0.1 * (pb.getFinalTime() - pb.getInitialState().getTime());
            double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
            double[] vecAbsoluteTolerance = { 1.0e-15, 1.0e-16 };
            double[] vecRelativeTolerance = { 1.0e-15, 1.0e-16 };

            ODEIntegrator integ = createIntegrator(minStep, maxStep,
                                                   vecAbsoluteTolerance, vecRelativeTolerance);
            TestProblemHandler handler = new TestProblemHandler(pb, integ);
            integ.addStepHandler(handler);
            integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedODEFormats.MINIMAL_STEPSIZE_REACHED_DURING_INTEGRATION,
                                miae.getSpecifier());
        }

    }

    @Test
    public abstract void testIncreasingTolerance();

    protected void doTestIncreasingTolerance(double factor, double epsilon) {

        int previousCalls = Integer.MAX_VALUE;
        for (int i = -12; i < -2; ++i) {
            TestProblem1 pb = new TestProblem1();
            double minStep = 0;
            double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
            double scalAbsoluteTolerance = FastMath.pow(10.0, i);
            double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

            ODEIntegrator integ = createIntegrator(minStep, maxStep,
                                                                  scalAbsoluteTolerance, scalRelativeTolerance);
            TestProblemHandler handler = new TestProblemHandler(pb, integ);
            integ.addStepHandler(handler);
            integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

            Assert.assertTrue(handler.getMaximalValueError() < (factor * scalAbsoluteTolerance));
            Assert.assertEquals(0, handler.getMaximalTimeError(), epsilon);

            int calls = pb.getCalls();
            Assert.assertEquals(integ.getEvaluations(), calls);
            Assert.assertTrue(calls <= previousCalls);
            previousCalls = calls;

        }

    }

    @Test
    public abstract void testEvents();

    protected void doTestEvents(final double epsilonMaxValue, final String name) {

      TestProblem4 pb = new TestProblem4();
      double minStep = 0;
      double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
      double scalAbsoluteTolerance = 1.0e-8;
      double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

      ODEIntegrator integ = createIntegrator(minStep, maxStep,
                                                            scalAbsoluteTolerance, scalRelativeTolerance);
      TestProblemHandler handler = new TestProblemHandler(pb, integ);
      integ.addStepHandler(handler);
      ODEEventHandler[] functions = pb.getEventsHandlers();
      double convergence = 1.0e-8 * maxStep;
      for (int l = 0; l < functions.length; ++l) {
          integ.addEventHandler(functions[l], Double.POSITIVE_INFINITY, convergence, 1000);
      }
      Assert.assertEquals(functions.length, integ.getEventHandlers().size());
      integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

      Assert.assertEquals(0, handler.getMaximalValueError(), epsilonMaxValue);
      Assert.assertEquals(0, handler.getMaximalTimeError(), convergence);
      Assert.assertEquals(12.0, handler.getLastTime(), convergence);
      Assert.assertEquals(name, integ.getName());
      integ.clearEventHandlers();
      Assert.assertEquals(0, integ.getEventHandlers().size());

    }

    @Test
    public void testEventsErrors() {
        try {
            final TestProblem1 pb = new TestProblem1();
            double minStep = 0;
            double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
            double scalAbsoluteTolerance = 1.0e-8;
            double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

            ODEIntegrator integ = createIntegrator(minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
            TestProblemHandler handler = new TestProblemHandler(pb, integ);
            integ.addStepHandler(handler);

            integ.addEventHandler(new ODEEventHandler() {
                public void init(ODEStateAndDerivative state0, double t) {
                }
                public Action eventOccurred(ODEStateAndDerivative state, boolean increasing) {
                    return Action.CONTINUE;
                }
                public double g(ODEStateAndDerivative state) {
                    double middle = 0.5 * (pb.getInitialState().getTime() + pb.getFinalTime());
                    double offset = state.getTime() - middle;
                    if (offset > 0) {
                        throw new LocalException();
                    }
                    return offset;
                }
                public ODEState resetState(ODEStateAndDerivative state) {
                    return state;
                }
            }, Double.POSITIVE_INFINITY, 1.0e-8 * maxStep, 1000);

            integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());
        } catch (LocalException le) {
            // expected
        }

    }

    @Test
    public void testEventsNoConvergence() {

        final TestProblem1 pb = new TestProblem1();
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

        ODEIntegrator integ = createIntegrator(minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);

        integ.addEventHandler(new ODEEventHandler() {
            public Action eventOccurred(ODEStateAndDerivative state, boolean increasing) {
                return Action.CONTINUE;
            }
            public double g(ODEStateAndDerivative state) {
                double middle = 0.5 * (pb.getInitialState().getTime() + pb.getFinalTime());
                double offset = state.getTime() - middle;
                return (offset > 0) ? offset + 0.5 : offset - 0.5;
            }
            public ODEState resetState(ODEStateAndDerivative state) {
                return state;
            }
        }, Double.POSITIVE_INFINITY, 1.0e-8 * maxStep, 3);

        try {
            integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mcee) {
            // Expected.
        }

    }

    @Test
    public void testSanityChecks() {
        TestProblem3 pb = new TestProblem3();
        try  {
            EmbeddedRungeKuttaIntegrator integrator = createIntegrator(0,
                                                                       pb.getFinalTime() - pb.getInitialState().getTime(),
                                                                       new double[4], new double[4]);
            integrator.integrate(new ExpandableODE(pb),
                                 new ODEState(pb.getInitialState().getTime(), new double[6]),
                                 pb.getFinalTime());
            Assert.fail("an exception should have been thrown");
        } catch(MathIllegalArgumentException ie) {
        }
        try  {
            EmbeddedRungeKuttaIntegrator integrator =
                            createIntegrator(0,
                                             pb.getFinalTime() - pb.getInitialState().getTime(),
                                             new double[2], new double[4]);
            integrator.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());
            Assert.fail("an exception should have been thrown");
        } catch(MathIllegalArgumentException ie) {
        }
        try  {
            EmbeddedRungeKuttaIntegrator integrator =
                            createIntegrator(0,
                                             pb.getFinalTime() - pb.getInitialState().getTime(),
                                             new double[4], new double[4]);
            integrator.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getInitialState().getTime());
            Assert.fail("an exception should have been thrown");
        } catch(MathIllegalArgumentException ie) {
        }
    }

    @Test
    public void testNullIntervalCheck()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        try {
            TestProblem1 pb = new TestProblem1();
            EmbeddedRungeKuttaIntegrator integrator = createIntegrator(0.0, 1.0, 1.0e-10, 1.0e-10);
            integrator.integrate(pb,
                                 new ODEState(0.0, new double[pb.getDimension()]),
                                 0.0);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedODEFormats.TOO_SMALL_INTEGRATION_INTERVAL,
                                miae.getSpecifier());
        }
    }

    @Test
    public abstract void testBackward();

    protected void doTestBackward(final double epsilonLast, final double epsilonMaxValue,
                                  final double epsilonMaxTime, final String name)
        throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem5 pb = new TestProblem5();
        double minStep = 0;
        double maxStep = FastMath.abs(pb.getFinalTime() - pb.getInitialState().getTime());
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

        EmbeddedRungeKuttaIntegrator integ = createIntegrator(minStep, maxStep,
                                                                      scalAbsoluteTolerance,
                                                                      scalRelativeTolerance);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

        Assert.assertEquals(0, handler.getLastError(),         epsilonLast);
        Assert.assertEquals(0, handler.getMaximalValueError(), epsilonMaxValue);
        Assert.assertEquals(0, handler.getMaximalTimeError(),  epsilonMaxTime);
        Assert.assertEquals(name, integ.getName());

    }

    @Test
    public abstract void testKepler();

    protected void doTestKepler(double epsilon) {

        final TestProblem3 pb  = new TestProblem3(0.9);
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
        double[] vecAbsoluteTolerance = { 1.0e-8, 1.0e-8, 1.0e-10, 1.0e-10 };
        double[] vecRelativeTolerance = { 1.0e-10, 1.0e-10, 1.0e-8, 1.0e-8 };

        ODEIntegrator integ = createIntegrator(minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
        integ.addStepHandler(new KeplerHandler(pb, epsilon));
        integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());
    }

    private static class KeplerHandler implements ODEStepHandler {
        private double maxError;
        private final TestProblem3 pb;
        private final double epsilon;
        public KeplerHandler(TestProblem3 pb, double epsilon) {
            this.pb      = pb;
            this.epsilon = epsilon;
            maxError     = 0;
        }
        public void init(ODEStateAndDerivative state0, double t) {
            maxError = 0;
        }
        public void handleStep(ODEStateInterpolator interpolator, boolean isLast) {

            ODEStateAndDerivative current = interpolator.getCurrentState();
            double[] theoreticalY  = pb.computeTheoreticalState(current.getTime());
            double dx = current.getPrimaryState()[0] - theoreticalY[0];
            double dy = current.getPrimaryState()[1] - theoreticalY[1];
            double error = dx * dx + dy * dy;
            if (error > maxError) {
                maxError = error;
            }
            if (isLast) {
                Assert.assertEquals(0.0, maxError, epsilon);
            }
        }
    }

    @Test
    public abstract void testMissedEndEvent();

    protected void doTestMissedEndEvent(final double epsilonY, final double epsilonT) {
        final double   t0     = 1878250320.0000029;
        final double   tEvent = 1878250379.9999986;
        final double[] k  = { 1.0e-4, 1.0e-5, 1.0e-6 };
        OrdinaryDifferentialEquation ode = new OrdinaryDifferentialEquation() {

            public int getDimension() {
                return k.length;
            }

            public double[] computeDerivatives(double t, double[] y) {
                final double[] yDot = new double[y.length];
                for (int i = 0; i < y.length; ++i) {
                    yDot[i] = k[i] * y[i];
                }
                return yDot;
            }
        };

        EmbeddedRungeKuttaIntegrator integrator = createIntegrator(0.0, 100.0, 1.0e-10, 1.0e-10);

        double[] y0   = new double[k.length];
        for (int i = 0; i < y0.length; ++i) {
            y0[i] = i + 1;
        }

        integrator.setInitialStepSize(60.0);
        ODEStateAndDerivative finalState = integrator.integrate(new ExpandableODE(ode), new ODEState(t0, y0), tEvent);
        Assert.assertEquals(tEvent, finalState.getTime(), epsilonT);
        double[] y = finalState.getPrimaryState();
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(y0[i] * FastMath.exp(k[i] * (finalState.getTime() - t0)), y[i], epsilonY);
        }

        integrator.setInitialStepSize(60.0);
        integrator.addEventHandler(new ODEEventHandler() {

            public double g(ODEStateAndDerivative s) {
                return s.getTime() - tEvent;
            }

            public Action eventOccurred(ODEStateAndDerivative s, boolean increasing) {
                Assert.assertEquals(tEvent, s.getTime(), epsilonT);
                return Action.CONTINUE;
            }
        }, Double.POSITIVE_INFINITY, 1.0e-20, 100);
        finalState = integrator.integrate(new ExpandableODE(ode), new ODEState(t0, y0), tEvent + 120);
        Assert.assertEquals(tEvent + 120, finalState.getTime(), epsilonT);
        y = finalState.getPrimaryState();
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(y0[i] * FastMath.exp(k[i] * (finalState.getTime() - t0)), y[i], epsilonY);
        }

    }

    @Test
    public void testTooLargeFirstStep() {

        AdaptiveStepsizeIntegrator integ =
                        createIntegrator(0, Double.POSITIVE_INFINITY, Double.NaN, Double.NaN);
        final double start = 0.0;
        final double end   = 0.001;
        OrdinaryDifferentialEquation equations = new OrdinaryDifferentialEquation() {

            public int getDimension() {
                return 1;
            }

            public double[] computeDerivatives(double t, double[] y) {
                Assert.assertTrue(t >= FastMath.nextAfter(start, Double.NEGATIVE_INFINITY));
                Assert.assertTrue(t <= FastMath.nextAfter(end,   Double.POSITIVE_INFINITY));
                return new double[] { -100.0 * y[0] };
            }

        };

        integ.setStepSizeControl(0, 1.0, 1.0e-6, 1.0e-8);
        integ.integrate(equations, new ODEState(start, new double[] { 1.0 }), end);

    }

    @Test
    public abstract void testVariableSteps();

    protected void doTestVariableSteps(final double min, final double max) {

        final TestProblem3 pb  = new TestProblem3(0.9);
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = scalAbsoluteTolerance;

        ODEIntegrator integ = createIntegrator(minStep, maxStep,
                                               scalAbsoluteTolerance,
                                               scalRelativeTolerance);
        integ.addStepHandler(new VariableHandler(min, max));
        double stopTime = integ.integrate(pb, pb.getInitialState(), pb.getFinalTime()).getTime();
        Assert.assertEquals(pb.getFinalTime(), stopTime, 1.0e-10);
    }

    @Test
    public abstract void testUnstableDerivative();

    protected void doTestUnstableDerivative(final double epsilon) {
        final StepProblem stepProblem = new StepProblem(0.0, 1.0, 2.0);
        ODEIntegrator integ = createIntegrator(0.1, 10, 1.0e-12, 0.0);
        integ.addEventHandler(stepProblem, 1.0, 1.0e-12, 1000);
        final ODEStateAndDerivative finalState =
                        integ.integrate(stepProblem, new ODEState(0.0, new double[] { 0.0 }), 10.0);
        Assert.assertEquals(8.0, finalState.getPrimaryState()[0], epsilon);
    }

    @Test
    public void testEventsScheduling() {

        OrdinaryDifferentialEquation sincos = new OrdinaryDifferentialEquation() {

            public int getDimension() {
                return 2;
            }

            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { y[1], -y[0] };
            }

        };

        SchedulingChecker sinChecker = new SchedulingChecker(0); // events at 0, PI, 2PI ...
        SchedulingChecker cosChecker = new SchedulingChecker(1); // events at PI/2, 3PI/2, 5PI/2 ...

        ODEIntegrator integ = createIntegrator(0.001, 1.0, 1.0e-12, 0.0);
        integ.addEventHandler(sinChecker, 0.01, 1.0e-7, 100);
        integ.addStepHandler(sinChecker);
        integ.addEventHandler(cosChecker, 0.01, 1.0e-7, 100);
        integ.addStepHandler(cosChecker);
        double   t0 = 0.5;
        double[] y0 = new double[] { FastMath.sin(t0), FastMath.cos(t0) };
        double   t  = 10.0;
        integ.integrate(sincos, new ODEState(t0, y0), t);

    }

    private static class SchedulingChecker implements ODEStepHandler, ODEEventHandler {

        int index;
        double tMin;

        public SchedulingChecker(int index) {
            this.index = index;
        }

        public void init(ODEStateAndDerivative s0, double t) {
            tMin = s0.getTime();
        }

        public void handleStep(ODEStateInterpolator interpolator, boolean isLast) {
            tMin = interpolator.getCurrentState().getTime();
        }

        public double g(ODEStateAndDerivative s) {
            // once a step has been handled by handleStep,
            // events checking should only refer to dates after the step
            Assert.assertTrue(s.getTime() >= tMin);
            return s.getPrimaryState()[index];
        }

        public Action eventOccurred(ODEStateAndDerivative s, boolean increasing) {
            return Action.RESET_STATE;
        }

        public ODEStateAndDerivative resetState(ODEStateAndDerivative s) {
            // in fact, we don't need to reset anything for the test
            return s;
        }

    }

    private static class VariableHandler implements ODEStepHandler {
        final double min;
        final double max;
        public VariableHandler(final double min, final double max) {
            this.min  = min;
            this.max  = max;
            firstTime = true;
            minStep   = 0;
            maxStep   = 0;
        }
        public void init(ODEStateAndDerivative s0, double t) {
            firstTime = true;
            minStep = 0;
            maxStep = 0;
        }
        public void handleStep(ODEStateInterpolator interpolator,
                               boolean isLast) {

            double step = FastMath.abs(interpolator.getCurrentState().getTime() -
                                       interpolator.getPreviousState().getTime());
            if (firstTime) {
                minStep   = FastMath.abs(step);
                maxStep   = minStep;
                firstTime = false;
            } else {
                if (step < minStep) {
                    minStep = step;
                }
                if (step > maxStep) {
                    maxStep = step;
                }
            }

            if (isLast) {
                Assert.assertEquals(min, minStep, 0.01 * min);
                Assert.assertEquals(max, maxStep, 0.01 * max);
            }
        }
        private boolean firstTime = true;
        private double  minStep = 0;
        private double  maxStep = 0;
    }

    @Test
    public void testWrongDerivative() {
        EmbeddedRungeKuttaIntegrator integrator =
                        createIntegrator(0.0, 1.0, 1.0e-10, 1.0e-10);
        OrdinaryDifferentialEquation equations =
                        new OrdinaryDifferentialEquation() {
            public double[] computeDerivatives(double t, double[] y) {
                if (t < -0.5) {
                    throw new LocalException();
                } else {
                    throw new RuntimeException("oops");
                }
            }
            public int getDimension() {
                return 1;
            }
        };

        try  {
            integrator.integrate(equations, new ODEState(-1.0, new double[1]), 0.0);
            Assert.fail("an exception should have been thrown");
        } catch(LocalException de) {
            // expected behavior
        }

        try  {
            integrator.integrate(equations, new ODEState(0.0, new double[1]), 1.0);
            Assert.fail("an exception should have been thrown");
        } catch(RuntimeException de) {
            // expected behavior
        }

    }

    @Test
    public abstract void testPartialDerivatives();

    protected void doTestPartialDerivatives(final double epsilonY,
                                            final double epsilonPartials) {

        double omega = 1.3;
        double t0    = 1.3;
        double[] y0  = new double[] { 3.0, 4.0 };
        double t     = 6.0;
        SinCos sinCos = new SinCos(omega);

        ExpandableODE       expandable   = new ExpandableODE(sinCos);
        VariationalEquation ve           = new VariationalEquation(expandable, sinCos);
        ODEState            initialState = ve.setUpInitialState(new ODEState(t0, y0));

        EmbeddedRungeKuttaIntegrator integrator =
                        createIntegrator(0.001 * (t - t0), t - t0, 1.0e-12, 1.0e-12);
        ODEStateAndDerivative finalState = integrator.integrate(expandable, initialState, t);

        // check that the final state contains both the state vector and its partial derivatives
        int n = sinCos.getDimension();
        int p = sinCos.getParametersNames().size();
        Assert.assertEquals(n,               finalState.getPrimaryStateDimension());
        Assert.assertEquals(n + n * (n + p), finalState.getCompleteStateDimension());

        // check values
        for (int i = 0; i < sinCos.getDimension(); ++i) {
            Assert.assertEquals(sinCos.theoreticalY(t)[i], finalState.getPrimaryState()[i], epsilonY);
        }

        // check derivatives
        double[][] dydy0 = ve.extractMainSetJacobian(finalState);
        for (int i = 0; i < dydy0.length; ++i) {
            for (int j = 0; j < dydy0[i].length; ++j) {
                Assert.assertEquals(sinCos.exactDyDy0(t)[i][j], dydy0[i][j], epsilonPartials);
            }
        }

        double[] dydp0 = ve.extractParameterJacobian(finalState, SinCos.OMEGA_PARAMETER);
        for (int i = 0; i < dydp0.length; ++i) {
            Assert.assertEquals(sinCos.exactDyDomega(t)[i], dydp0[i], epsilonPartials);
        }

    }

    @Test
    public abstract void testSecondaryEquations();

    protected void doTestSecondaryEquations(final double epsilonSinCos,
                                            final double epsilonLinear) {
        OrdinaryDifferentialEquation sinCos = new OrdinaryDifferentialEquation() {

            @Override
            public int getDimension() {
                return 2;
            }

            @Override
            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { y[1], -y[0] };
            }

        };

        SecondaryODE linear = new SecondaryODE() {

            @Override
            public int getDimension() {
                return 1;
            }

            @Override
            public double[] computeDerivatives(double t, double[] primary, double[] primaryDot, double[] secondary) {
                return new double[] { -1 };
            }

        };

        ExpandableODE expandable = new ExpandableODE(sinCos);
        expandable.addSecondaryEquations(linear);

        ODEIntegrator integrator = createIntegrator(0.001, 1.0, 1.0e-12, 1.0e-12);
        final double[] max = new double[2];
        integrator.addStepHandler(new ODEStepHandler() {
            @Override
            public void handleStep(ODEStateInterpolator interpolator, boolean isLast) {
                for (int i = 0; i <= 10; ++i) {
                    double tPrev = interpolator.getPreviousState().getTime();
                    double tCurr = interpolator.getCurrentState().getTime();
                    double t     = (tPrev * (10 - i) + tCurr * i) / 10;
                    ODEStateAndDerivative state = interpolator.getInterpolatedState(t);
                    Assert.assertEquals(2, state.getPrimaryStateDimension());
                    Assert.assertEquals(1, state.getNumberOfSecondaryStates());
                    Assert.assertEquals(2, state.getSecondaryStateDimension(0));
                    Assert.assertEquals(1, state.getSecondaryStateDimension(1));
                    Assert.assertEquals(3, state.getCompleteStateDimension());
                    max[0] = FastMath.max(max[0],
                                          FastMath.abs(FastMath.sin(t) - state.getPrimaryState()[0]));
                    max[0] = FastMath.max(max[0],
                                          FastMath.abs(FastMath.cos(t) - state.getPrimaryState()[1]));
                    max[1] = FastMath.max(max[1],
                                          FastMath.abs(1 - t - state.getSecondaryState(1)[0]));
                }
            }
        });

        double[] primary0 = new double[] { 0.0, 1.0 };
        double[][] secondary0 =  new double[][] { { 1.0 } };
        ODEState initialState = new ODEState(0.0, primary0, secondary0);

        ODEStateAndDerivative finalState =
                        integrator.integrate(expandable, initialState, 10.0);
        Assert.assertEquals(10.0, finalState.getTime(), 1.0e-12);
        Assert.assertEquals(0, max[0], epsilonSinCos);
        Assert.assertEquals(0, max[1], epsilonLinear);

    }

    @Test
    public void testNaNAppearing() {
        try {
            ODEIntegrator integ = createIntegrator(0.01, 1.0, 0.1, 0.1);
            integ.integrate(new OrdinaryDifferentialEquation() {
                public int getDimension() {
                    return 1;
                }
                public double[] computeDerivatives(double t, double[] y) {
                    return new double[] { FastMath.log(t) };
                }
            }, new ODEState(1.0, new double[] { 1.0 }), -1.0);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedODEFormats.NAN_APPEARING_DURING_INTEGRATION, mise.getSpecifier());
            Assert.assertTrue(((Double) mise.getParts()[0]).doubleValue() <= 0.0);
        }
    }

    private static class SinCos implements ODEJacobiansProvider {

        public static String OMEGA_PARAMETER = "omega";

        private final double omega;
        private       double r;
        private       double alpha;

        private double dRdY00;
        private double dRdY01;
        private double dAlphadOmega;
        private double dAlphadY00;
        private double dAlphadY01;

        protected SinCos(final double omega) {
            this.omega = omega;
        }

        public int getDimension() {
            return 2;
        }

        public void init(final double t0, final double[] y0,
                         final double finalTime) {

            // theoretical solution is y(t) = { r * sin(omega * t + alpha), r * cos(omega * t + alpha) }
            // so we retrieve alpha by identification from the initial state
            final double r2 = y0[0] * y0[0] + y0[1] * y0[1];

            this.r            = FastMath.sqrt(r2);
            this.dRdY00       = y0[0] / r;
            this.dRdY01       = y0[1] / r;

            this.alpha        = FastMath.atan2(y0[0], y0[1]) - t0 * omega;
            this.dAlphadOmega = -t0;
            this.dAlphadY00   =  y0[1] / r2;
            this.dAlphadY01   = -y0[0] / r2;

        }

        @Override
        public double[] computeDerivatives(final double t, final double[] y) {
            return new double[] {
                omega *  y[1],
                omega * -y[0]
            };
        }

        @Override
        public double[][] computeMainStateJacobian(final double t, final double[] y, final double[] yDot) {
            // this is the Jacobian of dYdot/dY
            return new double[][] {
                {   0.0,  omega },
                { -omega,  0.0  }
            };
        }

        @Override
        public double[] computeParameterJacobian(final double t, final double[] y, final double[] yDot,
                                                 final String paramName)
            throws MathIllegalArgumentException, MathIllegalStateException {
            if (!isSupported(paramName)) {
                throw new MathIllegalArgumentException(LocalizedODEFormats.UNKNOWN_PARAMETER, paramName);
            }
            // this is the Jacobian of dYdot/dOmega
            return new double[] {
                y[1],
               -y[0]
            };
        }

        public double[] theoreticalY(final double t) {
            final double theta = omega * t + alpha;
            return new double[] {
                r * FastMath.sin(theta), r * FastMath.cos(theta)
            };
        }

        public double[][] exactDyDy0(final double t) {

            // intermediate angle and state
            final double theta        = omega * t + alpha;
            final double sin          = FastMath.sin(theta);
            final double cos          = FastMath.cos(theta);
            final double y0           = r * sin;
            final double y1           = r * cos;

            return new double[][] {
                { dRdY00 * sin + y1 * dAlphadY00, dRdY01 * sin + y1 * dAlphadY01 },
                { dRdY00 * cos - y0 * dAlphadY00, dRdY01 * cos - y0 * dAlphadY01 }
            };

        }

        public double[] exactDyDomega(final double t) {

            // intermediate angle and state
            final double theta        = omega * t + alpha;
            final double sin          = FastMath.sin(theta);
            final double cos          = FastMath.cos(theta);
            final double y0           = r * sin;
            final double y1           = r * cos;

            return new double[] {
                y1 * (t + dAlphadOmega),
               -y0 * (t + dAlphadOmega)
            };

        }

        @Override
        public List<String> getParametersNames() {
            return Arrays.asList(OMEGA_PARAMETER);
        }

        @Override
        public boolean isSupported(final String name) {
            return OMEGA_PARAMETER.equals(name);
        }

    }

}
