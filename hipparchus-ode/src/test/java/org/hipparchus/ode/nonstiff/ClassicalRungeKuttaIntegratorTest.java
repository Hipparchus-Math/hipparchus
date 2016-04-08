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

package org.hipparchus.ode.nonstiff;


import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.TestProblem1;
import org.hipparchus.ode.TestProblem2;
import org.hipparchus.ode.TestProblem3;
import org.hipparchus.ode.TestProblem4;
import org.hipparchus.ode.TestProblem5;
import org.hipparchus.ode.TestProblem6;
import org.hipparchus.ode.TestProblemAbstract;
import org.hipparchus.ode.TestProblemHandler;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.ODEEventHandler;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.ode.sampling.StepInterpolatorTestUtils;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class ClassicalRungeKuttaIntegratorTest {

    @Test
    public void testMissedEndEvent()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        final double   t0     = 1878250320.0000029;
        final double   tEvent = 1878250379.9999986;
        final double[] k      = { 1.0e-4, 1.0e-5, 1.0e-6 };
        OrdinaryDifferentialEquation ode = new OrdinaryDifferentialEquation() {

            public int getDimension() {
                return k.length;
            }

            public double[] computeDerivatives(double t, double[] y) {
                double[] yDot = new double[y.length];
                for (int i = 0; i < y.length; ++i) {
                    yDot[i] = k[i] * y[i];
                }
                return yDot;
            }
        };

        ClassicalRungeKuttaIntegrator integrator = new ClassicalRungeKuttaIntegrator(60.0);

        double[] y0   = new double[k.length];
        for (int i = 0; i < y0.length; ++i) {
            y0[i] = i + 1;
        }
        double[] y    = new double[k.length];

        double finalT = integrator.integrate(ode, new ODEState(t0, y0), tEvent).getTime();
        Assert.assertEquals(tEvent, finalT, 5.0e-6);
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(y0[i] * FastMath.exp(k[i] * (finalT - t0)), y[i], 1.0e-9);
        }

        integrator.addEventHandler(new ODEEventHandler() {

            public double g(ODEStateAndDerivative s) {
                return s.getTime() - tEvent;
            }

            public Action eventOccurred(ODEStateAndDerivative s, boolean increasing) {
                Assert.assertEquals(tEvent, s.getTime(), 5.0e-6);
                return Action.CONTINUE;
            }
        }, Double.POSITIVE_INFINITY, 1.0e-20, 100);
        finalT = integrator.integrate(ode, new ODEState(t0, y0), tEvent + 120).getTime();
        Assert.assertEquals(tEvent + 120, finalT, 5.0e-6);
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(y0[i] * FastMath.exp(k[i] * (finalT - t0)), y[i], 1.0e-9);
        }

    }

    @Test
    public void testSanityChecks()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        try  {
            TestProblem1 pb = new TestProblem1();
            new ClassicalRungeKuttaIntegrator(0.01).integrate(pb,
                                                              new ODEState(0.0, new double[pb.getDimension()+10]),
                                                              1.0);
            Assert.fail("an exception should have been thrown");
        } catch(MathIllegalArgumentException ie) {
        }
        try  {
            TestProblem1 pb = new TestProblem1();
            new ClassicalRungeKuttaIntegrator(0.01).integrate(pb,
                                                              new ODEState(0.0, new double[pb.getDimension()]),
                                                              0.0);
            Assert.fail("an exception should have been thrown");
        } catch(MathIllegalArgumentException ie) {
        }
    }

    @Test
    public void testDecreasingSteps()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        for (TestProblemAbstract pb : new TestProblemAbstract[] {
                                                                 new TestProblem1(), new TestProblem2(), new TestProblem3(),
                                                                 new TestProblem4(), new TestProblem5(), new TestProblem6()
        }) {

            double previousValueError = Double.NaN;
            double previousTimeError = Double.NaN;
            for (int i = 4; i < 10; ++i) {

                double step = (pb.getFinalTime() - pb.getInitialTime()) * FastMath.pow(2.0, -i);

                ODEIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
                TestProblemHandler handler = new TestProblemHandler(pb, integ);
                integ.addStepHandler(handler);
                ODEEventHandler[] functions = pb.getEventsHandlers();
                for (int l = 0; l < functions.length; ++l) {
                    integ.addEventHandler(functions[l],
                                          Double.POSITIVE_INFINITY, 1.0e-6 * step, 1000);
                }
                Assert.assertEquals(functions.length, integ.getEventHandlers().size());
                double stopTime = integ.integrate(pb, pb.getInitialState(), pb.getFinalTime()).getTime();
                if (functions.length == 0) {
                    Assert.assertEquals(pb.getFinalTime(), stopTime, 1.0e-10);
                }

                double error = handler.getMaximalValueError();
                if (i > 4) {
                    Assert.assertTrue(error < 1.01 * FastMath.abs(previousValueError));
                }
                previousValueError = error;

                double timeError = handler.getMaximalTimeError();
                if (i > 4) {
                    Assert.assertTrue(timeError <= FastMath.abs(previousTimeError));
                }
                previousTimeError = timeError;

                integ.clearEventHandlers();
                Assert.assertEquals(0, integ.getEventHandlers().size());
            }

        }

    }

    @Test
    public void testSmallStep()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem1 pb = new TestProblem1();
        double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;

        ODEIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        Assert.assertTrue(handler.getLastError() < 2.0e-13);
        Assert.assertTrue(handler.getMaximalValueError() < 4.0e-12);
        Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
        Assert.assertEquals("classical Runge-Kutta", integ.getName());
    }

    @Test
    public void testBigStep()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem1 pb = new TestProblem1();
        double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.2;

        ODEIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        Assert.assertTrue(handler.getLastError() > 0.0004);
        Assert.assertTrue(handler.getMaximalValueError() > 0.005);
        Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);

    }

    @Test
    public void testBackward()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem5 pb = new TestProblem5();
        double step = FastMath.abs(pb.getFinalTime() - pb.getInitialTime()) * 0.001;

        ODEIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        Assert.assertTrue(handler.getLastError() < 5.0e-10);
        Assert.assertTrue(handler.getMaximalValueError() < 7.0e-10);
        Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
        Assert.assertEquals("classical Runge-Kutta", integ.getName());
    }

    @Test
    public void testKepler()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        final TestProblem3 pb  = new TestProblem3(0.9);
        double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.0003;

        ODEIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
        integ.addStepHandler(new KeplerHandler(pb));
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());
    }

    private static class KeplerHandler implements ODEStepHandler {
        public KeplerHandler(TestProblem3 pb) {
            this.pb = pb;
            maxError = 0;
        }
        public void init(ODEStateAndDerivative s0, double t) {
            maxError = 0;
        }
        public void handleStep(ODEStateInterpolator interpolator, boolean isLast)
                        throws MathIllegalStateException {

            double[] interpolatedY = interpolator.getCurrentState().getState();
            double[] theoreticalY  = pb.computeTheoreticalState(interpolator.getCurrentState().getTime());
            double dx = interpolatedY[0] - theoreticalY[0];
            double dy = interpolatedY[1] - theoreticalY[1];
            double error = dx * dx + dy * dy;
            if (error > maxError) {
                maxError = error;
            }
            if (isLast) {
                // even with more than 1000 evaluations per period,
                // RK4 is not able to integrate such an eccentric
                // orbit with a good accuracy
                Assert.assertTrue(maxError > 0.005);
            }
        }
        private double maxError = 0;
        private TestProblem3 pb;
    }

    @Test
    public void testStepSize() {
        final double step = 1.23456;
        ODEIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
        integ.addStepHandler(new ODEStepHandler() {
            public void handleStep(ODEStateInterpolator interpolator, boolean isLast) {
                if (! isLast) {
                    Assert.assertEquals(step,
                                        interpolator.getCurrentState().getTime() -
                                        interpolator.getPreviousState().getTime(),
                                        1.0e-12);
                }
            }
        });
        integ.integrate(new OrdinaryDifferentialEquation() {
            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { 1.0 };
            }
            public int getDimension() {
                return 1;
            }
        }, new ODEState(0.0, new double[] { 0.0 }), 5.0);
    }

    @Test
    public void testTooLargeFirstStep() {

        RungeKuttaIntegrator integ = new ClassicalRungeKuttaIntegrator(0.5);
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

        integ.integrate(equations, new ODEState(start, new double[] { 1.0 }), end);

    }

    @Test
    public void derivativesConsistency()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        TestProblem3 pb = new TestProblem3();
        double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;
        ClassicalRungeKuttaIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
        StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 6.6e-12);
    }

}
