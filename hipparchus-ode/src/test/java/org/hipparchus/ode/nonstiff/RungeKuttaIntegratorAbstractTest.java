/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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

package org.hipparchus.ode.nonstiff;


import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Random;

import org.hamcrest.MatcherAssert;
import org.hamcrest.Matchers;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.solvers.BracketedUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.DenseOutputModel;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.LocalizedODEFormats;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.SecondaryODE;
import org.hipparchus.ode.TestProblem1;
import org.hipparchus.ode.TestProblem2;
import org.hipparchus.ode.TestProblem3;
import org.hipparchus.ode.TestProblem4;
import org.hipparchus.ode.TestProblem5;
import org.hipparchus.ode.TestProblem6;
import org.hipparchus.ode.TestProblemAbstract;
import org.hipparchus.ode.TestProblemHandler;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.AdaptableInterval;
import org.hipparchus.ode.events.ODEEventDetector;
import org.hipparchus.ode.events.ODEEventHandler;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.ode.sampling.StepInterpolatorTestUtils;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public abstract class RungeKuttaIntegratorAbstractTest {

    protected abstract RungeKuttaIntegrator createIntegrator(double step);

    @Test
    public abstract void testMissedEndEvent();

    protected void doTestMissedEndEvent(final double epsilonT,
                                        final double epsilonY)
        throws MathIllegalArgumentException, MathIllegalStateException {
        final double   t0     = 1878250320.0000029;
        final double   tEvent = 1878250379.9999986;
        final double[] k      = new double[] { 1.0e-4, 1.0e-5, 1.0e-6 };
        OrdinaryDifferentialEquation ode = new OrdinaryDifferentialEquation() {

            public int getDimension() {
                return k.length;
            }

            public double[] computeDerivatives(double t, double[] y) {
                double[] yDot = new double[k.length];
                for (int i = 0; i < y.length; ++i) {
                    yDot[i] = k[i] * y[i];
                }
                return yDot;
            }
        };

        RungeKuttaIntegrator integrator = createIntegrator(60.0);

        double[] y0   = new double[k.length];
        for (int i = 0; i < y0.length; ++i) {
            y0[i] = i;
        }

        ODEStateAndDerivative result = integrator.integrate(new ExpandableODE(ode),
                                                            new ODEState(t0, y0),
                                                            tEvent);
        Assert.assertEquals(tEvent, result.getTime(), epsilonT);
        double[] y = result.getPrimaryState();
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(y0[i] * FastMath.exp(k[i] * (result.getTime() - t0)), y[i], epsilonY);
        }

        integrator.addEventDetector(new ODEEventDetector() {
            public AdaptableInterval getMaxCheckInterval() {
                return s-> Double.POSITIVE_INFINITY;
            }
            public int getMaxIterationCount() {
                return 100;
            }
            public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
                return new BracketingNthOrderBrentSolver(0, 1.0e-20, 0, 5);
            }
            public double g(ODEStateAndDerivative state) {
                return state.getTime() - tEvent;
            }
            public ODEEventHandler getHandler() {
                return (state, detector, increasing) -> {
                    Assert.assertEquals(tEvent, state.getTime(), epsilonT);
                    return Action.CONTINUE;
                };
            }
        });
        result = integrator.integrate(new ExpandableODE(ode),
                                      new ODEState(t0, y0),
                                      tEvent + 120);
        Assert.assertEquals(tEvent + 120, result.getTime(), epsilonT);
        y = result.getPrimaryState();
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(y0[i] * FastMath.exp(k[i] * (result.getTime() - t0)), y[i], epsilonY);
        }

    }

    @Test
    public abstract void testSanityChecks();

    protected void doTestSanityChecks() {
        RungeKuttaIntegrator integrator = createIntegrator(0.01);
        try  {
            TestProblem1 pb = new TestProblem1();
            integrator.integrate(new ExpandableODE(pb),
                                 new ODEState(0.0, new double[pb.getDimension() + 10]),
                                 1.0);
            Assert.fail("an exception should have been thrown");
        } catch(MathIllegalArgumentException ie) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, ie.getSpecifier());
        }
        try  {
            TestProblem1 pb = new TestProblem1();
            integrator.integrate(new ExpandableODE(pb),
                                 new ODEState(0.0, new double[pb.getDimension()]),
                                 0.0);
            Assert.fail("an exception should have been thrown");
        } catch(MathIllegalArgumentException ie) {
            Assert.assertEquals(LocalizedODEFormats.TOO_SMALL_INTEGRATION_INTERVAL, ie.getSpecifier());
        }
    }

    @Test
    public abstract void testDecreasingSteps();

    protected void doTestDecreasingSteps(final double safetyValueFactor,
                                         final double safetyTimeFactor,
                                         final double epsilonT)
        throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblemAbstract[] allProblems = new TestProblemAbstract[] {
            new TestProblem1(), new TestProblem2(), new TestProblem3(),
            new TestProblem4(), new TestProblem5(), new TestProblem6()
        };
        for (TestProblemAbstract pb :  allProblems) {

            double previousValueError = Double.NaN;
            double previousTimeError  = Double.NaN;
            for (int i = 4; i < 10; ++i) {

                double step = FastMath.scalb(pb.getFinalTime() - pb.getInitialState().getTime(), -i);

                RungeKuttaIntegrator integ = createIntegrator(step);
                TestProblemHandler handler = new TestProblemHandler(pb, integ);
                integ.addStepHandler(handler);
                double eventTol = 1.0e-6 * step;
                ODEEventDetector[] functions = pb.getEventDetectors(Double.POSITIVE_INFINITY, eventTol, 1000);
                for (int l = 0; l < functions.length; ++l) {
                    integ.addEventDetector(functions[l]);
                }
                Assert.assertEquals(functions.length, integ.getEventDetectors().size());
                ODEStateAndDerivative stop = integ.integrate(new ExpandableODE(pb),
                                                                     pb.getInitialState(),
                                                                     pb.getFinalTime());
                if (functions.length == 0) {
                    Assert.assertEquals(pb.getFinalTime(), stop.getTime(), epsilonT);
                }

                double error = handler.getMaximalValueError();
                if (i > 4) {
                    Assert.assertTrue(error < FastMath.abs(previousValueError * safetyValueFactor));
                }
                previousValueError = error;

                double timeError = handler.getMaximalTimeError();
                // can't expect time error to be less than event finding tolerance
                double timeTol = FastMath.max(eventTol, FastMath.abs(previousTimeError * safetyTimeFactor));
                if (i > 4) {
                    MatcherAssert.assertThat(
                            "Problem=" + pb + ", i=" + i + ", step=" + step,
                            timeError,
                            Matchers.lessThanOrEqualTo(timeTol));
                }
                previousTimeError = timeError;

                integ.clearEventDetectors();
                Assert.assertEquals(0, integ.getEventDetectors().size());
            }

        }

    }

    @Test
    public abstract void testSmallStep();

    protected void doTestSmallStep(final double epsilonLast,
                                   final double epsilonMaxValue,
                                   final double epsilonMaxTime,
                                   final String name) {

        TestProblem1 pb = new TestProblem1();
        double step = 0.001 * (pb.getFinalTime() - pb.getInitialState().getTime());

        RungeKuttaIntegrator integ = createIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

        Assert.assertEquals(0, handler.getLastError(),         epsilonLast);
        Assert.assertEquals(0, handler.getMaximalValueError(), epsilonMaxValue);
        Assert.assertEquals(0, handler.getMaximalTimeError(),  epsilonMaxTime);
        Assert.assertEquals(name, integ.getName());

    }

    @Test
    public abstract void testBigStep();

    protected void doTestBigStep(final double belowLast,
                                 final double belowMaxValue,
                                 final double epsilonMaxTime,
                                 final String name)
        throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem1 pb = new TestProblem1();
        double step = 0.2 * (pb.getFinalTime() - pb.getInitialState().getTime());

        RungeKuttaIntegrator integ = createIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

        Assert.assertTrue(handler.getLastError()         > belowLast);
        Assert.assertTrue(handler.getMaximalValueError() > belowMaxValue);
        Assert.assertEquals(0, handler.getMaximalTimeError(),  epsilonMaxTime);
        Assert.assertEquals(name, integ.getName());

    }

    @Test
    public abstract void testBackward();

    protected void doTestBackward(final double epsilonLast,
                                  final double epsilonMaxValue,
                                  final double epsilonMaxTime,
                                  final String name)
        throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem5 pb = new TestProblem5();
        double step = FastMath.abs(0.001 * (pb.getFinalTime() - pb.getInitialState().getTime()));

        RungeKuttaIntegrator integ = createIntegrator(step);
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

    protected void doTestKepler(double expectedMaxError, double epsilon)
        throws MathIllegalArgumentException, MathIllegalStateException {

        final TestProblem3 pb  = new TestProblem3(0.9);
        double step = 0.0003 * (pb.getFinalTime() - pb.getInitialState().getTime());

        RungeKuttaIntegrator integ = createIntegrator(step);
        integ.addStepHandler(new KeplerHandler(pb, expectedMaxError, epsilon));
        final ExpandableODE expandable = new ExpandableODE(pb);
        Assert.assertSame(pb, expandable.getPrimary());
        integ.integrate(expandable, pb.getInitialState(), pb.getFinalTime());
    }

    private static class KeplerHandler implements ODEStepHandler {
        private double maxError;
        private final TestProblem3 pb;
        private final double expectedMaxError;
        private final double epsilon;
        public KeplerHandler(TestProblem3 pb, double expectedMaxError, double epsilon) {
            this.pb               = pb;
            this.expectedMaxError = expectedMaxError;
            this.epsilon          = epsilon;
            this.maxError         = 0;
        }
        public void init(ODEStateAndDerivative state0, double t) {
            maxError = 0;
        }
        public void handleStep(ODEStateInterpolator interpolator) {

            ODEStateAndDerivative current = interpolator.getCurrentState();
            double[] theoreticalY  = pb.computeTheoreticalState(current.getTime());
            double dx = current.getPrimaryState()[0] - theoreticalY[0];
            double dy = current.getPrimaryState()[1] - theoreticalY[1];
            maxError = FastMath.max(maxError, dx * dx + dy * dy);
        }
        public void finish(ODEStateAndDerivative finalState) {
            Assert.assertEquals(expectedMaxError, maxError, epsilon);
        }
    }

    @Test
    public abstract void testStepSize();

    protected void doTestStepSize(final double epsilon)
        throws MathIllegalArgumentException, MathIllegalStateException {
        final double finalTime = 5.0;
        final double step = 1.23456;
        RungeKuttaIntegrator integ = createIntegrator(step);
        integ.addStepHandler(new ODEStepHandler() {
            public void handleStep(ODEStateInterpolator interpolator) {
                if (interpolator.getCurrentState().getTime() < finalTime - 0.001) {
                    Assert.assertEquals(step,
                                        interpolator.getCurrentState().getTime() - interpolator.getPreviousState().getTime(),
                                        epsilon);
                }
            }
        });
        integ.integrate(new ExpandableODE(new OrdinaryDifferentialEquation() {
            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { 1.0 };
            }
            public int getDimension() {
                return 1;
            }
        }), new ODEState(0, new double[1]), finalTime);
    }

    @Test
    public abstract void testSingleStep();

    protected void doTestSingleStep(final double epsilon) {

        final TestProblem3 pb  = new TestProblem3(0.9);
        double h = 0.0003 * (pb.getFinalTime() - pb.getInitialState().getTime());

        RungeKuttaIntegrator integ = createIntegrator(Double.NaN);
        double   t = pb.getInitialState().getTime();
        double[] y = pb.getInitialState().getPrimaryState();
        for (int i = 0; i < 100; ++i) {
            y  = integ.singleStep(pb, t, y, t + h);
            t += h;
        }
        double[] yth = pb.computeTheoreticalState(t);
        double dx    = y[0] - yth[0];
        double dy    = y[1] - yth[1];
        double error = dx * dx + dy * dy;
        Assert.assertEquals(0.0, error, epsilon);
    }

    @Test
    public abstract void testTooLargeFirstStep();

    protected void doTestTooLargeFirstStep() {

        RungeKuttaIntegrator integ = createIntegrator(0.5);
        final double   t0 = 0;
        final double[] y0 = new double[] { 1.0 };
        final double   t  = 0.001;
        OrdinaryDifferentialEquation equations = new OrdinaryDifferentialEquation() {

            public int getDimension() {
                return 1;
            }

            public double[] computeDerivatives(double t, double[] y) {
                Assert.assertTrue(t >= FastMath.nextAfter(t0, Double.NEGATIVE_INFINITY));
                Assert.assertTrue(t <= FastMath.nextAfter(t,  Double.POSITIVE_INFINITY));
                return new double[] { -100 * y[0] };
            }

        };

        integ.integrate(new ExpandableODE(equations), new ODEState(t0, y0), t);

    }

    @Test
    public abstract void testUnstableDerivative();

    protected void doTestUnstableDerivative(double epsilon) {
        final StepProblem stepProblem = new StepProblem(s -> 999.0, 1.0e+12, 1000000, 0.0, 1.0, 2.0).
                        withMaxCheck(1.0).
                        withMaxIter(1000).
                        withThreshold(1.0e-12);
        Assert.assertEquals(1.0,     stepProblem.getMaxCheckInterval().currentInterval(null), 1.0e-15);
        Assert.assertEquals(1000,    stepProblem.getMaxIterationCount());
        Assert.assertEquals(1.0e-12, stepProblem.getSolver().getAbsoluteAccuracy(), 1.0e-25);
        Assert.assertNotNull(stepProblem.getHandler());
        RungeKuttaIntegrator integ = createIntegrator(0.3);
      integ.addEventDetector(stepProblem);
      ODEStateAndDerivative result = integ.integrate(new ExpandableODE(stepProblem),
                                                     new ODEState(0, new double[1]),
                                                     10.0);
      Assert.assertEquals(8.0, result.getPrimaryState()[0], epsilon);
    }

    @Test
    public abstract void testDerivativesConsistency();

    protected void doTestDerivativesConsistency(double epsilon) {
        TestProblem3 pb = new TestProblem3();
        double step = 0.001 * (pb.getFinalTime() - pb.getInitialState().getTime());
        RungeKuttaIntegrator integ = createIntegrator(step);
        StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.001, 1.0e-10);
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
                // here, we compute only half of the derivative
                // we will compute the full derivatives by multiplying
                // the main equation from within the additional equation
                // it is not the proper way, but it is intended to check
                // additional equations *can* change main equation
                return new double[] { 0.5 * y[1], -0.5 * y[0] };
            }

        };

        SecondaryODE linear = new SecondaryODE() {

            @Override
            public int getDimension() {
                return 1;
            }

            @Override
            public double[] computeDerivatives(double t, double[] primary, double[] primaryDot, double[] secondary) {
                for (int i = 0; i < primaryDot.length; ++i) {
                    // this secondary equation also changes the primary state derivative
                    // a proper example of this is for example optimal control when
                    // the secondary equations handle co-state, which changes control,
                    // and the control changes the primary state
                    primaryDot[i] *= 2;
                }
                return new double[] { -1 };
            }

        };

        ExpandableODE expandable = new ExpandableODE(sinCos);
        expandable.addSecondaryEquations(linear);

        ODEIntegrator integrator = createIntegrator(0.001);
        final double[] max = new double[2];
        integrator.addStepHandler(new ODEStepHandler() {
            @Override
            public void handleStep(ODEStateInterpolator interpolator) {
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
            ODEIntegrator integ = createIntegrator(0.3);
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

    @Test
    public abstract void testSerialization();

    protected void doTestSerialization(int expectedSize, double tolerance) {
        try {
            TestProblem3 pb = new TestProblem3(0.9);
            double h = 0.0003 * (pb.getFinalTime() - pb.getInitialState().getTime());
            RungeKuttaIntegrator integ = createIntegrator(h);

            integ.addStepHandler(new DenseOutputModel());
            integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream    oos = new ObjectOutputStream(bos);
            for (ODEStepHandler handler : integ.getStepHandlers()) {
                oos.writeObject(handler);
            }

            Assert.assertTrue("size = " + bos.size (), bos.size () >  9 * expectedSize / 10);
            Assert.assertTrue("size = " + bos.size (), bos.size () < 11 * expectedSize / 10);

            ByteArrayInputStream  bis = new ByteArrayInputStream(bos.toByteArray());
            ObjectInputStream     ois = new ObjectInputStream(bis);
            DenseOutputModel cm  = (DenseOutputModel) ois.readObject();

            Random random = new Random(347588535632l);
            double maxError = 0.0;
            for (int i = 0; i < 1000; ++i) {
                double r = random.nextDouble();
                double time = r * pb.getInitialTime() + (1.0 - r) * pb.getFinalTime();
                double[] interpolatedY = cm.getInterpolatedState(time).getPrimaryState();
                double[] theoreticalY  = pb.computeTheoreticalState(time);
                double dx = interpolatedY[0] - theoreticalY[0];
                double dy = interpolatedY[1] - theoreticalY[1];
                double error = dx * dx + dy * dy;
                if (error > maxError) {
                    maxError = error;
                }
            }

            Assert.assertEquals(0, maxError, tolerance);

        } catch (IOException | ClassNotFoundException e) {
            Assert.fail(e.getLocalizedMessage());
        }

    }

    @Test
    public void testIssue250() {
        final double defaultStep = 60.;
        RungeKuttaIntegrator integrator = createIntegrator(defaultStep);
        Assert.assertEquals(defaultStep, integrator.getDefaultStep(), 0.);
    }

}
