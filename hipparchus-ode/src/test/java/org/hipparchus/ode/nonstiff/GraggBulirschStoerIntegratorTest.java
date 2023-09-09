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

import org.hamcrest.MatcherAssert;
import org.hamcrest.Matchers;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.LocalizedODEFormats;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.TestProblem1;
import org.hipparchus.ode.TestProblem3;
import org.hipparchus.ode.TestProblem4;
import org.hipparchus.ode.TestProblem5;
import org.hipparchus.ode.TestProblemAbstract;
import org.hipparchus.ode.TestProblemHandler;
import org.hipparchus.ode.events.ODEEventDetector;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.ode.sampling.StepInterpolatorTestUtils;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;


public class GraggBulirschStoerIntegratorTest {

    @Test(expected=MathIllegalArgumentException.class)
    public void testDimensionCheck() {
        TestProblem1 pb = new TestProblem1();
        AdaptiveStepsizeIntegrator integrator =
                        new GraggBulirschStoerIntegrator(0.0, 1.0, 1.0e-10, 1.0e-10);
        integrator.integrate(pb,
                             new ODEState(0.0, new double[pb.getDimension()+10]),
                             1.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testNullIntervalCheck() {
        TestProblem1 pb = new TestProblem1();
        GraggBulirschStoerIntegrator integrator =
                        new GraggBulirschStoerIntegrator(0.0, 1.0, 1.0e-10, 1.0e-10);
        integrator.integrate(pb,
                             new ODEState(0.0, new double[pb.getDimension()]),
                             0.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testMinStep() {

        TestProblem5 pb  = new TestProblem5();
        double minStep   = 0.1 * FastMath.abs(pb.getFinalTime() - pb.getInitialTime());
        double maxStep   = FastMath.abs(pb.getFinalTime() - pb.getInitialTime());
        double[] vecAbsoluteTolerance = { 1.0e-20, 1.0e-21 };
        double[] vecRelativeTolerance = { 1.0e-20, 1.0e-21 };

        ODEIntegrator integ =
                        new GraggBulirschStoerIntegrator(minStep, maxStep,
                                                         vecAbsoluteTolerance, vecRelativeTolerance);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

    }

    @Test
    public void testBackward() {

        TestProblem5 pb = new TestProblem5();
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

        ODEIntegrator integ = new GraggBulirschStoerIntegrator(minStep, maxStep,
                                                               scalAbsoluteTolerance,
                                                               scalRelativeTolerance);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        Assert.assertTrue(handler.getLastError() < 7.5e-9);
        Assert.assertTrue(handler.getMaximalValueError() < 8.1e-9);
        Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
        Assert.assertEquals("Gragg-Bulirsch-Stoer", integ.getName());
    }

    @Test
    public void testIncreasingTolerance() {

        int previousCalls = Integer.MAX_VALUE;
        for (int i = -12; i < -4; ++i) {
            TestProblem1 pb     = new TestProblem1();
            double minStep      = 0;
            double maxStep      = pb.getFinalTime() - pb.getInitialTime();
            double absTolerance = FastMath.pow(10.0, i);
            double relTolerance = absTolerance;

            ODEIntegrator integ =
                            new GraggBulirschStoerIntegrator(minStep, maxStep,
                                                             absTolerance, relTolerance);
            TestProblemHandler handler = new TestProblemHandler(pb, integ);
            integ.addStepHandler(handler);
            integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

            // the coefficients are only valid for this test
            // and have been obtained from trial and error
            // there is no general relation between local and global errors
            double ratio =  handler.getMaximalValueError() / absTolerance;
            Assert.assertTrue(ratio < 2.4);
            Assert.assertTrue(ratio > 0.02);
            Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);

            int calls = pb.getCalls();
            Assert.assertEquals(integ.getEvaluations(), calls);
            Assert.assertTrue(calls <= previousCalls);
            previousCalls = calls;

        }

    }

    @Test
    public void testIntegratorControls() {

        TestProblem3 pb = new TestProblem3(0.999);
        GraggBulirschStoerIntegrator integ =
                        new GraggBulirschStoerIntegrator(0, pb.getFinalTime() - pb.getInitialTime(),
                                                         1.0e-8, 1.0e-10);

        double errorWithDefaultSettings = getMaxError(integ, pb);

        // stability control
        integ.setStabilityCheck(true, 2, 1, 0.99);
        Assert.assertTrue(errorWithDefaultSettings < getMaxError(integ, pb));
        integ.setStabilityCheck(true, -1, -1, -1);

        integ.setControlFactors(0.5, 0.99, 0.1, 2.5);
        Assert.assertTrue(errorWithDefaultSettings < getMaxError(integ, pb));
        integ.setControlFactors(-1, -1, -1, -1);

        integ.setOrderControl(10, 0.7, 0.95);
        Assert.assertTrue(errorWithDefaultSettings < getMaxError(integ, pb));
        integ.setOrderControl(-1, -1, -1);

        integ.setInterpolationControl(true, 3);
        Assert.assertTrue(errorWithDefaultSettings < getMaxError(integ, pb));
        integ.setInterpolationControl(true, -1);

    }

    private double getMaxError(ODEIntegrator integrator, TestProblemAbstract pb) {
        TestProblemHandler handler = new TestProblemHandler(pb, integrator);
        integrator.addStepHandler(handler);
        integrator.integrate(pb, pb.getInitialState(), pb.getFinalTime());
        return handler.getMaximalValueError();
    }

    @Test
    public void testEvents() {

        TestProblem4 pb = new TestProblem4();
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        double scalAbsoluteTolerance = 1.0e-10;
        double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

        ODEIntegrator integ = new GraggBulirschStoerIntegrator(minStep, maxStep,
                                                               scalAbsoluteTolerance,
                                                               scalRelativeTolerance);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        // since state is approx. linear at g=0 need convergence <= (state tolerance) / 2.
        double convergence = 1.0e-11;
        ODEEventDetector[] functions = pb.getEventDetectors(Double.POSITIVE_INFINITY, convergence, 1000);
        for (int l = 0; l < functions.length; ++l) {
            integ.addEventDetector(functions[l]);
        }
        Assert.assertEquals(functions.length, integ.getEventDetectors().size());
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        MatcherAssert.assertThat(handler.getMaximalValueError(), Matchers.lessThan(2.5e-11));
        // integration error builds up by the last event,
        // so tolerance is slightly more than the convergence.
        Assert.assertEquals(0, handler.getMaximalTimeError(), 1.5 * convergence);
        Assert.assertEquals(12.0, handler.getLastTime(), convergence);
        integ.clearEventDetectors();
        Assert.assertEquals(0, integ.getEventDetectors().size());

    }

    @Test
    public void testKepler() {

        final TestProblem3 pb = new TestProblem3(0.9);
        double minStep        = 0;
        double maxStep        = pb.getFinalTime() - pb.getInitialTime();
        double absTolerance   = 1.0e-6;
        double relTolerance   = 1.0e-6;

        ODEIntegrator integ =
                        new GraggBulirschStoerIntegrator(minStep, maxStep,
                                                         absTolerance, relTolerance);
        integ.addStepHandler(new KeplerStepHandler(pb));
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        Assert.assertEquals(integ.getEvaluations(), pb.getCalls());
        Assert.assertTrue(pb.getCalls() < 2150);

    }

    @Test
    public void testVariableSteps() {

        final TestProblem3 pb = new TestProblem3(0.9);
        double minStep        = 0;
        double maxStep        = pb.getFinalTime() - pb.getInitialTime();
        double absTolerance   = 1.0e-8;
        double relTolerance   = 1.0e-8;
        ODEIntegrator integ =
                        new GraggBulirschStoerIntegrator(minStep, maxStep,
                                                         absTolerance, relTolerance);
        integ.addStepHandler(new VariableStepHandler());
        double stopTime = integ.integrate(pb, pb.getInitialState(), pb.getFinalTime()).getTime();
        Assert.assertEquals(pb.getFinalTime(), stopTime, 1.0e-10);
        Assert.assertEquals("Gragg-Bulirsch-Stoer", integ.getName());
    }

    @Test
    public void testTooLargeFirstStep() {

        AdaptiveStepsizeIntegrator integ =
                        new GraggBulirschStoerIntegrator(0, Double.POSITIVE_INFINITY, Double.NaN, Double.NaN);
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
    public void testUnstableDerivative() {
        final StepProblem stepProblem = new StepProblem(s -> 999.0, 1.0e+12, 1000000, 0.0, 1.0, 2.0).
                                        withMaxCheck(1.0).
                                        withMaxIter(1000).
                                        withThreshold(1.0e-12);
        Assert.assertEquals(1.0,     stepProblem.getMaxCheckInterval().currentInterval(null), 1.0e-15);
        Assert.assertEquals(1000,    stepProblem.getMaxIterationCount());
        Assert.assertEquals(1.0e-12, stepProblem.getSolver().getAbsoluteAccuracy(), 1.0e-25);
        Assert.assertNotNull(stepProblem.getHandler());
        ODEIntegrator integ =
                        new GraggBulirschStoerIntegrator(0.1, 10, 1.0e-12, 0.0);
        integ.addEventDetector(stepProblem);
        Assert.assertEquals(8.0,
                            integ.integrate(stepProblem, new ODEState(0.0, new double[] { 0.0 }), 10.0).getPrimaryState()[0],
                            1.0e-12);
    }

    @Test
    public void derivativesConsistency()
        throws MathIllegalArgumentException, MathIllegalStateException {
        TestProblem3 pb = new TestProblem3(0.9);
        double minStep   = 0;
        double maxStep   = pb.getFinalTime() - pb.getInitialTime();
        double absTolerance = 1.0e-8;
        double relTolerance = 1.0e-8;

        GraggBulirschStoerIntegrator integ =
                        new GraggBulirschStoerIntegrator(minStep, maxStep,
                                                         absTolerance, relTolerance);
        StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 5.9e-10);
    }

    @Test
    public void testIssue596() {
        ODEIntegrator integ = new GraggBulirschStoerIntegrator(1e-10, 100.0, 1e-7, 1e-7);
        integ.addStepHandler(interpolator -> {
            double t = interpolator.getCurrentState().getTime();
            double[] y = interpolator.getInterpolatedState(t).getPrimaryState();
            double[] yDot = interpolator.getInterpolatedState(t).getPrimaryDerivative();
            Assert.assertEquals(3.0 * t - 5.0, y[0], 1.0e-14);
            Assert.assertEquals(3.0, yDot[0], 1.0e-14);
        });
        double[] y = {4.0};
        double t0 = 3.0;
        double tend = 10.0;
        integ.integrate(new OrdinaryDifferentialEquation() {
            public int getDimension() {
                return 1;
            }

            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { 3.0 };
            }
        }, new ODEState(t0, y), tend);

    }

    @Test
    public void testNaNAppearing() {
        try {
            ODEIntegrator integ = new GraggBulirschStoerIntegrator(0.01, 100.0, 1.0e5, 1.0e5);
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

    private static class KeplerStepHandler implements ODEStepHandler {
        public KeplerStepHandler(TestProblem3 pb) {
            this.pb = pb;
        }
        public void init(ODEStateAndDerivative s0, double t) {
            nbSteps = 0;
            maxError = 0;
        }
        public void handleStep(ODEStateInterpolator interpolator) {

            ++nbSteps;
            for (int a = 1; a < 100; ++a) {

                double prev   = interpolator.getPreviousState().getTime();
                double curr   = interpolator.getCurrentState().getTime();
                double interp = ((100 - a) * prev + a * curr) / 100;

                double[] interpolatedY = interpolator.getInterpolatedState (interp).getPrimaryState();
                double[] theoreticalY  = pb.computeTheoreticalState(interp);
                double dx = interpolatedY[0] - theoreticalY[0];
                double dy = interpolatedY[1] - theoreticalY[1];
                double error = dx * dx + dy * dy;
                if (error > maxError) {
                    maxError = error;
                }
            }
        }
        public void finish(ODEStateAndDerivative finalState) {
            Assert.assertTrue(maxError < 2.7e-6);
            Assert.assertTrue(nbSteps < 80);
        }
        private int nbSteps;
        private double maxError;
        private TestProblem3 pb;
    }

    public static class VariableStepHandler implements ODEStepHandler {
        private boolean firstTime;
        private double  minStep;
        private double  maxStep;
        public VariableStepHandler() {
            firstTime = true;
            minStep = 0;
            maxStep = 0;
        }
        public void init(double t0, double[] y0, double t) {
            firstTime = true;
            minStep = 0;
            maxStep = 0;
        }
        public void handleStep(ODEStateInterpolator interpolator) {

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
        }
        public void finish(ODEStateAndDerivative finalState) {
            Assert.assertTrue(minStep < 8.2e-3);
            Assert.assertTrue(maxStep > 1.5);
        }
    }

}
