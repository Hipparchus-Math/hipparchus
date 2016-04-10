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
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.TestProblem1;
import org.hipparchus.ode.TestProblem3;
import org.hipparchus.ode.TestProblem4;
import org.hipparchus.ode.TestProblem5;
import org.hipparchus.ode.TestProblemHandler;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.ODEEventHandler;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;


public class DormandPrince853IntegratorTest {

    @Test
    public void testMissedEndEvent()
                    throws MathIllegalArgumentException, MathIllegalStateException {
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

        DormandPrince853Integrator integrator = new DormandPrince853Integrator(0.0, 100.0,
                                                                               1.0e-10, 1.0e-10);

        double[] y0   = new double[k.length];
        for (int i = 0; i < y0.length; ++i) {
            y0[i] = i + 1;
        }

        integrator.setInitialStepSize(60.0);
        ODEStateAndDerivative finalState = integrator.integrate(new ExpandableODE(ode), new ODEState(t0, y0), tEvent);
        Assert.assertEquals(tEvent, finalState.getTime(), 5.0e-6);
        double[] y = finalState.getState();
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(y0[i] * FastMath.exp(k[i] * (finalState.getTime() - t0)), y[i], 1.0e-9);
        }

        integrator.setInitialStepSize(60.0);
        integrator.addEventHandler(new ODEEventHandler() {

            public double g(ODEStateAndDerivative s) {
                return s.getTime() - tEvent;
            }

            public Action eventOccurred(ODEStateAndDerivative s, boolean increasing) {
                Assert.assertEquals(tEvent, s.getTime(), 5.0e-6);
                return Action.CONTINUE;
            }
        }, Double.POSITIVE_INFINITY, 1.0e-20, 100);
        finalState = integrator.integrate(new ExpandableODE(ode), new ODEState(t0, y0), tEvent + 120);
        Assert.assertEquals(tEvent + 120, finalState.getTime(), 5.0e-6);
        y = finalState.getState();
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(y0[i] * FastMath.exp(k[i] * (finalState.getTime() - t0)), y[i], 1.0e-9);
        }

    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testDimensionCheck()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        TestProblem1 pb = new TestProblem1();
        DormandPrince853Integrator integrator = new DormandPrince853Integrator(0.0, 1.0,
                                                                               1.0e-10, 1.0e-10);
        integrator.integrate(pb,
                             new ODEState(0.0, new double[pb.getDimension()+10]),
                             1.0);
        Assert.fail("an exception should have been thrown");
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testNullIntervalCheck()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        TestProblem1 pb = new TestProblem1();
        DormandPrince853Integrator integrator = new DormandPrince853Integrator(0.0, 1.0,
                                                                               1.0e-10, 1.0e-10);
        integrator.integrate(pb,
                             new ODEState(0.0, new double[pb.getDimension()]),
                             0.0);
        Assert.fail("an exception should have been thrown");
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testMinStep()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem1 pb = new TestProblem1();
        double minStep = 0.1 * (pb.getFinalTime() - pb.getInitialTime());
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        double[] vecAbsoluteTolerance = { 1.0e-15, 1.0e-16 };
        double[] vecRelativeTolerance = { 1.0e-15, 1.0e-16 };

        ODEIntegrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                             vecAbsoluteTolerance,
                                                             vecRelativeTolerance);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());
        Assert.fail("an exception should have been thrown");

    }

    @Test
    public void testIncreasingTolerance()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        int previousCalls = Integer.MAX_VALUE;
        AdaptiveStepsizeIntegrator integ =
                        new DormandPrince853Integrator(0, Double.POSITIVE_INFINITY,
                                                       Double.NaN, Double.NaN);
        for (int i = -12; i < -2; ++i) {
            TestProblem1 pb = new TestProblem1();
            double minStep = 0;
            double maxStep = pb.getFinalTime() - pb.getInitialTime();
            double scalAbsoluteTolerance = FastMath.pow(10.0, i);
            double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;
            integ.setStepSizeControl(minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);

            TestProblemHandler handler = new TestProblemHandler(pb, integ);
            integ.addStepHandler(handler);
            integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

            // the 1.3 factor is only valid for this test
            // and has been obtained from trial and error
            // there is no general relation between local and global errors
            Assert.assertTrue(handler.getMaximalValueError() < (1.3 * scalAbsoluteTolerance));
            Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);

            int calls = pb.getCalls();
            Assert.assertEquals(integ.getEvaluations(), calls);
            Assert.assertTrue(calls <= previousCalls);
            previousCalls = calls;

        }

    }

    @Test
    public void testTooLargeFirstStep()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        AdaptiveStepsizeIntegrator integ =
                        new DormandPrince853Integrator(0, Double.POSITIVE_INFINITY, Double.NaN, Double.NaN);
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
    public void testBackward()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem5 pb = new TestProblem5();
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

        ODEIntegrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                             scalAbsoluteTolerance,
                                                             scalRelativeTolerance);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        Assert.assertTrue(handler.getLastError() < 1.1e-7);
        Assert.assertTrue(handler.getMaximalValueError() < 1.1e-7);
        Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
        Assert.assertEquals("Dormand-Prince 8 (5, 3)", integ.getName());
    }

    @Test
    public void testEvents()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem4 pb = new TestProblem4();
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        double scalAbsoluteTolerance = 1.0e-9;
        double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

        ODEIntegrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                             scalAbsoluteTolerance,
                                                             scalRelativeTolerance);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        ODEEventHandler[] functions = pb.getEventsHandlers();
        double convergence = 1.0e-8 * maxStep;
        for (int l = 0; l < functions.length; ++l) {
            integ.addEventHandler(functions[l], Double.POSITIVE_INFINITY, convergence, 1000);
        }
        Assert.assertEquals(functions.length, integ.getEventHandlers().size());
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        Assert.assertEquals(0, handler.getMaximalValueError(), 2.1e-7);
        Assert.assertEquals(0, handler.getMaximalTimeError(), convergence);
        Assert.assertEquals(12.0, handler.getLastTime(), convergence);
        integ.clearEventHandlers();
        Assert.assertEquals(0, integ.getEventHandlers().size());

    }

    @Deprecated
    @Test
    public void testDeprecatedInterfaces()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        final DeprecatedODE pb = new DeprecatedODE();
        double minStep = 0;
        double maxStep = 1.0;
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = scalAbsoluteTolerance;

        ODEIntegrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                             scalAbsoluteTolerance,
                                                             scalRelativeTolerance);
        final DeprecatedStepHandler stepHandler = new DeprecatedStepHandler(integ);
        integ.addStepHandler(stepHandler);
        final DeprecatedEventHandler eventHandler = new DeprecatedEventHandler();
        integ.addEventHandler(eventHandler, 1.0, 1.0e-10, 100);
        double[] y = new double[1];
        double finalT = integ.integrate(pb, 0.0, y, 10.0, y);

        Assert.assertTrue(pb.initCalled);
        Assert.assertTrue(stepHandler.initCalled);
        Assert.assertTrue(stepHandler.lastStepSeen);
        Assert.assertTrue(eventHandler.initCalled);
        Assert.assertTrue(eventHandler.resetCalled);

        Assert.assertEquals(4.0, finalT, 1.0e-10);

    }

    @Deprecated
    private static class DeprecatedODE
    implements org.hipparchus.ode.FirstOrderDifferentialEquations {

        private boolean initCalled = false;

        public void init(double t0, double[] y0, double finalTime) {
            initCalled = true;
        }

        public int getDimension() {
            return 1;
        }

        public void computeDerivatives(double t, double[] y, double[] yDot) {
            yDot[0] = 1.0;
        }

    }

    @Deprecated
    private static class DeprecatedStepHandler
    implements org.hipparchus.ode.sampling.StepHandler {

        private final ODEIntegrator integrator;
        private boolean initCalled   = false;
        private boolean lastStepSeen = false;

        public DeprecatedStepHandler(final ODEIntegrator integrator) {
            this.integrator = integrator;
        }

        public void init(double t0, double[] y0, double t) {
            initCalled = true;
        }

        public void handleStep(org.hipparchus.ode.sampling.StepInterpolator interpolator, boolean isLast) {
            Assert.assertEquals(interpolator.getPreviousTime(),
                                integrator.getCurrentStepStart(),
                                1.0e-10);
            if (isLast) {
                lastStepSeen = true;
            }
        }

    }

    @Deprecated
    private static class DeprecatedEventHandler
    implements org.hipparchus.ode.events.EventHandler {

        private boolean initCalled = false;
        private boolean resetCalled = false;

        public void init(double t0, double[] y0, double t) {
            initCalled = true;
        }

        public double g(double t, double[] y) {
            return (t - 2.0) * (t - 4.0);
        }

        public Action eventOccurred(double t, double[] y, boolean increasing) {
            return t < 3 ? Action.RESET_STATE : Action.STOP;
        }

        public void resetState(double t, double[] y) {
            resetCalled = true;
        }

    }

    @Test
    public void testVariableSteps()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        final TestProblem3 pb  = new TestProblem3(0.9);
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = scalAbsoluteTolerance;

        ODEIntegrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                             scalAbsoluteTolerance,
                                                             scalRelativeTolerance);
        integ.addStepHandler(new VariableHandler());
        double stopTime = integ.integrate(pb, pb.getInitialState(), pb.getFinalTime()).getTime();
        Assert.assertEquals(pb.getFinalTime(), stopTime, 1.0e-10);
        Assert.assertEquals("Dormand-Prince 8 (5, 3)", integ.getName());
    }

    @Test
    public void testUnstableDerivative()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        final StepProblem stepProblem = new StepProblem(0.0, 1.0, 2.0);
        ODEIntegrator integ = new DormandPrince853Integrator(0.1, 10, 1.0e-12, 0.0);
        integ.addEventHandler(stepProblem, 1.0, 1.0e-12, 1000);
        double[] y = { Double.NaN };
        integ.integrate(stepProblem, new ODEState(0.0, new double[] { 0.0 }), 10.0);
        Assert.assertEquals(8.0, y[0], 1.0e-12);
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

        ODEIntegrator integ = new DormandPrince853Integrator(0.001, 1.0, 1.0e-12, 0.0);
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
            return s.getState()[index];
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
        public VariableHandler() {
            firstTime = true;
            minStep = 0;
            maxStep = 0;
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
                Assert.assertTrue(minStep < (1.0 / 100.0));
                Assert.assertTrue(maxStep > (1.0 / 2.0));
            }
        }
        private boolean firstTime = true;
        private double  minStep = 0;
        private double  maxStep = 0;
    }

}
