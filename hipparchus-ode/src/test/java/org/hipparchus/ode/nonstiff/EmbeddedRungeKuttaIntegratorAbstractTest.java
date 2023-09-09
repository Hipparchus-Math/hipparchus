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


import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.solvers.BracketedUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
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
import org.hipparchus.ode.TestProblem7;
import org.hipparchus.ode.TestProblem8;
import org.hipparchus.ode.TestProblemHandler;
import org.hipparchus.ode.VariationalEquation;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.AdaptableInterval;
import org.hipparchus.ode.events.ODEEventDetector;
import org.hipparchus.ode.events.ODEEventHandler;
import org.hipparchus.ode.events.ODEStepEndHandler;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.util.CombinatoricsUtils;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.SinCos;
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
      double convergence = 1.0e-8 * maxStep;
      ODEEventDetector[] functions = pb.getEventDetectors(Double.POSITIVE_INFINITY, convergence, 1000);
      for (int l = 0; l < functions.length; ++l) {
          integ.addEventDetector(functions[l]);
      }
      List<ODEEventDetector> detectors = new ArrayList<>(integ.getEventDetectors());
      Assert.assertEquals(functions.length, detectors.size());

      for (int i = 0; i < detectors.size(); ++i) {
          Assert.assertSame(functions[i], detectors.get(i).getHandler());
          Assert.assertEquals(Double.POSITIVE_INFINITY, detectors.get(i).getMaxCheckInterval().currentInterval(null), 1.0);
          Assert.assertEquals(convergence, detectors.get(i).getSolver().getAbsoluteAccuracy(), 1.0e-15 * convergence);
          Assert.assertEquals(1000, detectors.get(i).getMaxIterationCount());
      }

        integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

      Assert.assertEquals(0, handler.getMaximalValueError(), epsilonMaxValue);
      Assert.assertEquals(0, handler.getMaximalTimeError(), convergence);
      Assert.assertEquals(12.0, handler.getLastTime(), convergence);
      Assert.assertEquals(name, integ.getName());
      integ.clearEventDetectors();
      Assert.assertEquals(0, integ.getEventDetectors().size());

    }

    @Test
    public abstract void testStepEnd();

    protected void doTestStepEnd(final int expectedCount, final String name) {
        TestProblem4 pb = new TestProblem4();
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

        ODEIntegrator integ = createIntegrator(minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
        double convergence = 1.0e-8 * maxStep;
        ODEEventDetector[] functions = pb.getEventDetectors(Double.POSITIVE_INFINITY, convergence, 1000);
        for (int l = 0; l < functions.length; ++l) {
            integ.addEventDetector(functions[l]);
        }
        List<ODEEventDetector> detectors = new ArrayList<>(integ.getEventDetectors());
        Assert.assertEquals(functions.length, detectors.size());

        for (int i = 0; i < detectors.size(); ++i) {
            Assert.assertSame(functions[i], detectors.get(i).getHandler());
            Assert.assertEquals(Double.POSITIVE_INFINITY, detectors.get(i).getMaxCheckInterval().currentInterval(null), 1.0);
            Assert.assertEquals(convergence, detectors.get(i).getSolver().getAbsoluteAccuracy(), 1.0e-15 * convergence);
            Assert.assertEquals(1000, detectors.get(i).getMaxIterationCount());
        }

        final StepCounter counter = new StepCounter(expectedCount + 10, Action.STOP);
        integ.addStepEndHandler(counter);
        Assert.assertEquals(1, integ.getStepEndHandlers().size());
        integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

        Assert.assertEquals(expectedCount, counter.count);
        Assert.assertEquals(name, integ.getName());
        integ.clearEventDetectors();
        Assert.assertEquals(0, integ.getEventDetectors().size());
        integ.clearStepEndHandlers();
        Assert.assertEquals(0, integ.getStepEndHandlers().size());
    }

    @Test
    public abstract void testStopAfterStep();

    protected void doTestStopAfterStep(final int count, final double expectedTime) {
        TestProblem4 pb = new TestProblem4();
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

        ODEIntegrator integ = createIntegrator(minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
        double convergence = 1.0e-8 * maxStep;
        ODEEventDetector[] functions = pb.getEventDetectors(Double.POSITIVE_INFINITY, convergence, 1000);
        for (int l = 0; l < functions.length; ++l) {
            integ.addEventDetector(functions[l]);
        }
        List<ODEEventDetector> detectors = new ArrayList<>(integ.getEventDetectors());
        Assert.assertEquals(functions.length, detectors.size());

        for (int i = 0; i < detectors.size(); ++i) {
            Assert.assertSame(functions[i], detectors.get(i).getHandler());
            Assert.assertEquals(Double.POSITIVE_INFINITY, detectors.get(i).getMaxCheckInterval().currentInterval(null), 1.0);
            Assert.assertEquals(convergence, detectors.get(i).getSolver().getAbsoluteAccuracy(), 1.0e-15 * convergence);
            Assert.assertEquals(1000, detectors.get(i).getMaxIterationCount());
        }

        final StepCounter counter = new StepCounter(count, Action.STOP);
        integ.addStepEndHandler(counter);
        Assert.assertEquals(1, integ.getStepEndHandlers().size());
        ODEStateAndDerivative finalState = integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

        Assert.assertEquals(count, counter.count);
        Assert.assertEquals(expectedTime, finalState.getTime(), 1.0e-6);

    }

    @Test
    public abstract void testResetAfterStep();

    protected void doTestResetAfterStep(final int resetCount, final int expectedCount) {
        TestProblem4 pb = new TestProblem4();
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

        ODEIntegrator integ = createIntegrator(minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
        double convergence = 1.0e-8 * maxStep;
        ODEEventDetector[] functions = pb.getEventDetectors(Double.POSITIVE_INFINITY, convergence, 1000);
        for (int l = 0; l < functions.length; ++l) {
            integ.addEventDetector(functions[l]);
        }
        List<ODEEventDetector> detectors = new ArrayList<>(integ.getEventDetectors());
        Assert.assertEquals(functions.length, detectors.size());

        for (int i = 0; i < detectors.size(); ++i) {
            Assert.assertSame(functions[i], detectors.get(i).getHandler());
            Assert.assertEquals(Double.POSITIVE_INFINITY, detectors.get(i).getMaxCheckInterval().currentInterval(null), 1.0);
            Assert.assertEquals(convergence, detectors.get(i).getSolver().getAbsoluteAccuracy(), 1.0e-15 * convergence);
            Assert.assertEquals(1000, detectors.get(i).getMaxIterationCount());
        }

        final StepCounter counter = new StepCounter(resetCount, Action.RESET_STATE);
        integ.addStepEndHandler(counter);
        Assert.assertEquals(1, integ.getStepEndHandlers().size());
        ODEStateAndDerivative finalState = integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

        Assert.assertEquals(expectedCount, counter.count);
        Assert.assertEquals(12.0, finalState.getTime(), 1.0e-6); // this corresponds to the Stop event detector
        for (int i = 0; i < finalState.getPrimaryStateDimension(); ++i) {
            Assert.assertEquals(0.0, finalState.getPrimaryState()[i], 1.0e-15);
            Assert.assertEquals(0.0, finalState.getPrimaryDerivative()[i], 1.0e-15);
        }

    }

    private static class StepCounter implements ODEStepEndHandler {
        final int    max;
        final Action actionAtMax;
        int          count;
        StepCounter(final int max, final Action actionAtMax) {
            this.max         = max;
            this.actionAtMax = actionAtMax;
            this.count       = 0;
        }
        public Action stepEndOccurred(final ODEStateAndDerivative state, final boolean forward) {
            return ++count == max ? actionAtMax : Action.CONTINUE;
        }
        public ODEState resetState(ODEStateAndDerivative state) {
            return new ODEState(state.getTime(),
                                new double[state.getPrimaryStateDimension()]);
        }
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

            integ.addEventDetector(new ODEEventDetector() {
                public AdaptableInterval getMaxCheckInterval() {
                    return s-> Double.POSITIVE_INFINITY;
                }
                public int getMaxIterationCount() {
                    return 1000;
                }
                public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
                    return new BracketingNthOrderBrentSolver(0, 1.0e-8 * maxStep, 0, 5);
                }
                public ODEEventHandler getHandler() {
                    return (state, detector, increasing) -> Action.CONTINUE;
                }
                public double g(ODEStateAndDerivative state) {
                    double middle = 0.5 * (pb.getInitialState().getTime() + pb.getFinalTime());
                    double offset = state.getTime() - middle;
                    if (offset > 0) {
                        throw new LocalException();
                    }
                    return offset;
                }
            });

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

        integ.addEventDetector(new ODEEventDetector() {
            public AdaptableInterval getMaxCheckInterval() {
                return s-> Double.POSITIVE_INFINITY;
            }
            public int getMaxIterationCount() {
                return 3;
            }
            public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
                return new BracketingNthOrderBrentSolver(0, 1.0e-8 * maxStep, 0, 5);
            }
            public ODEEventHandler getHandler() {
                return (state, detector, increasing) -> Action.CONTINUE;
            }
            public double g(ODEStateAndDerivative state) {
                double middle = 0.5 * (pb.getInitialState().getTime() + pb.getFinalTime());
                double offset = state.getTime() - middle;
                return (offset > 0) ? offset + 0.5 : offset - 0.5;
            }
        });

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
        double minStep = 1.0e-10;
        double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
        double[] vecAbsoluteTolerance = { 1.0e-8, 1.0e-8, 1.0e-10, 1.0e-10 };
        double[] vecRelativeTolerance = { 1.0e-10, 1.0e-10, 1.0e-8, 1.0e-8 };

        AdaptiveStepsizeIntegrator integ = createIntegrator(minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);

        try {
            Method getStepSizeHelper = AdaptiveStepsizeIntegrator.class.getDeclaredMethod("getStepSizeHelper", (Class<?>[]) null);
            getStepSizeHelper.setAccessible(true);
            StepsizeHelper helper = (StepsizeHelper) getStepSizeHelper.invoke(integ, (Object[]) null);
            integ.setInitialStepSize(-999);
            Assert.assertEquals(-1.0, helper.getInitialStep(), 1.0e-10);
            integ.setInitialStepSize(+999);
            Assert.assertEquals(-1.0, helper.getInitialStep(), 1.0e-10);
        } catch (NoSuchMethodException | IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
            Assert.fail(e.getLocalizedMessage());
        }

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
        public void handleStep(ODEStateInterpolator interpolator) {

            ODEStateAndDerivative current = interpolator.getCurrentState();
            double[] theoreticalY  = pb.computeTheoreticalState(current.getTime());
            double dx = current.getPrimaryState()[0] - theoreticalY[0];
            double dy = current.getPrimaryState()[1] - theoreticalY[1];
            double error = dx * dx + dy * dy;
            if (error > maxError) {
                maxError = error;
            }
        }
        public void finish(ODEStateAndDerivative finalState) {
            Assert.assertEquals(0.0, maxError, epsilon);
        }
    }

    @Test
    public abstract void testTorqueFreeMotionOmegaOnly();

    protected void doTestTorqueFreeMotionOmegaOnly(double epsilon) {

        final TestProblem7 pb  = new TestProblem7();
        double minStep = 1.0e-10;
        double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
        double[] vecAbsoluteTolerance = { 1.0e-8, 1.0e-8, 1.0e-8 };
        double[] vecRelativeTolerance = { 1.0e-10, 1.0e-10, 1.0e-10 };

        EmbeddedRungeKuttaIntegrator integ = createIntegrator(minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
        integ.addStepHandler(new TorqueFreeOmegaOnlyHandler(pb, epsilon));
        integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());
    }

    private static class TorqueFreeOmegaOnlyHandler implements ODEStepHandler {
        private double maxError;
        private final TestProblem7 pb;
        private final double epsilon;
        public TorqueFreeOmegaOnlyHandler(TestProblem7 pb, double epsilon) {
            this.pb      = pb;
            this.epsilon = epsilon;
            maxError     = 0;
        }
        public void init(ODEStateAndDerivative state0, double t) {
            maxError = 0;
        }
        public void handleStep(ODEStateInterpolator interpolator) {

            ODEStateAndDerivative current = interpolator.getCurrentState();
            double[] theoreticalY  = pb.computeTheoreticalState(current.getTime());
            double do1   = current.getPrimaryState()[0] - theoreticalY[0];
            double do2   = current.getPrimaryState()[1] - theoreticalY[1];
            double do3   = current.getPrimaryState()[2] - theoreticalY[2];
            double error = do1 * do1 + do2 * do2 + do3 * do3;
            if (error > maxError) {
                maxError = error;
            }
        }
        public void finish(ODEStateAndDerivative finalState) {
            Assert.assertEquals(0.0, maxError, epsilon);
        }
    }

    @Test
    /** Compare that the analytical model and the numerical model that compute the  the quaternion in a torque-free configuration give same results.
     * This test is used to validate the results of the analytical model as defined by Landau & Lifchitz.
     */
    public abstract void testTorqueFreeMotion();

    protected void doTestTorqueFreeMotion(double epsilonOmega, double epsilonQ) {

        final double   t0          =  0.0;
        final double   t1          = 20.0;
        final Vector3D omegaBase   = new Vector3D(5.0, 0.0, 4.0);
        final Rotation rBase       = new Rotation(0.9, 0.437, 0.0, 0.0, true);
        final List<Double> inertiaBase = Arrays.asList(3.0 / 8.0, 1.0 / 2.0, 5.0 / 8.0);
        double minStep = 1.0e-10;
        double maxStep = t1 - t0;
        double[] vecAbsoluteTolerance = { 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14 };
        double[] vecRelativeTolerance = { 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14 };

        // prepare a stream of problems to integrate
        Stream<TestProblem8> problems = Stream.empty();

        // add all possible permutations of the base rotation rate
        problems = Stream.concat(problems,
                                 permute(omegaBase).
                                 map(omega -> new TestProblem8(t0, t1, omega, rBase,
                                                               inertiaBase.get(0), Vector3D.PLUS_I,
                                                               inertiaBase.get(1), Vector3D.PLUS_J,
                                                               inertiaBase.get(2), Vector3D.PLUS_K)));

        // add all possible permutations of the base rotation
        problems = Stream.concat(problems,
                                 permute(rBase).
                                 map(r -> new TestProblem8(t0, t1, omegaBase, r,
                                                           inertiaBase.get(0), Vector3D.PLUS_I,
                                                           inertiaBase.get(1), Vector3D.PLUS_J,
                                                           inertiaBase.get(2), Vector3D.PLUS_K)));

        // add all possible permutations of the base inertia
        problems = Stream.concat(problems,
                                 permute(inertiaBase).
                                 map(inertia -> new TestProblem8(t0, t1, omegaBase, rBase,
                                                                 inertia.get(0), Vector3D.PLUS_I,
                                                                 inertia.get(1), Vector3D.PLUS_J,
                                                                 inertia.get(2), Vector3D.PLUS_K)));

        problems.forEach(problem -> {
            EmbeddedRungeKuttaIntegrator integ = createIntegrator(minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);   
            integ.addStepHandler(new TorqueFreeHandler(problem, epsilonOmega, epsilonQ));
            integ.integrate(new ExpandableODE(problem), problem.getInitialState(), problem.getFinalTime());
        });

    }

    @Test
    public abstract void testTorqueFreeMotionIssue230();

    protected void doTestTorqueFreeMotionIssue230(double epsilonOmega, double epsilonQ) {

        double i1   = 3.0 / 8.0;
        Vector3D a1 = Vector3D.PLUS_I;
        double i2   = 5.0 / 8.0;
        Vector3D a2 = Vector3D.PLUS_K;
        double i3   = 1.0 / 2.0;
        Vector3D a3 = Vector3D.PLUS_J;
        Vector3D o0 = new Vector3D(5.0, 0.0, 4.0);
        double o1   = Vector3D.dotProduct(o0, a1);
        double o2   = Vector3D.dotProduct(o0, a2);
        double o3   = Vector3D.dotProduct(o0, a3);
        double e    = 0.5 * (i1 * o1 * o1 + i2 * o2 * o2 + i3 * o3 * o3);
        double r1   = FastMath.sqrt(2 * e * i1);
        double r3   = FastMath.sqrt(2 * e * i3);
        int n = 50;
        for (int i = 0; i < n; ++i) {
            SinCos sc = FastMath.sinCos(-0.5 * FastMath.PI * (i + 50) / 200);
            Vector3D om = new Vector3D(r1 * sc.cos() / i1, a1, r3 * sc.sin() / i3, a3);
            TestProblem8 problem = new TestProblem8(0, 20,
                                                    om,
                                                    new Rotation(0.9, 0.437, 0.0, 0.0, true),
                                                    i1, a1,
                                                    i2, a2,
                                                    i3, a3);

            double minStep = 1.0e-10;
            double maxStep = problem.getFinalTime() - problem.getInitialTime();
            double[] vecAbsoluteTolerance = { 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14 };
            double[] vecRelativeTolerance = { 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-14 };

            EmbeddedRungeKuttaIntegrator integ = createIntegrator(minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
            integ.addStepHandler(new TorqueFreeHandler(problem, epsilonOmega, epsilonQ));
            integ.integrate(new ExpandableODE(problem), problem.getInitialState(), problem.getFinalTime() * 0.1);

        }

    }

    /** Generate all permutations of vector coordinates.
     * @param v vector to permute
     * @return permuted vector
     */
    private Stream<Vector3D> permute(final Vector3D v) {
        return CombinatoricsUtils.
                        permutations(Arrays.asList(v.getX(), v.getY(), v.getZ())).
                        map(a -> new Vector3D(a.get(0), a.get(1), a.get(2)));
    }

    /** Generate all permutations of rotation coordinates.
     * @param r rotation to permute
     * @return permuted rotation
     */
    private Stream<Rotation> permute(final Rotation r) {
        return CombinatoricsUtils.
                        permutations(Arrays.asList(r.getQ0(), r.getQ1(), r.getQ2(), r.getQ3())).
                        map(a -> new Rotation(a.get(0), a.get(1), a.get(2), a.get(3), false));
    }

    /** Generate all permutations of a list.
     * @param list list to permute
     * @return permuted list
     */
    private Stream<List<Double>> permute(final List<Double> list) {
        return CombinatoricsUtils.permutations(list);
    }

    private static class TorqueFreeHandler implements ODEStepHandler {
        private double maxErrorOmega;
        private double maxErrorQ;
        private final TestProblem8 pb;
        private final double epsilonOmega;
        private final double epsilonQ;
        private double outputStep;
        private double current;

        public TorqueFreeHandler(TestProblem8 pb, double epsilonOmega, double epsilonQ) {
            this.pb           = pb;
            this.epsilonOmega = epsilonOmega;
            this.epsilonQ     = epsilonQ;
            maxErrorOmega     = 0;
            maxErrorQ         = 0;
            outputStep        = 0.01;
        }

        public void init(ODEStateAndDerivative state0, double t) {
            maxErrorOmega = 0;
            maxErrorQ     = 0;
            current       = state0.getTime() - outputStep;
        }

        public void handleStep(ODEStateInterpolator interpolator) {

            current += outputStep;
            while (interpolator.getPreviousState().getTime() <= current &&
                   interpolator.getCurrentState().getTime() > current) {
                ODEStateAndDerivative state = interpolator.getInterpolatedState(current);
                final double[] theoretical  = pb.computeTheoreticalState(state.getTime());
                final double errorOmega = Vector3D.distance(new Vector3D(state.getPrimaryState()[0],
                                                                         state.getPrimaryState()[1],
                                                                         state.getPrimaryState()[2]),
                                                            new Vector3D(theoretical[0],
                                                                         theoretical[1],
                                                                         theoretical[2]));
                maxErrorOmega = FastMath.max(maxErrorOmega, errorOmega);
                final double errorQ = Rotation.distance(new Rotation(state.getPrimaryState()[3],
                                                                     state.getPrimaryState()[4],
                                                                     state.getPrimaryState()[5],
                                                                     state.getPrimaryState()[6],
                                                                     true),
                                                        new Rotation(theoretical[3],
                                                                     theoretical[4],
                                                                     theoretical[5],
                                                                     theoretical[6],
                                                                     true));
                maxErrorQ = FastMath.max(maxErrorQ, errorQ);
                current += outputStep;

            }
        }

        public void finish(ODEStateAndDerivative finalState) {
            Assert.assertEquals(0.0, maxErrorOmega, epsilonOmega);
            Assert.assertEquals(0.0, maxErrorQ,     epsilonQ);
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
            public ODEEventHandler getHandler() {
                return (state, detector, increasing) -> Action.CONTINUE;
            }
            public double g(ODEStateAndDerivative s) {
                return s.getTime() - tEvent;
            }
        });
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
        final StepProblem stepProblem = new StepProblem(s -> 999.0, 1.0e+12, 1000000, 0.0, 1.0, 2.0).
                        withMaxCheck(1.0).
                        withMaxIter(1000).
                        withThreshold(1.0e-12);
        Assert.assertEquals(1.0,     stepProblem.getMaxCheckInterval().currentInterval(null), 1.0e-15);
        Assert.assertEquals(1000,    stepProblem.getMaxIterationCount());
        Assert.assertEquals(1.0e-12, stepProblem.getSolver().getAbsoluteAccuracy(), 1.0e-25);
        Assert.assertNotNull(stepProblem.getHandler());
        ODEIntegrator integ = createIntegrator(0.1, 10, 1.0e-12, 0.0);
        integ.addEventDetector(stepProblem);
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
        integ.addEventDetector(sinChecker);
        integ.addStepHandler(sinChecker);
        integ.addEventDetector(cosChecker);
        integ.addStepHandler(cosChecker);
        double   t0 = 0.5;
        double[] y0 = new double[] { FastMath.sin(t0), FastMath.cos(t0) };
        double   t  = 10.0;
        integ.integrate(sincos, new ODEState(t0, y0), t);

    }

    private static class SchedulingChecker implements ODEStepHandler, ODEEventDetector {

        int index;
        double tMin;

        public SchedulingChecker(int index) {
            this.index = index;
        }

        public AdaptableInterval getMaxCheckInterval() {
            return s-> 0.01;
        }

        public int getMaxIterationCount() {
            return 100;
        }

        public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
            return new BracketingNthOrderBrentSolver(0, 1.0e-7, 0, 5);
        }

        public void init(ODEStateAndDerivative s0, double t) {
            tMin = s0.getTime();
        }

        public void handleStep(ODEStateInterpolator interpolator) {
            tMin = interpolator.getCurrentState().getTime();
        }

        public double g(ODEStateAndDerivative s) {
            // once a step has been handled by handleStep,
            // events checking should only refer to dates after the step
            Assert.assertTrue(s.getTime() >= tMin);
            return s.getPrimaryState()[index];
        }

        public ODEEventHandler getHandler() {
            return new ODEEventHandler() {
                public Action eventOccurred(ODEStateAndDerivative s, ODEEventDetector detector, boolean increasing) {
                    return Action.RESET_STATE;
                }

                public ODEStateAndDerivative resetState(ODEEventDetector detector, ODEStateAndDerivative s) {
                    // in fact, we don't need to reset anything for the test
                    return s;
                }
            };
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
            Assert.assertEquals(min, minStep, 0.01 * min);
            Assert.assertEquals(max, maxStep, 0.01 * max);
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
        SinCosJacobiansProvider sinCos = new SinCosJacobiansProvider(omega);

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

        double[] dydp0 = ve.extractParameterJacobian(finalState, SinCosJacobiansProvider.OMEGA_PARAMETER);
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

    @Test
    public void testInfiniteIntegration() {
        ODEIntegrator integ = createIntegrator(0.01, 1.0, 0.1, 0.1);
        TestProblem1 pb = new TestProblem1();
        double convergence = 1e-6;
        integ.addEventDetector(new ODEEventDetector() {
            @Override
            public AdaptableInterval getMaxCheckInterval() {
                return s-> Double.POSITIVE_INFINITY;
            }
            @Override
            public int getMaxIterationCount() {
                return 1000;
            }
            @Override
            public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
                return new BracketingNthOrderBrentSolver(0, convergence, 0, 5);
            }
            @Override
            public double g(ODEStateAndDerivative state) {
                return state.getTime() - pb.getFinalTime();
            }
            @Override
            public ODEEventHandler getHandler() {
                return (state, detector, increasing) -> Action.STOP;
            }
        });
        ODEStateAndDerivative finalState = integ.integrate(pb, pb.getInitialState(), Double.POSITIVE_INFINITY);
        Assert.assertEquals(pb.getFinalTime(), finalState.getTime(), convergence);
    }

    private static class SinCosJacobiansProvider implements ODEJacobiansProvider {

        public static String OMEGA_PARAMETER = "omega";

        private final double omega;
        private       double r;
        private       double alpha;

        private double dRdY00;
        private double dRdY01;
        private double dAlphadOmega;
        private double dAlphadY00;
        private double dAlphadY01;

        protected SinCosJacobiansProvider(final double omega) {
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