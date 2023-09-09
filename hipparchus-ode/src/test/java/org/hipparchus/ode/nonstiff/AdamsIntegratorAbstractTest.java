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


import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.solvers.BracketedUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.AbstractIntegrator;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.LocalizedODEFormats;
import org.hipparchus.ode.MultistepIntegrator;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.SecondaryODE;
import org.hipparchus.ode.TestProblem1;
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
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public abstract class AdamsIntegratorAbstractTest {

    protected abstract AdamsIntegrator
    createIntegrator(final int nSteps, final double minStep, final double maxStep,
                     final double scalAbsoluteTolerance, final double scalRelativeTolerance);

    protected abstract AdamsIntegrator
    createIntegrator(final int nSteps, final double minStep, final double maxStep,
                     final double[] vecAbsoluteTolerance, final double[] vecRelativeTolerance);

    @Test(expected=MathIllegalArgumentException.class)
    public abstract void testMinStep();

    protected void doNbPointsTest() {
        try {
            createIntegrator(1, 1.0e-3, 1.0e+3, 1.0e-15, 1.0e-15);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedODEFormats.INTEGRATION_METHOD_NEEDS_AT_LEAST_TWO_PREVIOUS_POINTS,
                                miae.getSpecifier());
        }
        try {
            double[] vecAbsoluteTolerance = { 1.0e-15, 1.0e-16 };
            double[] vecRelativeTolerance = { 1.0e-15, 1.0e-16 };
            createIntegrator(1, 1.0e-3, 1.0e+3, vecAbsoluteTolerance, vecRelativeTolerance);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedODEFormats.INTEGRATION_METHOD_NEEDS_AT_LEAST_TWO_PREVIOUS_POINTS,
                                miae.getSpecifier());
        }
    }

    protected void doDimensionCheck() {
        TestProblem1 pb = new TestProblem1();

        double minStep = 0.1 * (pb.getFinalTime() - pb.getInitialState().getTime());
        double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
        double[] vecAbsoluteTolerance = { 1.0e-15, 1.0e-16 };
        double[] vecRelativeTolerance = { 1.0e-15, 1.0e-16 };

        ODEIntegrator integ = createIntegrator(4, minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

    }

    @Test
    public abstract void testIncreasingTolerance();

    protected void doTestIncreasingTolerance(double ratioMin, double ratioMax) {

        int previousCalls = Integer.MAX_VALUE;
        for (int i = -12; i < -2; ++i) {
            TestProblem1 pb = new TestProblem1();
            double minStep = 0;
            double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
            double scalAbsoluteTolerance = FastMath.pow(10.0, i);
            double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

            MultistepIntegrator integ = createIntegrator(4, minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
            int orderCorrection = integ instanceof AdamsBashforthIntegrator ? 0 : 1;
            Assert.assertEquals(FastMath.pow(2.0, 1.0 / (4 + orderCorrection)), integ.getMaxGrowth(), 1.0e-10);
            Assert.assertEquals(0.2, integ.getMinReduction(), 1.0e-10);
            Assert.assertEquals(4, integ.getNSteps());
            Assert.assertEquals(0.9, integ.getSafety(), 1.0e-10);
             Assert.assertTrue(integ.getStarterIntegrator() instanceof DormandPrince853Integrator);
            TestProblemHandler handler = new TestProblemHandler(pb, integ);
            integ.addStepHandler(handler);
            integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

            Assert.assertTrue(handler.getMaximalValueError() > ratioMin * scalAbsoluteTolerance);
            Assert.assertTrue(handler.getMaximalValueError() < ratioMax * scalAbsoluteTolerance);

            int calls = pb.getCalls();
            Assert.assertEquals(integ.getEvaluations(), calls);
            Assert.assertTrue(calls <= previousCalls);
            previousCalls = calls;

        }

    }

    @Test(expected = MathIllegalStateException.class)
    public abstract void exceedMaxEvaluations();

    protected void doExceedMaxEvaluations(final int max) {

        TestProblem1 pb  = new TestProblem1();
        double range = pb.getFinalTime() - pb.getInitialState().getTime();

        ODEIntegrator integ = createIntegrator(2, 0, range, 1.0e-12, 1.0e-12);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.setMaxEvaluations(max);
        integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

    }

    @Test
    public abstract void backward();

    protected void doBackward(final double epsilonLast,
                              final double epsilonMaxValue,
                              final double epsilonMaxTime,
                              final String name) {

        final double resetTime = -3.98;
        final TestProblem5 pb = new TestProblem5() {
            @Override
            public double[] getTheoreticalEventsTimes() {
                return new double[] { resetTime };
            }
        };
        final double range = pb.getFinalTime() - pb.getInitialState().getTime();

        AdamsIntegrator integ = createIntegrator(4, 0, range, 1.0e-12, 1.0e-12);
        ODEEventDetector event = new ODEEventDetector() {

            @Override
            public AdaptableInterval getMaxCheckInterval() {
                return s -> 0.5 * range;
            }

            @Override
            public int getMaxIterationCount() {
                return 100;
            }

            @Override
            public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
                return new BracketingNthOrderBrentSolver(0, 1.0e-6 * range, 0, 5);
            }

            @Override
            public ODEEventHandler getHandler() {
                return (state, detector, increasing) -> Action.RESET_STATE;
            }

            @Override
            public double g(ODEStateAndDerivative state) {
                return state.getTime() - resetTime;
            }

        };
        integ.addEventDetector(event);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

        Assert.assertEquals(0.0, handler.getLastError(), epsilonLast);
        Assert.assertEquals(0.0, handler.getMaximalValueError(), epsilonMaxValue);
        Assert.assertEquals(0, handler.getMaximalTimeError(), epsilonMaxTime);
        Assert.assertEquals(name, integ.getName());
    }

    @Test
    public abstract void polynomial();

    protected void doPolynomial(final int nLimit, final double epsilonBad, final double epsilonGood) {
        TestProblem6 pb = new TestProblem6();
        double range = FastMath.abs(pb.getFinalTime() - pb.getInitialState().getTime());

        for (int nSteps = 2; nSteps < 8; ++nSteps) {
            AdamsIntegrator integ = createIntegrator(nSteps, 1.0e-6 * range, 0.1 * range, 1.0e-4, 1.0e-4);
            ODEEventDetector event = new ODEEventDetector() {

                @Override
                public AdaptableInterval getMaxCheckInterval() {
                    return s -> 0.5 * range;
                }

                @Override
                public int getMaxIterationCount() {
                    return 100;
                }

                @Override
                public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
                    return new BracketingNthOrderBrentSolver(0, 1.0e-6 * range, 0, 5);
                }

                @Override
                public ODEEventHandler getHandler() {
                    return (state, detector, increasing) -> Action.RESET_STATE;
                }

                @Override
                public double g(ODEStateAndDerivative state) {
                    return state.getTime() - (pb.getInitialState().getTime() + 0.5 * range);
                }

            };
            integ.addEventDetector(event);
            integ.setStarterIntegrator(new PerfectStarter(pb, nSteps));
            TestProblemHandler handler = new TestProblemHandler(pb, integ);
            integ.addStepHandler(handler);
            integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());
            if (nSteps < nLimit) {
                Assert.assertTrue(handler.getMaximalValueError() > epsilonBad);
            } else {
                Assert.assertTrue(handler.getMaximalValueError() < epsilonGood);
            }
        }

    }

    @Test
    public void testNaNAppearing() {
        try {
            ODEIntegrator integ = createIntegrator(8, 0.01, 1.0, 0.1, 0.1);
            integ.integrate(new OrdinaryDifferentialEquation() {
                public int getDimension() {
                    return 1;
                }
                public double[] computeDerivatives(double t, double[] y) {
                    return new double[] { FastMath.log(t) };
                }
            }, new ODEState(10.0, new double[] { 1.0 }), -1.0);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedODEFormats.NAN_APPEARING_DURING_INTEGRATION, mise.getSpecifier());
            Assert.assertTrue(((Double) mise.getParts()[0]).doubleValue() <= 0.0);
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

        ODEIntegrator integrator = createIntegrator(6, 0.001, 1.0, 1.0e-12, 1.0e-12);
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

    @Test(expected=MathIllegalStateException.class)
    public abstract void testStartFailure();

    protected void doTestStartFailure() {
        TestProblem1 pb = new TestProblem1();
        double minStep = 0.0001 * (pb.getFinalTime() - pb.getInitialState().getTime());
        double maxStep = pb.getFinalTime() - pb.getInitialState().getTime();
        double scalAbsoluteTolerance = 1.0e-6;
        double scalRelativeTolerance = 1.0e-7;

        MultistepIntegrator integ = createIntegrator(6, minStep, maxStep,
                                                     scalAbsoluteTolerance, scalRelativeTolerance);
        integ.setStarterIntegrator(new DormandPrince853Integrator(maxStep * 0.5, maxStep, 0.1, 0.1));
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(new ExpandableODE(pb), pb.getInitialState(), pb.getFinalTime());

    }

    private static class PerfectStarter extends AbstractIntegrator {

        private final PerfectInterpolator interpolator;
        private final int nbSteps;

        public PerfectStarter(final TestProblemAbstract problem, final int nbSteps) {
            super("perfect-starter");
            this.interpolator = new PerfectInterpolator(problem);
            this.nbSteps      = nbSteps;
        }

        public ODEStateAndDerivative integrate(ExpandableODE equations,
                                               ODEState initialState, double finalTime) {
            double tStart = initialState.getTime() + 0.01 * (finalTime - initialState.getTime());
            getEvaluationsCounter().increment(nbSteps);
            interpolator.setCurrentTime(initialState.getTime());
            for (int i = 0; i < nbSteps; ++i) {
                double tK = ((nbSteps - 1 - (i + 1)) * initialState.getTime() + (i + 1) * tStart) /
                            (nbSteps - 1);
                interpolator.setPreviousTime(interpolator.getCurrentTime());
                interpolator.setCurrentTime(tK);
                for (ODEStepHandler handler : getStepHandlers()) {
                    handler.handleStep(interpolator);
                    if (i == nbSteps - 1) {
                        handler.finish(interpolator.getCurrentState());
                    }
                }
            }
            return interpolator.getInterpolatedState(tStart);
        }

    }

    private static class PerfectInterpolator implements ODEStateInterpolator {
        private static final long serialVersionUID = 20160417L;
        private final TestProblemAbstract problem;
        private double previousTime;
        private double currentTime;

        public PerfectInterpolator(final TestProblemAbstract problem) {
            this.problem = problem;
        }

        public void setPreviousTime(double previousTime) {
            this.previousTime = previousTime;
        }

        public void setCurrentTime(double currentTime) {
            this.currentTime = currentTime;
        }

        public double getCurrentTime() {
            return currentTime;
        }

        public boolean isForward() {
            return problem.getFinalTime() >= problem.getInitialState().getTime();
        }

        public ODEStateAndDerivative getPreviousState() {
            return getInterpolatedState(previousTime);
        }

        @Override
        public boolean isPreviousStateInterpolated() {
            return false;
        }

        public ODEStateAndDerivative getCurrentState() {
            return getInterpolatedState(currentTime);
        }

        @Override
        public boolean isCurrentStateInterpolated() {
            return false;
        }

        public ODEStateAndDerivative getInterpolatedState(double time) {
            double[] y    = problem.computeTheoreticalState(time);
            double[] yDot = problem.computeDerivatives(time, y);
            return new ODEStateAndDerivative(time, y, yDot);
        }

    }

}
