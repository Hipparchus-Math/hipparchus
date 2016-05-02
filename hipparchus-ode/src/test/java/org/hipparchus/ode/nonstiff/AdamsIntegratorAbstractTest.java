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

package org.hipparchus.ode.nonstiff;


import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.ode.AbstractIntegrator;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.MultistepIntegrator;
import org.hipparchus.ode.TestProblem1;
import org.hipparchus.ode.TestProblem5;
import org.hipparchus.ode.TestProblem6;
import org.hipparchus.ode.TestProblemAbstract;
import org.hipparchus.ode.TestProblemHandler;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
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

            ODEIntegrator integ = createIntegrator(4, minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
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

        TestProblem5 pb = new TestProblem5();
        double range = pb.getFinalTime() - pb.getInitialState().getTime();

        AdamsIntegrator integ = createIntegrator(4, 0, range, 1.0e-12, 1.0e-12);
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
                    handler.handleStep(interpolator, i == nbSteps - 1);
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

        public ODEStateAndDerivative getCurrentState() {
            return getInterpolatedState(currentTime);
        }

        public ODEStateAndDerivative getInterpolatedState(double time) {
            double[] y    = problem.computeTheoreticalState(time);
            double[] yDot = problem.computeDerivatives(time, y);
            return new ODEStateAndDerivative(time, y, yDot);
        }

    }

}
