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
import org.hipparchus.ode.AbstractIntegrator;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.TestProblem1;
import org.hipparchus.ode.TestProblem5;
import org.hipparchus.ode.TestProblem6;
import org.hipparchus.ode.TestProblemAbstract;
import org.hipparchus.ode.TestProblemHandler;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class AdamsMoultonIntegratorTest {

    @Test(expected=MathIllegalArgumentException.class)
    public void dimensionCheck()
        throws MathIllegalArgumentException, MathIllegalStateException {
        TestProblem1 pb = new TestProblem1();
        ODEIntegrator integ =
            new AdamsMoultonIntegrator(2, 0.0, 1.0, 1.0e-10, 1.0e-10);
        integ.integrate(pb,
                        new ODEState(0.0, new double[pb.getDimension()+10]),
                        1.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testMinStep()
            throws MathIllegalArgumentException, MathIllegalStateException {

          TestProblem1 pb = new TestProblem1();
          double minStep = 0.1 * (pb.getFinalTime() - pb.getInitialTime());
          double maxStep = pb.getFinalTime() - pb.getInitialTime();
          double[] vecAbsoluteTolerance = { 1.0e-15, 1.0e-16 };
          double[] vecRelativeTolerance = { 1.0e-15, 1.0e-16 };

          ODEIntegrator integ = new AdamsMoultonIntegrator(4, minStep, maxStep,
                                                           vecAbsoluteTolerance,
                                                           vecRelativeTolerance);
          TestProblemHandler handler = new TestProblemHandler(pb, integ);
          integ.addStepHandler(handler);
          integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

    }

    @Test
    public void testIncreasingTolerance()
            throws MathIllegalArgumentException, MathIllegalStateException {

        int previousCalls = Integer.MAX_VALUE;
        for (int i = -12; i < -2; ++i) {
            TestProblem1 pb = new TestProblem1();
            double minStep = 0;
            double maxStep = pb.getFinalTime() - pb.getInitialTime();
            double scalAbsoluteTolerance = FastMath.pow(10.0, i);
            double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

            ODEIntegrator integ = new AdamsMoultonIntegrator(4, minStep, maxStep,
                                                                    scalAbsoluteTolerance,
                                                                    scalRelativeTolerance);
            TestProblemHandler handler = new TestProblemHandler(pb, integ);
            integ.addStepHandler(handler);
            integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

            // the 0.45 and 8.69 factors are only valid for this test
            // and has been obtained from trial and error
            // there is no general relation between local and global errors
            Assert.assertTrue(handler.getMaximalValueError() > (0.45 * scalAbsoluteTolerance));
            Assert.assertTrue(handler.getMaximalValueError() < (8.69 * scalAbsoluteTolerance));
            Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-16);

            int calls = pb.getCalls();
            Assert.assertEquals(integ.getEvaluations(), calls);
            Assert.assertTrue(calls <= previousCalls);
            previousCalls = calls;

        }

    }

    @Test(expected = MathIllegalStateException.class)
    public void exceedMaxEvaluations()
            throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem1 pb  = new TestProblem1();
        double range = pb.getFinalTime() - pb.getInitialTime();

        AdamsMoultonIntegrator integ = new AdamsMoultonIntegrator(2, 0, range, 1.0e-12, 1.0e-12);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.setMaxEvaluations(650);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

    }

    @Test
    public void backward()
            throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem5 pb = new TestProblem5();
        double range = FastMath.abs(pb.getFinalTime() - pb.getInitialTime());

        ODEIntegrator integ = new AdamsMoultonIntegrator(4, 0, range, 1.0e-12, 1.0e-12);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        Assert.assertTrue(handler.getLastError() < 3.0e-9);
        Assert.assertTrue(handler.getMaximalValueError() < 3.0e-9);
        Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-16);
        Assert.assertEquals("Adams-Moulton", integ.getName());
    }

    @Test
    public void polynomial()
            throws MathIllegalArgumentException, MathIllegalStateException {
        TestProblem6 pb = new TestProblem6();
        double range = FastMath.abs(pb.getFinalTime() - pb.getInitialTime());

        for (int nSteps = 2; nSteps < 8; ++nSteps) {
            AdamsMoultonIntegrator integ =
                new AdamsMoultonIntegrator(nSteps, 1.0e-6 * range, 0.1 * range, 1.0e-5, 1.0e-5);
            integ.setStarterIntegrator(new PerfectStarter(pb, nSteps));
            TestProblemHandler handler = new TestProblemHandler(pb, integ);
            integ.addStepHandler(handler);
            integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());
            if (nSteps < 5) {
                Assert.assertTrue(handler.getMaximalValueError() > 2.2e-05);
            } else {
                Assert.assertTrue(handler.getMaximalValueError() < 1.1e-11);
            }
        }

    }

    private static class PerfectStarter extends AbstractIntegrator {

        private PerfectInterpolator interpolator;
        private final int nbSteps;

        public PerfectStarter(final TestProblemAbstract problem, final int nbSteps) {
            super("perfect");
            this.interpolator = new PerfectInterpolator(problem);
            this.nbSteps      = nbSteps;
        }

        public ODEStateAndDerivative integrate(ExpandableODE equations, ODEState initialState, double finalTime) {
            double tStart = initialState.getTime() + 0.01 * (finalTime - initialState.getTime());
            getEvaluationsCounter().increment(nbSteps);
            for (int i = 0; i < nbSteps; ++i) {
                double tK = ((nbSteps - 1 - (i + 1)) * initialState.getTime() + (i + 1) * tStart) / (nbSteps - 1);
                interpolator.setPreviousTime(interpolator.getCurrentTime());
                interpolator.setCurrentTime(tK);
                for (ODEStepHandler handler : getStepHandlers()) {
                    handler.handleStep(interpolator, i == nbSteps - 1);
                }
            }
            return interpolator.getCurrentState();
        }

    }

    private static class PerfectInterpolator implements ODEStateInterpolator {

        private static final long serialVersionUID = 20160403L;
        private final TestProblemAbstract problem;
        private double previousTime;
        private double currentTime;

        public PerfectInterpolator(final TestProblemAbstract problem) {
            this.problem      = problem;
            this.previousTime = problem.getInitialTime();
            this.currentTime  = problem.getInitialTime();
        }

        public void setPreviousTime(double time) {
            previousTime = time;
        }

        public double getCurrentTime() {
            return currentTime;
        }

        public void setCurrentTime(double time) {
            currentTime = time;
        }

        public ODEStateAndDerivative getPreviousState() {
            return getInterpolatedState(previousTime);
        }

        public ODEStateAndDerivative getCurrentState() {
            return getInterpolatedState(currentTime);
        }

        public ODEStateAndDerivative getInterpolatedState(final double time) {
            double[] y    = problem.computeTheoreticalState(time);
            double[] yDot = problem.computeDerivatives(time, y);
            return new ODEStateAndDerivative(time, y, yDot);
        }

        public boolean isForward() {
            return problem.getFinalTime() > problem.getInitialTime();
        }

    }

}
