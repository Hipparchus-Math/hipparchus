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

package org.hipparchus.ode.sampling;

import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.TestProblem3;
import org.hipparchus.ode.TestProblemAbstract;
import org.hipparchus.ode.nonstiff.DormandPrince54Integrator;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;


public class StepNormalizerTest {

    @Test
    public void testBoundariesDefault() {
        final TestProblemAbstract pb = new TestProblem3();
        final double range = pb.getFinalTime() - pb.getInitialTime();
        final double stepSize = range / 5.5;
        doTestBoundaries(pb, null, stepSize,
                         pb.getInitialTime(),
                         pb.getInitialTime() + FastMath.floor(range / stepSize) * stepSize);
    }

    @Test
    public void testBoundariesNeither() {
        final TestProblemAbstract pb = new TestProblem3();
        final double range = pb.getFinalTime() - pb.getInitialTime();
        final double stepSize = range / 5.5;
        doTestBoundaries(pb, StepNormalizerBounds.NEITHER, stepSize,
                         pb.getInitialTime() + stepSize,
                         pb.getInitialTime() + FastMath.floor(range / stepSize) * stepSize);
    }

    @Test
    public void testBoundariesFirst() {
        final TestProblemAbstract pb = new TestProblem3();
        final double range = pb.getFinalTime() - pb.getInitialTime();
        final double stepSize = range / 5.5;
        doTestBoundaries(pb, StepNormalizerBounds.FIRST, stepSize,
                         pb.getInitialTime(),
                         pb.getInitialTime() + FastMath.floor(range / stepSize) * stepSize);
    }

    @Test
    public void testBoundariesLast() {
        final TestProblemAbstract pb = new TestProblem3();
        final double range = pb.getFinalTime() - pb.getInitialTime();
        final double stepSize = range / 5.5;
        doTestBoundaries(pb, StepNormalizerBounds.LAST, stepSize,
                         pb.getInitialTime() + stepSize,
                         pb.getFinalTime());
    }

    @Test
    public void testBoundariesBoth() {
        final TestProblemAbstract pb = new TestProblem3();
        final double range = pb.getFinalTime() - pb.getInitialTime();
        final double stepSize = range / 5.5;
        doTestBoundaries(pb, StepNormalizerBounds.BOTH, stepSize,
                         pb.getInitialTime(),
                         pb.getFinalTime());
    }

    @Test
    public void testBeforeEnd() {
        final TestProblemAbstract pb = new TestProblem3();
        final double range = pb.getFinalTime() - pb.getInitialTime();
        final double stepSize = range / 10.5;
        doTestBoundaries(pb, null, stepSize,
                         pb.getInitialTime(),
                         pb.getFinalTime() - range / 21.0);
    }

    @Test
    public void testModeForwardMultiples() {
        doTestStepsAtIntegerTimes(new TestProblem3(), StepNormalizerMode.MULTIPLES,
                                  2.0, 7.5, 2.5, 7.5, 4.0);
    }

    @Test
    public void testModeForwardIncrement() {
        doTestStepsAtIntegerTimes(new TestProblem3(), StepNormalizerMode.INCREMENT,
                                  2.0, 7.5, 2.5, 7.5, 3.5);
    }

    @Test
    public void testModeBackwardMultiples() {
        doTestStepsAtIntegerTimes(new TestProblem3(), StepNormalizerMode.MULTIPLES,
                                  2.0, 2.5, 7.5, 2.5, 6.0);
    }

    @Test
    public void testModeBackwardIncrement() {
        doTestStepsAtIntegerTimes(new TestProblem3(), StepNormalizerMode.INCREMENT,
                                  2.0, 2.5, 7.5, 2.5, 6.5);
    }

    private void doTestBoundaries(final TestProblemAbstract pb,
                                  final StepNormalizerBounds bounds, final double stepSize,
                                  final double expectedFirst, final double expectedLast) {
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        ODEIntegrator integ = new DormandPrince54Integrator(minStep, maxStep, 10.e-8, 1.0e-8);
        final Checker checker = new Checker();
        if (bounds == null) {            
            integ.addStepHandler(new StepNormalizer(stepSize, checker));
        } else {
            integ.addStepHandler(new StepNormalizer(stepSize, checker, bounds));
        }
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());
        Assertions.assertEquals(expectedFirst, checker.firstTime, 1.0e-10);
        Assertions.assertEquals(expectedLast,  checker.lastTime,  1.0e-10);
    }

    private void doTestStepsAtIntegerTimes(final TestProblemAbstract pb,
                                           final StepNormalizerMode mode, final double stepSize,
                                           final double t0, final double t1,
                                           final double expectedFirst,
                                           final double expectedLast) {
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        ODEIntegrator integ = new DormandPrince54Integrator(minStep, maxStep, 10.e-8, 1.0e-8);
        final Checker checker = new Checker();
        integ.addStepHandler(new StepNormalizer(stepSize, checker, mode));
        integ.integrate(new OrdinaryDifferentialEquation() {
            public int getDimension() { return 1; }
            public double[] computeDerivatives(double t, double[] y) { return y; }
        }, new ODEState(t0, new double[1]), t1);
        Assertions.assertEquals(expectedFirst, checker.firstTime, 1.0e-10);
        Assertions.assertEquals(expectedLast, checker.lastTime,  1.0e-10);
    }

    private static class Checker implements ODEFixedStepHandler {

        private double firstTime;
        private double lastTime;

        public void init(final ODEStateAndDerivative initialState, final double finalTime) {
            firstTime = Double.NaN;
            lastTime  = Double.NaN;
        }

        public void handleStep(ODEStateAndDerivative s, boolean isLast) {
            if (Double.isNaN(firstTime)) {
                firstTime = s.getTime();
            }
            if (isLast) {
                lastTime = s.getTime();
            }
        }

    }

}
