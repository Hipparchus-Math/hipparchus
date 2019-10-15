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
package org.hipparchus.ode.sampling;

import org.hipparchus.RealFieldElement;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEIntegrator;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.ode.TestFieldProblem3;
import org.hipparchus.ode.TestFieldProblemAbstract;
import org.hipparchus.ode.nonstiff.DormandPrince54FieldIntegrator;
import org.hipparchus.util.Decimal64;
import org.hipparchus.util.Decimal64Field;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;


public class FieldStepNormalizerTest {

    @Test
    public void testBoundariesDefault() {
        final TestFieldProblemAbstract<Decimal64> pb = new TestFieldProblem3<Decimal64>(new Decimal64(0.9));
        final Decimal64 range = pb.getFinalTime().subtract(pb.getInitialTime());
        final Decimal64 stepSize = range.divide(5.5);
        doTestBoundaries(pb, null, stepSize,
                         pb.getInitialTime(),
                         pb.getInitialTime().add(FastMath.floor(range.divide(stepSize)).multiply(stepSize)));
    }

    @Test
    public void testBoundariesNeither() {
        final TestFieldProblemAbstract<Decimal64> pb = new TestFieldProblem3<Decimal64>(new Decimal64(0.9));
        final Decimal64 range = pb.getFinalTime().subtract(pb.getInitialTime());
        final Decimal64 stepSize = range.divide(5.5);
        doTestBoundaries(pb, StepNormalizerBounds.NEITHER, stepSize,
                         pb.getInitialTime().add(stepSize),
                         pb.getInitialTime().add(FastMath.floor(range.divide(stepSize)).multiply(stepSize)));
    }

    @Test
    public void testBoundariesFirst() {
        final TestFieldProblemAbstract<Decimal64> pb = new TestFieldProblem3<Decimal64>(new Decimal64(0.9));
        final Decimal64 range = pb.getFinalTime().subtract(pb.getInitialTime());
        final Decimal64 stepSize = range.divide(5.5);
        doTestBoundaries(pb, StepNormalizerBounds.FIRST, stepSize,
                         pb.getInitialTime(),
                         pb.getInitialTime().add(FastMath.floor(range.divide(stepSize)).multiply(stepSize)));
    }

    @Test
    public void testBoundariesLast() {
        final TestFieldProblemAbstract<Decimal64> pb = new TestFieldProblem3<Decimal64>(new Decimal64(0.9));
        final Decimal64 range = pb.getFinalTime().subtract(pb.getInitialTime());
        final Decimal64 stepSize = range.divide(5.5);
        doTestBoundaries(pb, StepNormalizerBounds.LAST, stepSize,
                         pb.getInitialTime().add(stepSize),
                         pb.getFinalTime());
    }

    @Test
    public void testBoundariesBoth() {
        final TestFieldProblemAbstract<Decimal64> pb = new TestFieldProblem3<Decimal64>(new Decimal64(0.9));
        final Decimal64 range = pb.getFinalTime().subtract(pb.getInitialTime());
        final Decimal64 stepSize = range.divide(5.5);
        doTestBoundaries(pb, StepNormalizerBounds.BOTH, stepSize,
                         pb.getInitialTime(),
                         pb.getFinalTime());
    }

    @Test
    public void testBeforeEnd() {
        final TestFieldProblemAbstract<Decimal64> pb = new TestFieldProblem3<Decimal64>(new Decimal64(0.9));
        final Decimal64 range = pb.getFinalTime().subtract(pb.getInitialTime());
        final Decimal64 stepSize = range.divide(10.5);
        doTestBoundaries(pb, null, stepSize,
                         pb.getInitialTime(),
                         pb.getFinalTime().subtract(range.divide(21.0)));
    }

    @Test
    public void testModeForwardMultiples() {
        doTestStepsAtIntegerTimes(StepNormalizerMode.MULTIPLES, 2.0, 7.5, 2.5, 7.5, 4.0);
    }

    @Test
    public void testModeForwardIncrement() {
        doTestStepsAtIntegerTimes(StepNormalizerMode.INCREMENT, 2.0, 7.5, 2.5, 7.5, 3.5);
    }

    @Test
    public void testModeBackwardMultiples() {
        doTestStepsAtIntegerTimes(StepNormalizerMode.MULTIPLES, 2.0, 2.5, 7.5, 2.5, 6.0);
    }

    @Test
    public void testModeBackwardIncrement() {
        doTestStepsAtIntegerTimes(StepNormalizerMode.INCREMENT, 2.0, 2.5, 7.5, 2.5, 6.5);
    }

    private <T extends RealFieldElement<T>> void doTestBoundaries(final TestFieldProblemAbstract<T> pb,
                                                                  final StepNormalizerBounds bounds,
                                                                  final T stepSize,
                                                                  final T expectedFirst,
                                                                  final T expectedLast) {
        double minStep = 0;
        double maxStep = pb.getFinalTime().getReal() - pb.getInitialTime().getReal();
        FieldODEIntegrator<T> integ =
                        new DormandPrince54FieldIntegrator<>(stepSize.getField(),
                                                             minStep, maxStep, 10.e-8, 1.0e-8);
        final Checker<T> checker = new Checker<>();
        if (bounds == null) {            
            integ.addStepHandler(new FieldStepNormalizer<>(stepSize.getReal(), checker));
        } else {
            integ.addStepHandler(new FieldStepNormalizer<>(stepSize.getReal(), checker, bounds));
        }
        integ.integrate(new FieldExpandableODE<>(pb), pb.getInitialState(), pb.getFinalTime());
        Assert.assertEquals(expectedFirst.getReal(), checker.firstTime.getReal(), 1.0e-10);
        Assert.assertEquals(expectedLast.getReal(),  checker.lastTime.getReal(),  1.0e-10);
    }

    private void doTestStepsAtIntegerTimes(final StepNormalizerMode mode,
                                           final double stepSize,
                                           final double t0, final double t1,
                                           final double expectedFirst,
                                           final double expectedLast) {
        final TestFieldProblemAbstract<Decimal64> pb = new TestFieldProblem3<Decimal64>(new Decimal64(0.9));
        double minStep = 0;
        double maxStep = pb.getFinalTime().getReal() - pb.getInitialTime().getReal();
        FieldODEIntegrator<Decimal64> integ =
                        new DormandPrince54FieldIntegrator<>(Decimal64Field.getInstance(),
                                                             minStep, maxStep, 10.e-8, 1.0e-8);
        final Checker<Decimal64> checker = new Checker<>();
        integ.addStepHandler(new FieldStepNormalizer<Decimal64>(stepSize, checker, mode));
        integ.integrate(new FieldExpandableODE<>(new FieldOrdinaryDifferentialEquation<Decimal64>() {
            public int getDimension() { return 1; }
            public Decimal64[] computeDerivatives(Decimal64 t, Decimal64[] y) { return y; }
        }), new FieldODEState<>(new Decimal64(t0), new Decimal64[] { new Decimal64(0) }), new Decimal64(t1));
        Assert.assertEquals(expectedFirst, checker.firstTime.getReal(), 1.0e-10);
        Assert.assertEquals(expectedLast, checker.lastTime.getReal(),  1.0e-10);
    }

    private static class Checker<T extends RealFieldElement<T>> implements FieldODEFixedStepHandler<T> {

        private T firstTime = null;
        private T lastTime = null;

        public void handleStep(FieldODEStateAndDerivative<T> s, boolean isLast) {
            if (firstTime == null) {
                firstTime = s.getTime();
            }
            if (isLast) {
                lastTime = s.getTime();
            }
        }

    }

}
