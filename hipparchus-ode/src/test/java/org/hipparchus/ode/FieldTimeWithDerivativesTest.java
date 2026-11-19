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
package org.hipparchus.ode;

import org.hipparchus.analysis.differentiation.Gradient;
import org.hipparchus.analysis.differentiation.GradientField;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.ode.nonstiff.ClassicalRungeKuttaFieldIntegrator;
import org.hipparchus.ode.nonstiff.DormandPrince54FieldIntegrator;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

/** Field integration where the time itself carries derivatives.
 * <p>
 * The time checks of the integrator are about the time value only;
 * they must not depend on the derivatives carried by the time.
 * </p>
 */
class FieldTimeWithDerivativesTest {

    /** y' = -y. */
    private static class Decay implements FieldOrdinaryDifferentialEquation<Gradient> {
        @Override
        public int getDimension() {
            return 1;
        }
        @Override
        public Gradient[] computeDerivatives(final Gradient t, final Gradient[] y) {
            return new Gradient[] { y[0].negate() };
        }
    }

    private FieldODEStateAndDerivative<Gradient> integrate(final double seed, final double t0, final double t1) {
        final GradientField field = GradientField.getField(2);
        final DormandPrince54FieldIntegrator<Gradient> integrator =
                        new DormandPrince54FieldIntegrator<>(field, 1.0e-6, 1.0, 1.0e-10, 1.0e-10);
        final FieldODEState<Gradient> initial =
                        new FieldODEState<>(new Gradient(t0, seed, 0.0), new Gradient[] { field.getOne() });
        return integrator.integrate(new Decay(), initial, new Gradient(t1, 0.0, seed));
    }

    @Test
    void testTooSmallIntervalDetectedWhateverTimeDerivatives() {
        for (final double seed : new double[] { 0.0, 1.0, 1.0e3 }) {
            assertThrows(MathIllegalArgumentException.class, () -> integrate(seed, 1.0, 1.0 + 1.0e-15));
        }
    }

    @Test
    void testFinalStateIndependentOfTimeDerivatives() {
        final FieldODEStateAndDerivative<Gradient> ref = integrate(0.0, 0.0, 2.0);
        assertEquals(FastMath.exp(-2.0), ref.getPrimaryState()[0].getValue(), 1.0e-9);
        for (final double seed : new double[] { 1.0, 1.0e3 }) {
            final FieldODEStateAndDerivative<Gradient> s = integrate(seed, 0.0, 2.0);
            assertEquals(2.0, s.getTime().getValue(), 0.0);
            assertEquals(ref.getPrimaryState()[0].getValue(), s.getPrimaryState()[0].getValue(),
                         1.0e-15 * ref.getPrimaryState()[0].getValue());
        }
    }

    @Test
    void testLastStepDetectedWhateverTimeDerivatives() {
        final int refSteps = countSteps(0.0);
        for (final double seed : new double[] { 1.0, 1.0e3 }) {
            assertEquals(refSteps, countSteps(seed));
        }
    }

    /** Count steps of a fixed-step integration whose times carry derivatives. */
    private int countSteps(final double seed) {
        final GradientField field = GradientField.getField(2);
        final ClassicalRungeKuttaFieldIntegrator<Gradient> integrator =
                        new ClassicalRungeKuttaFieldIntegrator<>(field, new Gradient(0.1, 1.0e-3 * seed, 0.0));
        final int[] steps = new int[1];
        integrator.addStepHandler(interpolator -> ++steps[0]);
        final FieldODEState<Gradient> initial =
                        new FieldODEState<>(new Gradient(0.0, seed, 0.0), new Gradient[] { field.getOne() });
        integrator.integrate(new Decay(), initial, new Gradient(1.0, 0.0, seed));
        return steps[0];
    }

}
