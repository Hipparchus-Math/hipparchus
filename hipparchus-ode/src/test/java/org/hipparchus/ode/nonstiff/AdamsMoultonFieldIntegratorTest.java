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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;

class AdamsMoultonFieldIntegratorTest extends AdamsFieldIntegratorAbstractTest {

    protected <T extends CalculusFieldElement<T>> AdamsFieldIntegrator<T>
    createIntegrator(Field<T> field, final int nSteps, final double minStep, final double maxStep,
                     final double scalAbsoluteTolerance, final double scalRelativeTolerance) {
        return new AdamsMoultonFieldIntegrator<>(field, nSteps, minStep, maxStep,
                scalAbsoluteTolerance, scalRelativeTolerance);
    }

    protected <T extends CalculusFieldElement<T>> AdamsFieldIntegrator<T>
    createIntegrator(Field<T> field, final int nSteps, final double minStep, final double maxStep,
                     final double[] vecAbsoluteTolerance, final double[] vecRelativeTolerance) {
        return new AdamsMoultonFieldIntegrator<>(field, nSteps, minStep, maxStep,
                vecAbsoluteTolerance, vecRelativeTolerance);
    }

    @Test
    void testNbPoints() {
        doNbPointsTest();
    }

    @Test
    public void testMinStep() {
        assertThrows(MathIllegalArgumentException.class, () -> doDimensionCheck(Binary64Field.getInstance()));
    }

    @Test
    public void testIncreasingTolerance() {
        // the 0.45 and 8.69 factors are only valid for this test
        // and has been obtained from trial and error
        // there are no general relationship between local and global errors
        doTestIncreasingTolerance(Binary64Field.getInstance(), 0.45, 8.69);
    }

    @Test
    public void exceedMaxEvaluations() {
        assertThrows(MathIllegalStateException.class, () -> doExceedMaxEvaluations(Binary64Field.getInstance(), 650));
    }

    @Test
    public void backward() {
        doBackward(Binary64Field.getInstance(), 3.0e-9, 3.0e-9, 1.0e-16, "Adams-Moulton");
    }

    @Test
    public void polynomial() {
        doPolynomial(Binary64Field.getInstance(), 5, 2.2e-05, 2.0e-11);
    }

    @Test
    public void testSecondaryEquations() {
        doTestSecondaryEquations(Binary64Field.getInstance(), 1.9e-11, 1.1e-14);
    }

    @Test
    public void testStartFailure() {
        assertThrows(MathIllegalStateException.class, () -> doTestStartFailure(Binary64Field.getInstance()));
    }

    @Test
    public void testOrder() {
        Assertions.assertEquals(1.97, estimateOrderestOrder(2).getReal(), 0.01);
        Assertions.assertEquals(2.95, estimateOrderestOrder(3).getReal(), 0.01);
        Assertions.assertEquals(3.92, estimateOrderestOrder(4).getReal(), 0.01);
        Assertions.assertEquals(4.92, estimateOrderestOrder(5).getReal(), 0.01);
        Assertions.assertEquals(5.89, estimateOrderestOrder(6).getReal(), 0.01);
    }

    private Binary64 estimateOrderestOrder(final int n) {
        FieldOrdinaryDifferentialEquation<Binary64> ode = new FieldOrdinaryDifferentialEquation<>() {

            /** {@inheritDoc} */
            @Override
            public int getDimension() {
                return 1;
            }

            /** {@inheritDoc} */
            @Override
            public Binary64[] computeDerivatives(final Binary64 t, final Binary64[] y) {
                return new Binary64[] {
                        y[0].multiply(FastMath.cos(t))
                };
            }
        };

        final FieldODEState<Binary64> initialState =
                new FieldODEState<>(Binary64.ZERO, new Binary64[] { Binary64.ONE });
        Binary64 prev = Binary64.NAN;
        Binary64 curr = Binary64.NAN;
        for (int k = 2; k <= 7; k++) {
            double h = 1.0 / FastMath.scalb(1.0, k);
            AdamsMoultonFieldIntegrator<Binary64> am =
                    new AdamsMoultonFieldIntegrator<>(Binary64Field.getInstance(), n, h, h, 1e3, 1e3);
            final Binary64 integrated = am.integrate(ode, initialState, Binary64.ONE).getPrimaryState()[0];
            prev = curr;
            curr = FastMath.abs(integrated.subtract(FastMath.exp(FastMath.sin(Binary64.ONE))));
        }
        return FastMath.log(prev.divide(curr)).divide(FastMath.log(new Binary64(2.0)));

    }

}
