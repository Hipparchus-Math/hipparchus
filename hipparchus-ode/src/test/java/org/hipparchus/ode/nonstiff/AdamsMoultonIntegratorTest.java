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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;

class AdamsMoultonIntegratorTest extends AdamsIntegratorAbstractTest {

    protected AdamsIntegrator
    createIntegrator(final int nSteps, final double minStep, final double maxStep,
                     final double scalAbsoluteTolerance, final double scalRelativeTolerance) {
        return new AdamsMoultonIntegrator(nSteps, minStep, maxStep,
                                          scalAbsoluteTolerance, scalRelativeTolerance);
    }

    protected AdamsIntegrator
    createIntegrator(final int nSteps, final double minStep, final double maxStep,
                     final double[] vecAbsoluteTolerance, final double[] vecRelativeTolerance) {
        return new AdamsMoultonIntegrator(nSteps, minStep, maxStep,
                                          vecAbsoluteTolerance, vecRelativeTolerance);
    }

    @Test
    void testNbPoints() {
        doNbPointsTest();
    }

    @Test
    public void testMinStep() {
        assertThrows(MathIllegalArgumentException.class, this::doDimensionCheck);
    }

    @Test
    public void testIncreasingTolerance() {
        // the 0.45 and 8.69 factors are only valid for this test
        // and has been obtained from trial and error
        // there are no general relationship between local and global errors
        doTestIncreasingTolerance(0.45, 8.69);
    }

    @Test
    public void exceedMaxEvaluations() {
        assertThrows(MathIllegalStateException.class, () -> doExceedMaxEvaluations(650));
    }

    @Test
    public void backward() {
        doBackward(3.0e-9, 3.0e-9, 1.0e-16, "Adams-Moulton");
    }

    @Test
    public void polynomial() {
        doPolynomial(5, 2.2e-05, 2.0e-11);
    }

    @Test
    public void testSecondaryEquations() {
        doTestSecondaryEquations(1.9e-11, 7.2e-15);
    }

    @Test
    public void testStartFailure() {
        assertThrows(MathIllegalStateException.class, this::doTestStartFailure);
    }

    @Test
    public void testOrder() {
        Assertions.assertEquals(1.97, estimateOrderestOrder(2), 0.01);
        Assertions.assertEquals(2.95, estimateOrderestOrder(3), 0.01);
        Assertions.assertEquals(3.92, estimateOrderestOrder(4), 0.01);
        Assertions.assertEquals(4.92, estimateOrderestOrder(5), 0.01);
        Assertions.assertEquals(5.89, estimateOrderestOrder(6), 0.01);
    }

    private double estimateOrderestOrder(final int n) {
        OrdinaryDifferentialEquation ode = new OrdinaryDifferentialEquation() {

            /** {@inheritDoc} */
            @Override
            public int getDimension() {
                return 1;
            }

            /** {@inheritDoc} */
            @Override
            public double[] computeDerivatives(final double t, final double[] y) {
                return new double[] {
                        y[0] * FastMath.cos(t)
                };
            }
        };

        final ODEState initialState =  new ODEState(0.0, new double[] { 1.0 });
        double prev = Double.NaN;
        double curr = Double.NaN;
        for (int k = 2; k <= 7; k++) {
            double h = 1.0 / FastMath.scalb(1.0, k);
            AdamsMoultonIntegrator am = new AdamsMoultonIntegrator(n, h, h, 1e3, 1e3);
            final double integrated = am.integrate(ode, initialState, 1.0).getPrimaryState()[0];
            prev = curr;
            curr = FastMath.abs(integrated - FastMath.exp(FastMath.sin(1.0)));
        }
        return FastMath.log(prev / curr) / FastMath.log(2);

    }

}
