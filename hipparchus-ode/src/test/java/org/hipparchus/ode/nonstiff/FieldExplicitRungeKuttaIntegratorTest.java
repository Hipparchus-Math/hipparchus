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

import org.hipparchus.Field;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexField;
import org.hipparchus.ode.*;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.SinCos;
import org.junit.Assert;
import org.junit.Test;

public class FieldExplicitRungeKuttaIntegratorTest {

    @Test
    public void testFraction() {
        // GIVEN
        final int p = 2;
        final int q = 3;
        final Field<Complex> field = ComplexField.getInstance();
        // WHEN
        final Complex actualFraction = FieldExplicitRungeKuttaIntegrator.fraction(field, p, q);
        // THEN
        final Complex expectedFraction = FieldExplicitRungeKuttaIntegrator.fraction(field, (double) p, (double) q);
        Assert.assertEquals(expectedFraction.getReal(), actualFraction.getReal(), 0);
    }

    @Test
    public void testVersusNonField() {
        // GIVEN
        final TestFieldEquations testFieldEquations = new TestFieldEquations();
        final FieldExpandableODE<Complex> fieldExpandableODE = new FieldExpandableODE<>(testFieldEquations);

        final ComplexField field = ComplexField.getInstance();
        final Complex t0 = Complex.ZERO;
        final Complex[] y0 = new Complex[] { Complex.ONE, Complex.ZERO };
        final Complex h = new Complex(0.1);
        final FieldButcherArrayProvider<Complex> fieldExplicitRungeKutta = new ThreeEighthesFieldIntegrator<>(field,
                Complex.NaN);

        // WHEN
        final Complex[][] yDotK = MathArrays.buildArray(field, 4, 2);
        yDotK[0] = fieldExpandableODE.computeDerivatives(t0, y0);
        FieldExplicitRungeKuttaIntegrator.applyInternalButcherWeights(fieldExpandableODE, t0, y0, h,
                fieldExplicitRungeKutta.getA(), fieldExplicitRungeKutta.getC(), yDotK);
        final Complex[] actualState = FieldExplicitRungeKuttaIntegrator.applyExternalButcherWeights(y0, yDotK, h,
                fieldExplicitRungeKutta.getB());

        // THEN
        final double[] y0Real = new double[] { y0[0].getReal(), y0[1].getReal() };
        final double[][] yDotKReal = new double[yDotK.length][yDotK[0].length];
        final ExpandableODE expandableODE = new ExpandableODE(new TestEquations());
        yDotKReal[0] = expandableODE.computeDerivatives(t0.getReal(), y0Real);
        final ButcherArrayProvider explicitRungeKutta = new ThreeEighthesIntegrator(Double.NaN);
        ExplicitRungeKuttaIntegrator.applyInternalButcherWeights(expandableODE, t0.getReal(), y0Real, h.getReal(),
                explicitRungeKutta.getA(), explicitRungeKutta.getC(), yDotKReal);
        final double[] expectedState = ExplicitRungeKuttaIntegrator.applyExternalButcherWeights(y0Real, yDotKReal, h.getReal(),
                explicitRungeKutta.getB());
        for (int i = 0; i < expectedState.length; i++) {
            Assert.assertEquals(expectedState[i], actualState[i].getReal(), 0);
        }
    }

    @Test
    public void testRealCoefficientsVersusField() {
        // GIVEN
        final TestFieldEquations testFieldEquations = new TestFieldEquations();
        final FieldExpandableODE<Complex> fieldExpandableODE = new FieldExpandableODE<>(testFieldEquations);
        final ComplexField field = ComplexField.getInstance();
        final Complex t0 = Complex.ZERO;
        final Complex[] y0 = new Complex[] { Complex.ONE, Complex.ZERO };
        final Complex h = new Complex(0.1);
        final FieldExplicitRungeKuttaIntegrator<Complex> fieldExplicitRungeKutta = new ThreeEighthesFieldIntegrator<>(field,
                Complex.NaN);

        // WHEN
        final Complex[][] yDotK = MathArrays.buildArray(field, 4, 2);
        yDotK[0] = fieldExpandableODE.computeDerivatives(t0, y0);
        FieldExplicitRungeKuttaIntegrator.applyInternalButcherWeights(fieldExpandableODE, t0, y0, h,
                fieldExplicitRungeKutta.getA(), fieldExplicitRungeKutta.getC(), yDotK);
        final Complex[] actualState = FieldExplicitRungeKuttaIntegrator.applyExternalButcherWeights(y0, yDotK, h,
                fieldExplicitRungeKutta.getB());

        // THEN
        FieldExplicitRungeKuttaIntegrator.applyInternalButcherWeights(fieldExpandableODE, t0, y0, h,
                fieldExplicitRungeKutta.getRealA(), fieldExplicitRungeKutta.getRealC(), yDotK);
        final Complex[] expectedState = FieldExplicitRungeKuttaIntegrator.applyExternalButcherWeights(y0, yDotK, h,
                fieldExplicitRungeKutta.getRealB());
        for (int i = 0; i < expectedState.length; i++) {
            Assert.assertEquals(expectedState[i], actualState[i]);
        }
    }

    private static class TestFieldEquations implements FieldOrdinaryDifferentialEquation<Complex> {

        @Override
        public int getDimension() {
            return 2;
        }

        @Override
        public Complex[] computeDerivatives(Complex t, Complex[] y) {

            final FieldSinCos<Complex> sinCos = FastMath.sinCos(t);
            return new Complex[] { sinCos.sin(), sinCos.cos().multiply(y[0]) };
        }

    }

    private static class TestEquations implements OrdinaryDifferentialEquation {

        @Override
        public int getDimension() {
            return 2;
        }

        @Override
        public double[] computeDerivatives(double t, double[] y) {
            final SinCos sinCos = FastMath.sinCos(t);
            return new double[] { sinCos.sin(), sinCos.cos() * y[0] };
        }
    } 
    
}
