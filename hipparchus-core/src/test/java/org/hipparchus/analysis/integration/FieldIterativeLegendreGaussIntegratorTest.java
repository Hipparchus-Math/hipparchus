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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.analysis.integration;

import java.util.Random;

import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.Decimal64;
import org.hipparchus.util.Decimal64Field;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;


public class FieldIterativeLegendreGaussIntegratorTest {

    @Test
    public void testSinFunction() {
        BaseAbstractFieldUnivariateIntegrator<Decimal64> integrator
            = new IterativeLegendreFieldGaussIntegrator<>(Decimal64Field.getInstance(), 5, 1.0e-14, 1.0e-10, 2, 15);

        Decimal64 min = new Decimal64(0);
        Decimal64 max = new Decimal64(FastMath.PI);
        double expected = 2;
        double tolerance = FastMath.max(integrator.getAbsoluteAccuracy(),
                                        FastMath.abs(expected * integrator.getRelativeAccuracy()));
        double result = integrator.integrate(10000, x -> x.sin(), min, max).getReal();
        Assert.assertEquals(expected, result, tolerance);

        min = new Decimal64(-FastMath.PI/3);
        max = new Decimal64(0);
        expected = -0.5;
        tolerance = FastMath.max(integrator.getAbsoluteAccuracy(),
                                 FastMath.abs(expected * integrator.getRelativeAccuracy()));
        result = integrator.integrate(10000, x -> x.sin(), min, max).getReal();
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testQuinticFunction() {
        CalculusFieldUnivariateFunction<Decimal64> f =
                        t -> t.subtract(1).multiply(t.subtract(0.5)).multiply(t).multiply(t.add(0.5)).multiply(t.add(1));
        FieldUnivariateIntegrator<Decimal64> integrator =
                new IterativeLegendreFieldGaussIntegrator<>(Decimal64Field.getInstance(), 3,
                                                            BaseAbstractUnivariateIntegrator.DEFAULT_RELATIVE_ACCURACY,
                                                            BaseAbstractUnivariateIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                                                            BaseAbstractUnivariateIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                                                            64);
        Decimal64 min = new Decimal64(0);
        Decimal64 max = new Decimal64(1);
        double expected = -1.0/48;
        double result = integrator.integrate(10000, f, min, max).getReal();
        Assert.assertEquals(expected, result, 1.0e-16);

        min = new Decimal64(0);
        max = new Decimal64(0.5);
        expected = 11.0/768;
        result = integrator.integrate(10000, f, min, max).getReal();
        Assert.assertEquals(expected, result, 1.0e-16);

        min = new Decimal64(-1);
        max = new Decimal64(4);
        expected = 2048/3.0 - 78 + 1.0/48;
        result = integrator.integrate(10000, f, min, max).getReal();
        Assert.assertEquals(expected, result, 4.0e-16 * expected);
    }

    @Test
    public void testExactIntegration() {
        Random random = new Random(86343623467878363l);
        for (int n = 2; n < 6; ++n) {
            IterativeLegendreFieldGaussIntegrator<Decimal64> integrator =
                new IterativeLegendreFieldGaussIntegrator<>(Decimal64Field.getInstance(), n,
                                                            BaseAbstractUnivariateIntegrator.DEFAULT_RELATIVE_ACCURACY,
                                                            BaseAbstractUnivariateIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                                                            BaseAbstractUnivariateIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                                                            64);

            // an n points Gauss-Legendre integrator integrates 2n-1 degree polynoms exactly
            for (int degree = 0; degree <= 2 * n - 1; ++degree) {
                for (int i = 0; i < 10; ++i) {
                    double[] coeff = new double[degree + 1];
                    for (int k = 0; k < coeff.length; ++k) {
                        coeff[k] = 2 * random.nextDouble() - 1;
                    }
                    PolynomialFunction p = new PolynomialFunction(coeff);
                    double result = integrator.integrate(10000,
                                                         p.toCalculusFieldUnivariateFunction(Decimal64Field.getInstance()),
                                                         new Decimal64(-5.0), new Decimal64(15.0)).getReal();
                    double reference = exactIntegration(p, -5.0, 15.0);
                    Assert.assertEquals(n + " " + degree + " " + i, reference, result, 1.0e-12 * (1.0 + FastMath.abs(reference)));
                }
            }

        }
    }

    // Cf. MATH-995
    @Test
    public void testNormalDistributionWithLargeSigma() {
        final Decimal64 sigma  = new Decimal64(1000);
        final Decimal64 mean   = new Decimal64(0);
        final Decimal64 factor = sigma.multiply(FastMath.sqrt(2 * FastMath.PI)).reciprocal();
        final Decimal64 i2s2   = sigma.multiply(sigma).reciprocal().multiply(0.5);
        final CalculusFieldUnivariateFunction<Decimal64> normal =
                        x -> x.subtract(mean).multiply(x.subtract(mean)).multiply(i2s2).negate().exp().multiply(factor);

        final double tol = 1e-2;
        final IterativeLegendreFieldGaussIntegrator<Decimal64> integrator =
            new IterativeLegendreFieldGaussIntegrator<>(Decimal64Field.getInstance(), 5, tol, tol);

        final Decimal64 a = new Decimal64(-5000);
        final Decimal64 b = new Decimal64(5000);
        final double s = integrator.integrate(50, normal, a, b).getReal();
        Assert.assertEquals(1, s, 1e-5);
    }

    @Test
    public void testIssue464() {
        final Decimal64 value = new Decimal64(0.2);
        CalculusFieldUnivariateFunction<Decimal64> f =
                        x -> (x.getReal() >= 0 && x.getReal() <= 5) ? value : Decimal64.ZERO;
        IterativeLegendreFieldGaussIntegrator<Decimal64> gauss
            = new IterativeLegendreFieldGaussIntegrator<>(Decimal64Field.getInstance(), 5, 3, 100);

        // due to the discontinuity, integration implies *many* calls
        double maxX = 0.32462367623786328;
        Assert.assertEquals(maxX * value.getReal(),
                            gauss.integrate(Integer.MAX_VALUE, f, new Decimal64(-10), new Decimal64(maxX)).getReal(),
                            1.0e-7);
        Assert.assertTrue(gauss.getEvaluations() > 37000000);
        Assert.assertTrue(gauss.getIterations() < 30);

        // setting up limits prevents such large number of calls
        try {
            gauss.integrate(1000, f, new Decimal64(-10), new Decimal64(maxX));
            Assert.fail("expected MathIllegalStateException");
        } catch (MathIllegalStateException tmee) {
            // expected
            Assert.assertEquals(LocalizedCoreFormats.MAX_COUNT_EXCEEDED, tmee.getSpecifier());
            Assert.assertEquals(1000, ((Integer) tmee.getParts()[0]).intValue());
        }

        // integrating on the two sides should be simpler
        double sum1 = gauss.integrate(1000, f, new Decimal64(-10), new Decimal64(0)).getReal();
        int eval1   = gauss.getEvaluations();
        double sum2 = gauss.integrate(1000, f, new Decimal64(0), new Decimal64(maxX)).getReal();
        int eval2   = gauss.getEvaluations();
        Assert.assertEquals(maxX * value.getReal(), sum1 + sum2, 1.0e-7);
        Assert.assertTrue(eval1 + eval2 < 200);

    }

    private double exactIntegration(PolynomialFunction p, double a, double b) {
        final double[] coeffs = p.getCoefficients();
        double yb = coeffs[coeffs.length - 1] / coeffs.length;
        double ya = yb;
        for (int i = coeffs.length - 2; i >= 0; --i) {
            yb = yb * b + coeffs[i] / (i + 1);
            ya = ya * a + coeffs[i] / (i + 1);
        }
        return yb * b - ya * a;
    }
}
