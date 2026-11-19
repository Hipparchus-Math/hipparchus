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
package org.hipparchus.analysis.integration;

import org.hipparchus.analysis.QuinticFunction;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.function.Gaussian;
import org.hipparchus.analysis.function.Sin;
import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;


class IterativeLegendreGaussIntegratorTest {

    @Test
    void testSinFunction() {
        UnivariateFunction f = new Sin();
        BaseAbstractUnivariateIntegrator integrator
            = new IterativeLegendreGaussIntegrator(5, 1.0e-14, 1.0e-10, 2, 15);
        double min, max, expected, result, tolerance;

        min = 0; max = FastMath.PI; expected = 2;
        tolerance = FastMath.max(integrator.getAbsoluteAccuracy(),
                             FastMath.abs(expected * integrator.getRelativeAccuracy()));
        result = integrator.integrate(10000, f, min, max);
        assertEquals(expected, result, tolerance);

        min = -FastMath.PI/3; max = 0; expected = -0.5;
        tolerance = FastMath.max(integrator.getAbsoluteAccuracy(),
                FastMath.abs(expected * integrator.getRelativeAccuracy()));
        result = integrator.integrate(10000, f, min, max);
        assertEquals(expected, result, tolerance);
    }

    @Test
    void testQuinticFunction() {
        UnivariateFunction f = new QuinticFunction();
        UnivariateIntegrator integrator =
                new IterativeLegendreGaussIntegrator(3,
                                                     BaseAbstractUnivariateIntegrator.DEFAULT_RELATIVE_ACCURACY,
                                                     BaseAbstractUnivariateIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                                                     BaseAbstractUnivariateIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                                                     64);
        double min, max, expected, result;

        min = 0; max = 1; expected = -1.0/48;
        result = integrator.integrate(10000, f, min, max);
        assertEquals(expected, result, 1.0e-16);

        min = 0; max = 0.5; expected = 11.0/768;
        result = integrator.integrate(10000, f, min, max);
        assertEquals(expected, result, 1.0e-16);

        min = -1; max = 4; expected = 2048/3.0 - 78 + 1.0/48;
        result = integrator.integrate(10000, f, min, max);
        assertEquals(expected, result, 1.0e-16);
    }

    @Test
    void testExactIntegration() {
        Random random = new Random(86343623467878363l);
        for (int n = 2; n < 6; ++n) {
            IterativeLegendreGaussIntegrator integrator =
                new IterativeLegendreGaussIntegrator(n,
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
                    double result    = integrator.integrate(10000, p, -5.0, 15.0);
                    double reference = exactIntegration(p, -5.0, 15.0);
                    assertEquals(reference, result, 1.0e-12 * (1.0 + FastMath.abs(reference)), n + " " + degree + " " + i);
                }
            }

        }
    }

    // Cf. MATH-995
    @Test
    void testNormalDistributionWithLargeSigma() {
        final double sigma = 1000;
        final double mean = 0;
        final double factor = 1 / (sigma * FastMath.sqrt(2 * FastMath.PI));
        final UnivariateFunction normal = new Gaussian(factor, mean, sigma);

        final double tol = 1e-2;
        final IterativeLegendreGaussIntegrator integrator =
            new IterativeLegendreGaussIntegrator(5, tol, tol);

        final double a = -5000;
        final double b = 5000;
        final double s = integrator.integrate(50, normal, a, b);
        assertEquals(1, s, 1e-5);
    }

    @Test
    void testIssue464() {
        final double value = 0.2;
        UnivariateFunction f = x -> (x >= 0 && x <= 5) ? value : 0.0;
        IterativeLegendreGaussIntegrator gauss
            = new IterativeLegendreGaussIntegrator(5, 3, 100);

        // due to the discontinuity, integration implies *many* calls
        double maxX = 0.32462367623786328;
        assertEquals(maxX * value, gauss.integrate(Integer.MAX_VALUE, f, -10, maxX), 1.0e-7);
        assertTrue(gauss.getEvaluations() > 37000000);
        assertTrue(gauss.getIterations() < 30);

        // setting up limits prevents such large number of calls
        try {
            gauss.integrate(1000, f, -10, maxX);
            fail("expected MathIllegalStateException");
        } catch (MathIllegalStateException tmee) {
            // expected
            assertEquals(LocalizedCoreFormats.MAX_COUNT_EXCEEDED, tmee.getSpecifier());
            assertEquals(1000, ((Integer) tmee.getParts()[0]).intValue());
        }

        // integrating on the two sides should be simpler
        double sum1 = gauss.integrate(1000, f, -10, 0);
        int eval1   = gauss.getEvaluations();
        double sum2 = gauss.integrate(1000, f, 0, maxX);
        int eval2   = gauss.getEvaluations();
        assertEquals(maxX * value, sum1 + sum2, 1.0e-7);
        assertTrue(eval1 + eval2 < 200);

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
