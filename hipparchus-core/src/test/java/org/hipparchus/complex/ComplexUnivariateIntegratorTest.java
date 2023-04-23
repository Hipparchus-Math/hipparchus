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
package org.hipparchus.complex;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.analysis.integration.IterativeLegendreGaussIntegrator;
import org.hipparchus.analysis.polynomials.FieldPolynomialFunction;
import org.hipparchus.util.MathUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;


public class ComplexUnivariateIntegratorTest {

    private ComplexUnivariateIntegrator integrator;

    @Test
    public void testZero() {
        final Complex start = new Complex(-1.75,   4.0);
        final Complex end   = new Complex( 1.5,  -12.0);
        UnitTestUtils.assertEquals(Complex.ZERO,
                                   integrator.integrate(1000, z -> Complex.ZERO, start, end),
                                   1.0e-15);
    }

    @Test
    public void testIdentity() {
        final Complex end = new Complex( 1.5, -12.0);
        UnitTestUtils.assertEquals(end.multiply(end).multiply(0.5),
                                   integrator.integrate(1000, z -> z, Complex.ZERO, end),
                                   1.0e-15);
    }

    @Test
    public void testPolynomialStraightPath() {
        final FieldPolynomialFunction<Complex> polynomial =
                        new FieldPolynomialFunction<>(new Complex[] {
                            new Complex(1.25, 2.0), new Complex(-3.25, 0.125), new Complex(0.0, 3.0)
                        });
        final Complex start = new Complex(-1.75,   4.0);
        final Complex end   = new Complex( 1.5,  -12.0);
        UnitTestUtils.assertEquals(polynomial.integrate(start, end),
                                   integrator.integrate(1000, polynomial, start, end),
                                   1.0e-15);
    }

    @Test
    public void testPolynomialPolylinePath() {
        final FieldPolynomialFunction<Complex> polynomial =
                        new FieldPolynomialFunction<>(new Complex[] {
                            new Complex(1.25, 2.0), new Complex(-3.25, 0.125), new Complex(0.0, 3.0)
                        });
        final Complex z0 = new Complex(-1.75,  4.0);
        final Complex z1 = new Complex( 1.00,  3.0);
        final Complex z2 = new Complex( 6.00,  0.5);
        final Complex z3 = new Complex( 6.00, -6.5);
        final Complex z4 = new Complex( 1.5, -12.0);
        UnitTestUtils.assertEquals(polynomial.integrate(z0, z4),
                                   integrator.integrate(1000, polynomial, z0, z1, z2, z3, z4),
                                   1.0e-15);
    }

    @Test
    public void testAroundPole() {
        final Complex pole = new Complex(-2.0, -1.0);
        final CalculusFieldUnivariateFunction<Complex> f = z -> z.subtract(pole).reciprocal();
        final Complex z0 = new Complex( 1,  0);
        final Complex z1 = new Complex(-1,  2);
        final Complex z2 = new Complex(-3,  2);
        final Complex z3 = new Complex(-5,  0);
        final Complex z4 = new Complex(-5, -2);
        final Complex z5 = new Complex(-4, -4);
        final Complex z6 = new Complex(-1, -4);
        final Complex z7 = new Complex( 1, -2);
        final Complex z8 = new Complex( 1,  0);
        UnitTestUtils.assertEquals(new Complex(0.0, MathUtils.TWO_PI),
                                   integrator.integrate(1000, f, z0, z1, z2, z3, z4, z5, z6, z7, z8),
                                   1.0e-15);
    }

    @Test
    public void testAroundRoot() {
        final Complex pole = new Complex(-2.0, -1.0);
        final CalculusFieldUnivariateFunction<Complex> f = z -> z.subtract(pole);
        final Complex z0 = new Complex( 1,  0);
        final Complex z1 = new Complex(-1,  2);
        final Complex z2 = new Complex(-3,  2);
        final Complex z3 = new Complex(-5,  0);
        final Complex z4 = new Complex(-5, -2);
        final Complex z5 = new Complex(-4, -4);
        final Complex z6 = new Complex(-1, -4);
        final Complex z7 = new Complex( 1, -2);
        final Complex z8 = new Complex( 1,  0);
        UnitTestUtils.assertEquals(Complex.ZERO,
                                   integrator.integrate(1000, f, z0, z1, z2, z3, z4, z5, z6, z7, z8),
                                   1.0e-15);
    }

    @Before
    public void setUp() {
        integrator = new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                          1.0e-12,
                                                                                          1.0e-12));
    }

    @After
    public void tearDown() {
        integrator = null;
    }

}
