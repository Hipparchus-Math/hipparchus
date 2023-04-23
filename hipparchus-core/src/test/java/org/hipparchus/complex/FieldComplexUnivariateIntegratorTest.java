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
import org.hipparchus.analysis.integration.IterativeLegendreFieldGaussIntegrator;
import org.hipparchus.analysis.polynomials.FieldPolynomialFunction;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.MathUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;


public class FieldComplexUnivariateIntegratorTest {

    private FieldComplexUnivariateIntegrator<Binary64> integrator;
    private FieldComplex<Binary64> zero = FieldComplex.getZero(Binary64Field.getInstance());

    private FieldComplex<Binary64> buildComplex(final double r, final double i) {
        return new FieldComplex<>(new Binary64(r), new Binary64(i));
    }

    @Test
    public void testZero() {
        final FieldComplex<Binary64> start = buildComplex(-1.75,   4.0);
        final FieldComplex<Binary64> end   = buildComplex( 1.5,  -12.0);
        UnitTestUtils.assertEquals(zero,
                                   integrator.integrate(1000, z -> zero, start, end),
                                   1.0e-15);
    }

    @Test
    public void testIdentity() {
        final FieldComplex<Binary64> end = buildComplex( 1.5, -12.0);
        UnitTestUtils.assertEquals(end.multiply(end).multiply(0.5),
                                   integrator.integrate(1000, z -> z, zero, end),
                                   1.0e-15);
    }

    @Test
    public void testPolynomialStraightPath() {
        @SuppressWarnings("unchecked")
        final FieldPolynomialFunction<FieldComplex<Binary64>> polynomial =
                        new FieldPolynomialFunction<>(new FieldComplex[] {
                            buildComplex(1.25, 2.0), buildComplex(-3.25, 0.125), buildComplex(0.0, 3.0)
                        });
        final FieldComplex<Binary64> start = buildComplex(-1.75,   4.0);
        final FieldComplex<Binary64> end   = buildComplex( 1.5,  -12.0);
        UnitTestUtils.assertEquals(polynomial.integrate(start, end),
                                   integrator.integrate(1000, polynomial, start, end),
                                   1.0e-15);
    }

    @SuppressWarnings("unchecked")
    @Test
    public void testPolynomialPolylinePath() {
        final FieldPolynomialFunction<FieldComplex<Binary64>> polynomial =
                        new FieldPolynomialFunction<>(new FieldComplex[] {
                            buildComplex(1.25, 2.0), buildComplex(-3.25, 0.125), buildComplex(0.0, 3.0)
                        });
        final FieldComplex<Binary64> z0 = buildComplex(-1.75,  4.0);
        final FieldComplex<Binary64> z1 = buildComplex( 1.00,  3.0);
        final FieldComplex<Binary64> z2 = buildComplex( 6.00,  0.5);
        final FieldComplex<Binary64> z3 = buildComplex( 6.00, -6.5);
        final FieldComplex<Binary64> z4 = buildComplex( 1.5, -12.0);
        UnitTestUtils.assertEquals(polynomial.integrate(z0, z4),
                                   integrator.integrate(1000, polynomial, z0, z1, z2, z3, z4),
                                   1.0e-15);
    }

    @SuppressWarnings("unchecked")
    @Test
    public void testAroundPole() {
        final FieldComplex<Binary64> pole = buildComplex(-2.0, -1.0);
        final CalculusFieldUnivariateFunction<FieldComplex<Binary64>> f = z -> z.subtract(pole).reciprocal();
        final FieldComplex<Binary64> z0 = buildComplex( 1,  0);
        final FieldComplex<Binary64> z1 = buildComplex(-1,  2);
        final FieldComplex<Binary64> z2 = buildComplex(-3,  2);
        final FieldComplex<Binary64> z3 = buildComplex(-5,  0);
        final FieldComplex<Binary64> z4 = buildComplex(-5, -2);
        final FieldComplex<Binary64> z5 = buildComplex(-4, -4);
        final FieldComplex<Binary64> z6 = buildComplex(-1, -4);
        final FieldComplex<Binary64> z7 = buildComplex( 1, -2);
        final FieldComplex<Binary64> z8 = buildComplex( 1,  0);
        UnitTestUtils.assertEquals(buildComplex(0.0, MathUtils.TWO_PI),
                                   integrator.integrate(1000, f, z0, z1, z2, z3, z4, z5, z6, z7, z8),
                                   1.0e-15);
    }

    @SuppressWarnings("unchecked")
    @Test
    public void testAroundRoot() {
        final FieldComplex<Binary64> pole = buildComplex(-2.0, -1.0);
        final CalculusFieldUnivariateFunction<FieldComplex<Binary64>> f = z -> z.subtract(pole);
        final FieldComplex<Binary64> z0 = buildComplex( 1,  0);
        final FieldComplex<Binary64> z1 = buildComplex(-1,  2);
        final FieldComplex<Binary64> z2 = buildComplex(-3,  2);
        final FieldComplex<Binary64> z3 = buildComplex(-5,  0);
        final FieldComplex<Binary64> z4 = buildComplex(-5, -2);
        final FieldComplex<Binary64> z5 = buildComplex(-4, -4);
        final FieldComplex<Binary64> z6 = buildComplex(-1, -4);
        final FieldComplex<Binary64> z7 = buildComplex( 1, -2);
        final FieldComplex<Binary64> z8 = buildComplex( 1,  0);
        UnitTestUtils.assertEquals(zero,
                                   integrator.integrate(1000, f, z0, z1, z2, z3, z4, z5, z6, z7, z8),
                                   1.0e-15);
    }

    @Before
    public void setUp() {
        integrator = new FieldComplexUnivariateIntegrator<>(new IterativeLegendreFieldGaussIntegrator<>(Binary64Field.getInstance(),
                                                                                                        24,
                                                                                                        1.0e-12,
                                                                                                        1.0e-12));
    }

    @After
    public void tearDown() {
        integrator = null;
    }

}
