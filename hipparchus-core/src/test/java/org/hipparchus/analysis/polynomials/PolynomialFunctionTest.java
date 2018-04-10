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
package org.hipparchus.analysis.polynomials;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.util.Decimal64;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * Tests the PolynomialFunction implementation of a UnivariateFunction.
 *
 */
public final class PolynomialFunctionTest {
    /** Error tolerance for tests */
    protected double tolerance = 1e-12;

    /**
     * tests the value of a constant polynomial.
     *
     * <p>value of this is 2.5 everywhere.</p>
     */
    @Test
    public void testConstants() {
        double[] c = { 2.5 };
        PolynomialFunction f = new PolynomialFunction(c);

        // verify that we are equal to c[0] at several (nonsymmetric) places
        Assert.assertEquals(f.value(0), c[0], tolerance);
        Assert.assertEquals(f.value(-1), c[0], tolerance);
        Assert.assertEquals(f.value(-123.5), c[0], tolerance);
        Assert.assertEquals(f.value(3), c[0], tolerance);
        Assert.assertEquals(f.value(456.89), c[0], tolerance);

        Assert.assertEquals(f.degree(), 0);
        Assert.assertEquals(f.polynomialDerivative().value(0), 0, tolerance);

        Assert.assertEquals(f.polynomialDerivative().polynomialDerivative().value(0), 0, tolerance);
    }

    /**
     * tests the value of a linear polynomial.
     *
     * <p>This will test the function f(x) = 3*x - 1.5</p>
     * <p>This will have the values
     *  <tt>f(0) = -1.5, f(-1) = -4.5, f(-2.5) = -9,
     *      f(0.5) = 0, f(1.5) = 3</tt> and {@code f(3) = 7.5}
     * </p>
     */
    @Test
    public void testLinear() {
        double[] c = { -1.5, 3 };
        PolynomialFunction f = new PolynomialFunction(c);

        // verify that we are equal to c[0] when x=0
        Assert.assertEquals(f.value(new Decimal64(0)).getReal(), c[0], tolerance);

        // now check a few other places
        Assert.assertEquals(-4.5, f.value(new Decimal64(-1)).getReal(), tolerance);
        Assert.assertEquals(-9, f.value(new Decimal64(-2.5)).getReal(), tolerance);
        Assert.assertEquals(0, f.value(new Decimal64(0.5)).getReal(), tolerance);
        Assert.assertEquals(3, f.value(new Decimal64(1.5)).getReal(), tolerance);
        Assert.assertEquals(7.5, f.value(new Decimal64(3)).getReal(), tolerance);

        Assert.assertEquals(f.degree(), 1);

        Assert.assertEquals(f.polynomialDerivative().polynomialDerivative().value(0), 0, tolerance);
    }

    /**
     * Tests a second order polynomial.
     * <p> This will test the function f(x) = 2x^2 - 3x -2 = (2x+1)(x-2)</p>
     */
    @Test
    public void testQuadratic() {
        double[] c = { -2, -3, 2 };
        PolynomialFunction f = new PolynomialFunction(c);

        // verify that we are equal to c[0] when x=0
        Assert.assertEquals(f.value(0), c[0], tolerance);

        // now check a few other places
        Assert.assertEquals(0, f.value(-0.5), tolerance);
        Assert.assertEquals(0, f.value(2), tolerance);
        Assert.assertEquals(-2, f.value(1.5), tolerance);
        Assert.assertEquals(7, f.value(-1.5), tolerance);
        Assert.assertEquals(265.5312, f.value(12.34), tolerance);
    }

    /**
     * This will test the quintic function
     *   f(x) = x^2(x-5)(x+3)(x-1) = x^5 - 3x^4 -13x^3 + 15x^2</p>
     */
    @Test
    public void testQuintic() {
        double[] c = { 0, 0, 15, -13, -3, 1 };
        PolynomialFunction f = new PolynomialFunction(c);

        // verify that we are equal to c[0] when x=0
        Assert.assertEquals(f.value(0), c[0], tolerance);

        // now check a few other places
        Assert.assertEquals(0, f.value(5), tolerance);
        Assert.assertEquals(0, f.value(1), tolerance);
        Assert.assertEquals(0, f.value(-3), tolerance);
        Assert.assertEquals(54.84375, f.value(-1.5), tolerance);
        Assert.assertEquals(-8.06637, f.value(1.3), tolerance);

        Assert.assertEquals(f.degree(), 5);
    }

    /**
     * tests the firstDerivative function by comparison
     *
     * <p>This will test the functions
     * {@code f(x) = x^3 - 2x^2 + 6x + 3, g(x) = 3x^2 - 4x + 6}
     * and {@code h(x) = 6x - 4}
     */
    @Test
    public void testfirstDerivativeComparison() {
        double[] f_coeff = { 3, 6, -2, 1 };
        double[] g_coeff = { 6, -4, 3 };
        double[] h_coeff = { -4, 6 };

        PolynomialFunction f = new PolynomialFunction(f_coeff);
        PolynomialFunction g = new PolynomialFunction(g_coeff);
        PolynomialFunction h = new PolynomialFunction(h_coeff);

        // compare f' = g
        Assert.assertEquals(f.polynomialDerivative().value(0), g.value(0), tolerance);
        Assert.assertEquals(f.polynomialDerivative().value(1), g.value(1), tolerance);
        Assert.assertEquals(f.polynomialDerivative().value(100), g.value(100), tolerance);
        Assert.assertEquals(f.polynomialDerivative().value(4.1), g.value(4.1), tolerance);
        Assert.assertEquals(f.polynomialDerivative().value(-3.25), g.value(-3.25), tolerance);

        // compare g' = h
        Assert.assertEquals(g.polynomialDerivative().value(FastMath.PI), h.value(FastMath.PI), tolerance);
        Assert.assertEquals(g.polynomialDerivative().value(FastMath.E),  h.value(FastMath.E),  tolerance);
    }

    @Test
    public void testString() {
        PolynomialFunction p = new PolynomialFunction(new double[] { -5, 3, 1 });
        checkPolynomial(p, "-5 + 3 x + x^2");
        checkPolynomial(new PolynomialFunction(new double[] { 0, -2, 3 }),
                        "-2 x + 3 x^2");
        checkPolynomial(new PolynomialFunction(new double[] { 1, -2, 3 }),
                      "1 - 2 x + 3 x^2");
        checkPolynomial(new PolynomialFunction(new double[] { 0,  2, 3 }),
                       "2 x + 3 x^2");
        checkPolynomial(new PolynomialFunction(new double[] { 1,  2, 3 }),
                     "1 + 2 x + 3 x^2");
        checkPolynomial(new PolynomialFunction(new double[] { 1,  0, 3 }),
                     "1 + 3 x^2");
        checkPolynomial(new PolynomialFunction(new double[] { 0 }),
                     "0");
    }

    @Test
    public void testAddition() {
        PolynomialFunction p1 = new PolynomialFunction(new double[] { -2, 1 });
        PolynomialFunction p2 = new PolynomialFunction(new double[] { 2, -1, 0 });
        checkNullPolynomial(p1.add(p2));

        p2 = p1.add(p1);
        checkPolynomial(p2, "-4 + 2 x");

        p1 = new PolynomialFunction(new double[] { 1, -4, 2 });
        p2 = new PolynomialFunction(new double[] { -1, 3, -2 });
        p1 = p1.add(p2);
        Assert.assertEquals(1, p1.degree());
        checkPolynomial(p1, "-x");
    }

    @Test
    public void testSubtraction() {
        PolynomialFunction p1 = new PolynomialFunction(new double[] { -2, 1 });
        checkNullPolynomial(p1.subtract(p1));

        PolynomialFunction p2 = new PolynomialFunction(new double[] { -2, 6 });
        p2 = p2.subtract(p1);
        checkPolynomial(p2, "5 x");

        p1 = new PolynomialFunction(new double[] { 1, -4, 2 });
        p2 = new PolynomialFunction(new double[] { -1, 3, 2 });
        p1 = p1.subtract(p2);
        Assert.assertEquals(1, p1.degree());
        checkPolynomial(p1, "2 - 7 x");
    }

    @Test
    public void testMultiplication() {
        PolynomialFunction p1 = new PolynomialFunction(new double[] { -3, 2 });
        PolynomialFunction p2 = new PolynomialFunction(new double[] { 3, 2, 1 });
        checkPolynomial(p1.multiply(p2), "-9 + x^2 + 2 x^3");

        p1 = new PolynomialFunction(new double[] { 0, 1 });
        p2 = p1;
        for (int i = 2; i < 10; ++i) {
            p2 = p2.multiply(p1);
            checkPolynomial(p2, "x^" + i);
        }
    }

    @Test
    public void testSerial() {
        PolynomialFunction p2 = new PolynomialFunction(new double[] { 3, 2, 1 });
        Assert.assertEquals(p2, UnitTestUtils.serializeAndRecover(p2));
    }

    /**
     * tests the firstDerivative function by comparison
     *
     * <p>This will test the functions
     * {@code f(x) = x^3 - 2x^2 + 6x + 3, g(x) = 3x^2 - 4x + 6}
     * and {@code h(x) = 6x - 4}
     */
    @Test
    public void testMath341() {
        double[] f_coeff = { 3, 6, -2, 1 };
        double[] g_coeff = { 6, -4, 3 };
        double[] h_coeff = { -4, 6 };

        PolynomialFunction f = new PolynomialFunction(f_coeff);
        PolynomialFunction g = new PolynomialFunction(g_coeff);
        PolynomialFunction h = new PolynomialFunction(h_coeff);

        // compare f' = g
        Assert.assertEquals(f.polynomialDerivative().value(0), g.value(0), tolerance);
        Assert.assertEquals(f.polynomialDerivative().value(1), g.value(1), tolerance);
        Assert.assertEquals(f.polynomialDerivative().value(100), g.value(100), tolerance);
        Assert.assertEquals(f.polynomialDerivative().value(4.1), g.value(4.1), tolerance);
        Assert.assertEquals(f.polynomialDerivative().value(-3.25), g.value(-3.25), tolerance);

        // compare g' = h
        Assert.assertEquals(g.polynomialDerivative().value(FastMath.PI), h.value(FastMath.PI), tolerance);
        Assert.assertEquals(g.polynomialDerivative().value(FastMath.E),  h.value(FastMath.E),  tolerance);
    }

    @Test
    public void testAntiDerivative() {
        // 1 + 2x + 3x^2
        final double[] coeff = {1, 2, 3};
        final PolynomialFunction p = new PolynomialFunction(coeff);
        // x + x^2 + x^3
        final double[] aCoeff = {0, 1, 1, 1};
        Assert.assertArrayEquals(aCoeff, p.antiDerivative().getCoefficients(), Double.MIN_VALUE);
    }

    @Test
    public void testAntiDerivativeConstant() {
        final double[] coeff = {2};
        final PolynomialFunction p = new PolynomialFunction(coeff);
        final double[] aCoeff = {0, 2};
        Assert.assertArrayEquals(aCoeff, p.antiDerivative().getCoefficients(), Double.MIN_VALUE);
    }

    @Test
    public void testAntiDerivativeZero() {
        final double[] coeff = {0};
        final PolynomialFunction p = new PolynomialFunction(coeff);
        final double[] aCoeff = {0};
        Assert.assertArrayEquals(aCoeff, p.antiDerivative().getCoefficients(), Double.MIN_VALUE);
    }

    @Test
    public void testAntiDerivativeRandom() {
        final RandomDataGenerator ran = new RandomDataGenerator(1000);
        double[] coeff = null;
        PolynomialFunction p = null;
        int d = 0;
        for (int i = 0; i < 20; i++) {
            d = ran.nextInt(1, 50);
            coeff = new double[d];
            for (int j = 0; j < d; j++) {
                coeff[j] = ran.nextUniform(-100, 1000);
            }
            p = new PolynomialFunction(coeff);
            checkInverseDifferentiation(p);
        }
    }

    @Test
    public void testIntegrate() {
        // -x^2
        final double[] coeff = {0, 0, -1};
        final PolynomialFunction p = new PolynomialFunction(coeff);
        Assert.assertEquals(-2d/3d, p.integrate(-1, 1),Double.MIN_VALUE);

        // x(x-1)(x+1) - should integrate to 0 over [-1,1]
        final PolynomialFunction p2 = new PolynomialFunction(new double[] {0, 1}).
                multiply(new PolynomialFunction(new double[]{-1,1})).
                         multiply(new PolynomialFunction(new double[] {1, 1}));
        Assert.assertEquals(0, p2.integrate(-1, 1), Double.MIN_VALUE);
    }

    @Test(expected = MathIllegalArgumentException.class)
    public void testIntegrateInfiniteBounds() {
        final PolynomialFunction p = new PolynomialFunction(new double[]{1});
        p.integrate(0, Double.POSITIVE_INFINITY);
    }

    @Test(expected = MathIllegalArgumentException.class)
    public void testIntegrateBadInterval() {
        final PolynomialFunction p = new PolynomialFunction(new double[]{1});
        p.integrate(0, -1);
    }

    public void checkPolynomial(PolynomialFunction p, String reference) {
        Assert.assertEquals(reference, p.toString());
    }

    private void checkInverseDifferentiation(PolynomialFunction p) {
        Assert.assertArrayEquals(p.getCoefficients(),
                                 p.antiDerivative().polynomialDerivative().getCoefficients(),
                                 1e-12);
    }

    private void checkNullPolynomial(PolynomialFunction p) {
        for (double coefficient : p.getCoefficients()) {
            Assert.assertEquals(0, coefficient, 1e-15);
        }
    }
}
