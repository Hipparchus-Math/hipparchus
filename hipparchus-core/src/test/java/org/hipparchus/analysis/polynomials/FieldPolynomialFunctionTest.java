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
package org.hipparchus.analysis.polynomials;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.analysis.differentiation.Gradient;
import org.hipparchus.analysis.differentiation.UnivariateDerivative1;
import org.hipparchus.analysis.differentiation.UnivariateDerivative2;
import org.hipparchus.complex.Complex;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Tests the FieldPolynomialFunction implementation of a UnivariateFunction.
 *
 */
final class FieldPolynomialFunctionTest {
    /** Error tolerance for tests */
    protected double tolerance = 1e-12;

    /**
     * tests the value of a constant polynomial.
     *
     * <p>value of this is 2.5 everywhere.</p>
     */
    @Test
    void testConstants() {
        double c0 = 2.5;
        FieldPolynomialFunction<Binary64> f = buildD64(c0);

        // verify that we are equal to c[0] at several (nonsymmetric) places
        assertEquals(c0, f.value(0).getReal(), tolerance);
        assertEquals(c0, f.value(-1).getReal(), tolerance);
        assertEquals(c0, f.value(-123.5).getReal(), tolerance);
        assertEquals(c0, f.value(3).getReal(), tolerance);
        assertEquals(c0, f.value(new Binary64(456.89)).getReal(), tolerance);

        assertEquals(0, f.degree());
        assertEquals(0, f.polynomialDerivative().value(0).getReal(), tolerance);

        assertEquals(0, f.polynomialDerivative().polynomialDerivative().value(0).getReal(), tolerance);
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
    void testLinear() {
       FieldPolynomialFunction<Binary64> f = buildD64(-1.5, 3);

        // verify that we are equal to c[0] when x=0
        assertEquals(-1.5, f.value(new Binary64(0)).getReal(), tolerance);

        // now check a few other places
        assertEquals(-4.5, f.value(new Binary64(-1)).getReal(), tolerance);
        assertEquals(-9, f.value(new Binary64(-2.5)).getReal(), tolerance);
        assertEquals(0, f.value(new Binary64(0.5)).getReal(), tolerance);
        assertEquals(3, f.value(new Binary64(1.5)).getReal(), tolerance);
        assertEquals(7.5, f.value(new Binary64(3)).getReal(), tolerance);

        assertEquals(1, f.degree());

        assertEquals(0, f.polynomialDerivative().polynomialDerivative().value(0).getReal(), tolerance);
    }

    /**
     * Tests a second order polynomial.
     * <p> This will test the function f(x) = 2x^2 - 3x -2 = (2x+1)(x-2)</p>
     */
    @Test
    void testQuadratic() {
        FieldPolynomialFunction<Binary64> f = buildD64(-2, -3, 2);

        // verify that we are equal to c[0] when x=0
        assertEquals(-2, f.value(0).getReal(), tolerance);

        // now check a few other places
        assertEquals(0, f.value(-0.5).getReal(), tolerance);
        assertEquals(0, f.value(2).getReal(), tolerance);
        assertEquals(-2, f.value(1.5).getReal(), tolerance);
        assertEquals(7, f.value(-1.5).getReal(), tolerance);
        assertEquals(265.5312, f.value(12.34).getReal(), tolerance);
    }

    /**
     * This will test the quintic function
     *   f(x) = x^2(x-5)(x+3)(x-1) = x^5 - 3x^4 -13x^3 + 15x^2</p>
     */
    @Test
    void testQuintic() {
        FieldPolynomialFunction<Binary64> f = buildD64(0, 0, 15, -13, -3, 1);

        // verify that we are equal to c[0] when x=0
        assertEquals(0, f.value(0).getReal(), tolerance);

        // now check a few other places
        assertEquals(0, f.value(5).getReal(), tolerance);
        assertEquals(0, f.value(1).getReal(), tolerance);
        assertEquals(0, f.value(-3).getReal(), tolerance);
        assertEquals(54.84375, f.value(-1.5).getReal(), tolerance);
        assertEquals(-8.06637, f.value(1.3).getReal(), tolerance);

        assertEquals(5, f.degree());
    }

    /**
     * tests the firstDerivative function by comparison
     *
     * <p>This will test the functions
     * {@code f(x) = x^3 - 2x^2 + 6x + 3, g(x) = 3x^2 - 4x + 6}
     * and {@code h(x) = 6x - 4}
     */
    @Test
    void testfirstDerivativeComparison() {
        double[] f_coeff = { 3, 6, -2, 1 };
        double[] g_coeff = { 6, -4, 3 };
        double[] h_coeff = { -4, 6 };

        FieldPolynomialFunction<Binary64> f = buildD64(f_coeff);
        FieldPolynomialFunction<Binary64> g = buildD64(g_coeff);
        FieldPolynomialFunction<Binary64> h = buildD64(h_coeff);

        // compare f' = g
        assertEquals(f.polynomialDerivative().value(0).getReal(), g.value(0).getReal(), tolerance);
        assertEquals(f.polynomialDerivative().value(1).getReal(), g.value(1).getReal(), tolerance);
        assertEquals(f.polynomialDerivative().value(100).getReal(), g.value(100).getReal(), tolerance);
        assertEquals(f.polynomialDerivative().value(4.1).getReal(), g.value(4.1).getReal(), tolerance);
        assertEquals(f.polynomialDerivative().value(-3.25).getReal(), g.value(-3.25).getReal(), tolerance);

        // compare g' = h
        assertEquals(g.polynomialDerivative().value(FastMath.PI).getReal(), h.value(FastMath.PI).getReal(), tolerance);
        assertEquals(g.polynomialDerivative().value(FastMath.E).getReal(),  h.value(FastMath.E).getReal(),  tolerance);
    }

    @Test
    void testAddition() {
        FieldPolynomialFunction<Binary64> p1 = buildD64( -2, 1 );
        FieldPolynomialFunction<Binary64> p2 = buildD64( 2, -1, 0 );
        checkNullPolynomial(p1.add(p2));

        p2 = p1.add(p1);
        checkCoeffs(Double.MIN_VALUE, p2, -4, 2);

        p1 = buildD64( 1, -4, 2 );
        p2 = buildD64( -1, 3, -2 );
        p1 = p1.add(p2);
        assertEquals(1, p1.degree());
        checkCoeffs(Double.MIN_VALUE, p1, 0, -1);
    }

    @Test
    void testSubtraction() {
        FieldPolynomialFunction<Binary64> p1 = buildD64( -2, 1 );
        checkNullPolynomial(p1.subtract(p1));

        FieldPolynomialFunction<Binary64> p2 = buildD64( -2, 6 );
        p2 = p2.subtract(p1);
        checkCoeffs(Double.MIN_VALUE, p2, 0, 5);

        p1 = buildD64( 1, -4, 2 );
        p2 = buildD64( -1, 3, 2 );
        p1 = p1.subtract(p2);
        assertEquals(1, p1.degree());
        checkCoeffs(Double.MIN_VALUE, p1, 2, -7);
    }

    @Test
    void testMultiplication() {
        FieldPolynomialFunction<Binary64> p1 = buildD64( -3, 2 );
        FieldPolynomialFunction<Binary64> p2 = buildD64( 3, 2, 1 );
        checkCoeffs(Double.MIN_VALUE, p1.multiply(p2), -9, 0, 1, 2);

        p1 = buildD64( 0, 1 );
        p2 = p1;
        for (int i = 2; i < 10; ++i) {
            p2 = p2.multiply(p1);
            double[] c = new double[i + 1];
            c[i] = 1;
            checkCoeffs(Double.MIN_VALUE, p2, c);
        }
    }

    /**
     * tests the firstDerivative function by comparison
     *
     * <p>This will test the functions
     * {@code f(x) = x^3 - 2x^2 + 6x + 3, g(x) = 3x^2 - 4x + 6}
     * and {@code h(x) = 6x - 4}
     */
    @Test
    void testMath341() {
        double[] f_coeff = { 3, 6, -2, 1 };
        double[] g_coeff = { 6, -4, 3 };
        double[] h_coeff = { -4, 6 };

        FieldPolynomialFunction<Binary64> f = buildD64(f_coeff);
        FieldPolynomialFunction<Binary64> g = buildD64(g_coeff);
        FieldPolynomialFunction<Binary64> h = buildD64(h_coeff);

        // compare f' = g
        assertEquals(f.polynomialDerivative().value(0).getReal(), g.value(0).getReal(), tolerance);
        assertEquals(f.polynomialDerivative().value(1).getReal(), g.value(1).getReal(), tolerance);
        assertEquals(f.polynomialDerivative().value(100).getReal(), g.value(100).getReal(), tolerance);
        assertEquals(f.polynomialDerivative().value(4.1).getReal(), g.value(4.1).getReal(), tolerance);
        assertEquals(f.polynomialDerivative().value(-3.25).getReal(), g.value(-3.25).getReal(), tolerance);

        // compare g' = h
        assertEquals(g.polynomialDerivative().value(FastMath.PI).getReal(), h.value(FastMath.PI).getReal(), tolerance);
        assertEquals(g.polynomialDerivative().value(FastMath.E).getReal(),  h.value(FastMath.E).getReal(),  tolerance);
    }

    @Test
    void testAntiDerivative() {
        // 1 + 2x + 3x^2
        final double[] coeff = {1, 2, 3};
        final FieldPolynomialFunction<Binary64> p = buildD64(coeff);
        // x + x^2 + x^3
        checkCoeffs(Double.MIN_VALUE, p.antiDerivative(), 0, 1, 1, 1);
    }

    @Test
    void testAntiDerivativeConstant() {
        final double[] coeff = {2};
        final FieldPolynomialFunction<Binary64> p = buildD64(coeff);
        checkCoeffs(Double.MIN_VALUE, p.antiDerivative(), 0, 2);
    }

    @Test
    void testAntiDerivativeZero() {
        final double[] coeff = {0};
        final FieldPolynomialFunction<Binary64> p = buildD64(coeff);
        checkCoeffs(Double.MIN_VALUE, p.antiDerivative(), 0);
    }

    @Test
    void testAntiDerivativeRandom() {
        final RandomDataGenerator ran = new RandomDataGenerator(1000);
        double[] coeff = null;
        FieldPolynomialFunction<Binary64> p = null;
        int d = 0;
        for (int i = 0; i < 20; i++) {
            d = ran.nextInt(1, 50);
            coeff = new double[d];
            for (int j = 0; j < d; j++) {
                coeff[j] = ran.nextUniform(-100, 1000);
            }
            p = buildD64(coeff);
            checkInverseDifferentiation(p);
        }
    }

    @Test
    void testIntegrate() {
        // -x^2
        final double[] coeff = {0, 0, -1};
        final FieldPolynomialFunction<Binary64> p = buildD64(coeff);
        assertEquals(-2d/3d, p.integrate(-1, 1).getReal(),Double.MIN_VALUE);

        // x(x-1)(x+1) - should integrate to 0 over [-1,1]
        final FieldPolynomialFunction<Binary64> p2 = buildD64(0, 1).
                                                      multiply(buildD64(-1, 1)).
                                                      multiply(buildD64(1, 1));
        assertEquals(0, p2.integrate(-1, 1).getReal(), Double.MIN_VALUE);
    }

    @Test
    void testIntegrateInfiniteBounds() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final FieldPolynomialFunction<Binary64> p = buildD64(1);
            p.integrate(0, Double.POSITIVE_INFINITY);
        });
    }

    @Test
    void testIntegrateBadInterval() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final FieldPolynomialFunction<Binary64> p = buildD64(1);
            p.integrate(0, -1);
        });
    }

    @Test
    void testIssue259WithComplex() {
        final double nonZeroParameterForQuadraticCoeff = -1.;
        final Complex[] coefficients = buildImaginaryCoefficients(1., 2., nonZeroParameterForQuadraticCoeff);
        templateIssue259DisappearingCoefficients(coefficients);
    }

    @Test
    void testIssue259WithGradient() {
        final double nonZeroParameterForQuadraticCoeff = -1.;
        final Gradient[] coefficients = buildUnivariateGradientCoefficients(1.,2., nonZeroParameterForQuadraticCoeff);
        templateIssue259DisappearingCoefficients(coefficients);
    }

    @Test
    void testIssue259WithDerivativeStructure() {
        final double nonZeroParameterForQuadraticCoeff = -1.;
        final DerivativeStructure[] coefficients = buildUnivariateDSCoefficients(1.,2., nonZeroParameterForQuadraticCoeff);
        templateIssue259DisappearingCoefficients(coefficients);
    }

    @Test
    void testIssue259WithUnivariateDerivative1() {
        final double nonZeroParameterForQuadraticCoeff = -1.;
        final UnivariateDerivative1[] coefficients = buildUnivariateDerivative1Coefficients(1.,2., nonZeroParameterForQuadraticCoeff);
        templateIssue259DisappearingCoefficients(coefficients);
    }

    @Test
    void testIssue259WithUnivariateDerivative2() {
        final double nonZeroParameterForQuadraticCoeff = -1.;
        final UnivariateDerivative2[] coefficients = buildUnivariateDerivative2Coefficients(1.,2., nonZeroParameterForQuadraticCoeff);
        templateIssue259DisappearingCoefficients(coefficients);
    }

    private <T extends CalculusFieldElement<T>> void templateIssue259DisappearingCoefficients(T[] coefficients) {
        final FieldPolynomialFunction<T> polynomialFunction = new FieldPolynomialFunction<>(coefficients);
        assertEquals(coefficients.length, polynomialFunction.getCoefficients().length);
        for (int i = 0; i < coefficients.length; i++) {
            assertEquals(coefficients[i], polynomialFunction.getCoefficients()[i]);
        }
    }

    private <T extends CalculusFieldElement<T>> void checkInverseDifferentiation(FieldPolynomialFunction<T> p) {
        final T[] c0 = p.getCoefficients();
        final T[] c1 = p.antiDerivative().polynomialDerivative().getCoefficients();
        assertEquals(c0.length, c1.length);
        for (int i = 0; i < c0.length; ++i) {
            assertEquals(c0[i].getReal(), c1[i].getReal(), 1e-12);
        }
    }

    private <T extends CalculusFieldElement<T>> void checkCoeffs(final double tolerance, final FieldPolynomialFunction<T> p,
                                                                 final double... ref) {
        final T[] c = p.getCoefficients();
        assertEquals(ref.length, c.length);
        for (int i = 0; i < ref.length; ++i) {
            assertEquals(ref[i], c[i].getReal(), tolerance);
        }
    }

    private <T extends CalculusFieldElement<T>> void checkNullPolynomial(FieldPolynomialFunction<T> p) {
        for (T coefficient : p.getCoefficients()) {
            assertEquals(0, coefficient.getReal(), 1e-15);
        }
    }

    private FieldPolynomialFunction<Binary64> buildD64(double...c) {
        Binary64[] array = new Binary64[c.length];
        for (int i = 0; i < c.length; ++i) {
            array[i] = new Binary64(c[i]);
        }
        return new FieldPolynomialFunction<>(array);
    }

    private Complex[] buildImaginaryCoefficients(double...c) {
        Complex[] array = new Complex[c.length];
        for (int i = 0; i < c.length; ++i) {
            array[i] = new Complex(0., c[i]);
        }
        return array;
    }

    private DerivativeStructure[] buildUnivariateDSCoefficients(double...c) {
        DerivativeStructure[] array = new DerivativeStructure[c.length];
        final DSFactory factory = new DSFactory(1, 1);
        for (int i = 0; i < c.length; ++i) {
            array[i] = factory.variable(0, 0.).multiply(c[i]);
        }
        return array;
    }

    private Gradient[] buildUnivariateGradientCoefficients(double...c) {
        Gradient[] array = new Gradient[c.length];
        for (int i = 0; i < c.length; ++i) {
            array[i] = Gradient.variable(1, 0, 0.).multiply(c[i]);
        }
        return array;
    }

    private UnivariateDerivative1[] buildUnivariateDerivative1Coefficients(double...c) {
        UnivariateDerivative1[] array = new UnivariateDerivative1[c.length];
        for (int i = 0; i < c.length; ++i) {
            array[i] = new UnivariateDerivative1(0., c[i]);
        }
        return array;
    }

    private UnivariateDerivative2[] buildUnivariateDerivative2Coefficients(double...c) {
        UnivariateDerivative2[] array = new UnivariateDerivative2[c.length];
        for (int i = 0; i < c.length; ++i) {
            array[i] = new UnivariateDerivative2(0., c[i], 0.);
        }
        return array;
    }

}
