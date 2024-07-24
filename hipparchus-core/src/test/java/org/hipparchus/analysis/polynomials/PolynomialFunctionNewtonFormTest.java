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

import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Binary64;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;


/**
 * Test case for Newton form of polynomial function.
 * <p>
 * The small tolerance number is used only to account for round-off errors.
 *
 */
final class PolynomialFunctionNewtonFormTest {

    /**
     * Test of polynomial for the linear function.
     */
    @Test
    void testLinearFunction() {
        PolynomialFunctionNewtonForm p;
        double[] coefficients;
        double z;
        double expected;
        double result;
        double tolerance = 1E-12;

        // p(x) = 1.5x - 4 = 2 + 1.5(x-4)
        double[] a = { 2.0, 1.5 };
        double[] c = { 4.0 };
        p = new PolynomialFunctionNewtonForm(a, c);

        z = 2.0; expected = -1.0; result = p.value(z);
        assertEquals(expected, result, tolerance);

        z = 4.5; expected = 2.75; result = p.value(z);
        assertEquals(expected, result, tolerance);

        z = 6.0; expected = 5.0; result = p.value(new Binary64(z)).getReal();
        assertEquals(expected, result, tolerance);

        assertEquals(1, p.degree());

        coefficients = p.getCoefficients();
        assertEquals(2, coefficients.length);
        assertEquals(-4.0, coefficients[0], tolerance);
        assertEquals(1.5, coefficients[1], tolerance);
    }

    /**
     * Test of polynomial for the quadratic function.
     */
    @Test
    void testQuadraticFunction() {
        PolynomialFunctionNewtonForm p;
        double[] coefficients;
        double z;
        double expected;
        double result;
        double tolerance = 1E-12;

        // p(x) = 2x^2 + 5x - 3 = 4 + 3(x-1) + 2(x-1)(x+2)
        double[] a = { 4.0, 3.0, 2.0 };
        double[] c = { 1.0, -2.0 };
        p = new PolynomialFunctionNewtonForm(a, c);

        z = 1.0; expected = 4.0; result = p.value(z);
        assertEquals(expected, result, tolerance);

        z = 2.5; expected = 22.0; result = p.value(z);
        assertEquals(expected, result, tolerance);

        z = -2.0; expected = -5.0; result = p.value(z);
        assertEquals(expected, result, tolerance);

        assertEquals(2, p.degree());

        coefficients = p.getCoefficients();
        assertEquals(3, coefficients.length);
        assertEquals(-3.0, coefficients[0], tolerance);
        assertEquals(5.0, coefficients[1], tolerance);
        assertEquals(2.0, coefficients[2], tolerance);
    }

    /**
     * Test of polynomial for the quintic function.
     */
    @Test
    void testQuinticFunction() {
        PolynomialFunctionNewtonForm p;
        double[] coefficients;
        double z;
        double expected;
        double result;
        double tolerance = 1E-12;

        // p(x) = x^5 - x^4 - 7x^3 + x^2 + 6x
        //      = 6x - 6x^2 -6x^2(x-1) + x^2(x-1)(x+1) + x^2(x-1)(x+1)(x-2)
        double[] a = { 0.0, 6.0, -6.0, -6.0, 1.0, 1.0 };
        double[] c = { 0.0, 0.0, 1.0, -1.0, 2.0 };
        p = new PolynomialFunctionNewtonForm(a, c);

        z = 0.0; expected = 0.0; result = p.value(z);
        assertEquals(expected, result, tolerance);

        z = -2.0; expected = 0.0; result = p.value(z);
        assertEquals(expected, result, tolerance);

        z = 4.0; expected = 360.0; result = p.value(z);
        assertEquals(expected, result, tolerance);

        assertEquals(5, p.degree());

        coefficients = p.getCoefficients();
        assertEquals(6, coefficients.length);
        assertEquals(0.0, coefficients[0], tolerance);
        assertEquals(6.0, coefficients[1], tolerance);
        assertEquals(1.0, coefficients[2], tolerance);
        assertEquals(-7.0, coefficients[3], tolerance);
        assertEquals(-1.0, coefficients[4], tolerance);
        assertEquals(1.0, coefficients[5], tolerance);
    }

    /**
     * Test for derivatives.
     */
    @Test
    void testDerivative() {

        // x^3 = 0 * [1] + 1 * [x] + 3 * [x(x-1)] + 1 * [x(x-1)(x-2)]
        PolynomialFunctionNewtonForm p =
                new PolynomialFunctionNewtonForm(new double[] { 0, 1, 3, 1 },
                                                 new double[] { 0, 1, 2 });

        double eps = 2.0e-14;
        DSFactory factory = new DSFactory(1, 4);
        for (double t = 0.0; t < 10.0; t += 0.1) {
            DerivativeStructure x = factory.variable(0, t);
            DerivativeStructure y = p.value(x);
            assertEquals(t * t * t,   y.getValue(),              eps * t * t * t);
            assertEquals(3.0 * t * t, y.getPartialDerivative(1), eps * 3.0 * t * t);
            assertEquals(6.0 * t,     y.getPartialDerivative(2), eps * 6.0 * t);
            assertEquals(6.0,         y.getPartialDerivative(3), eps * 6.0);
            assertEquals(0.0,         y.getPartialDerivative(4), eps);
        }

    }

    /**
     * Test of parameters for the polynomial.
     */
    @Test
    void testParameters() {

        try {
            // bad input array length
            double[] a = { 1.0 };
            double[] c = { 2.0 };
            new PolynomialFunctionNewtonForm(a, c);
            fail("Expecting MathIllegalArgumentException - bad input array length");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            // mismatch input arrays
            double[] a = { 1.0, 2.0, 3.0, 4.0 };
            double[] c = { 4.0, 3.0, 2.0, 1.0 };
            new PolynomialFunctionNewtonForm(a, c);
            fail("Expecting MathIllegalArgumentException - mismatch input arrays");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }
}
