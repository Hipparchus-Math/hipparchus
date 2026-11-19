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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.Binary64;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Tests the PolynomialSplineFunction implementation.
 *
 */
class PolynomialSplineFunctionTest {

    /** Error tolerance for tests */
    protected double tolerance = 1.0e-12;

    /**
     * Quadratic polynomials used in tests:
     *
     * x^2 + x            [-1, 0)
     * x^2 + x + 2        [0, 1)
     * x^2 + x + 4        [1, 2)
     *
     * Defined so that evaluation using PolynomialSplineFunction evaluation
     * algorithm agrees at knot point boundaries.
     */
    protected PolynomialFunction[] polynomials = {
        new PolynomialFunction(new double[] {0d, 1d, 1d}),
        new PolynomialFunction(new double[] {2d, 1d, 1d}),
        new PolynomialFunction(new double[] {4d, 1d, 1d})
    };

    /** Knot points  */
    protected double[] knots = {-1, 0, 1, 2};

    /** Derivative of test polynomials -- 2x + 1  */
    protected PolynomialFunction dp =
        new PolynomialFunction(new double[] {1d, 2d});


    @Test
    void testConstructor() {
        PolynomialSplineFunction spline =
            new PolynomialSplineFunction(knots, polynomials);
        assertArrayEquals(knots, spline.getKnots());
        assertEquals(1d, spline.getPolynomials()[0].getCoefficients()[2], 0);
        assertEquals(3, spline.getN());

        try { // too few knots
            new PolynomialSplineFunction(new double[] {0}, polynomials);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try { // too many knots
            new PolynomialSplineFunction(new double[] {0,1,2,3,4}, polynomials);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try { // knots not increasing
            new PolynomialSplineFunction(new double[] {0,1, 3, 2}, polynomials);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testValues() {
        PolynomialSplineFunction spline =
            new PolynomialSplineFunction(knots, polynomials);
        PolynomialSplineFunction dSpline = spline.polynomialSplineDerivative();

        /**
         * interior points -- spline value at x should equal p(x - knot)
         * where knot is the largest knot point less than or equal to x and p
         * is the polynomial defined over the knot segment to which x belongs.
         */
        double x = -1;
        int index = 0;
        for (int i = 0; i < 10; i++) {
           x+=0.25;
           index = findKnot(knots, x);
           assertEquals(polynomials[index].value(x - knots[index]), spline.value(x), tolerance, "spline function evaluation failed for x=" + x);
           assertEquals(dp.value(x - knots[index]), dSpline.value(x), tolerance, "spline derivative evaluation failed for x=" + x);
        }

        // knot points -- centering should zero arguments
        for (int i = 0; i < 3; i++) {
            assertEquals(polynomials[i].value(0), spline.value(knots[i]), tolerance, "spline function evaluation failed for knot=" + knots[i]);
            assertEquals(dp.value(0), dSpline.value(new Binary64(knots[i])).getReal(), tolerance, "spline function evaluation failed for knot=" + knots[i]);
        }

        try { //outside of domain -- under min
            x = spline.value(-1.5);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try { //outside of domain -- over max
            x = spline.value(2.5);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testIsValidPoint() {
        final PolynomialSplineFunction spline =
            new PolynomialSplineFunction(knots, polynomials);
        final double xMin = knots[0];
        final double xMax = knots[knots.length - 1];

        double x;

        x = xMin;
        assertTrue(spline.isValidPoint(x));
        // Ensure that no exception is thrown.
        spline.value(x);

        x = xMax;
        assertTrue(spline.isValidPoint(x));
        // Ensure that no exception is thrown.
        spline.value(x);

        final double xRange = xMax - xMin;
        x = xMin + xRange / 3.4;
        assertTrue(spline.isValidPoint(x));
        // Ensure that no exception is thrown.
        spline.value(x);

        final double small = 1e-8;
        x = xMin - small;
        assertFalse(spline.isValidPoint(x));
        // Ensure that an exception would have been thrown.
        try {
            spline.value(x);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException expected) {}
    }

    /**
     *  Do linear search to find largest knot point less than or equal to x.
     *  Implementation does binary search.
     */
     protected int findKnot(double[] knots, double x) {
         if (x < knots[0] || x >= knots[knots.length -1]) {
             throw new MathIllegalArgumentException(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE,
                                                    x, knots[0], knots[knots.length -1]);
         }
         for (int i = 0; i < knots.length; i++) {
             if (knots[i] > x) {
                 return i - 1;
             }
         }
         throw new MathIllegalStateException(LocalizedCoreFormats.ILLEGAL_STATE);
     }
}

