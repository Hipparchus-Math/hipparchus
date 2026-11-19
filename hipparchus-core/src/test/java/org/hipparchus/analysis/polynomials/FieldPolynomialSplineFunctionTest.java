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

import org.hipparchus.analysis.differentiation.Gradient;
import org.hipparchus.analysis.differentiation.GradientField;
import org.hipparchus.analysis.interpolation.LinearInterpolator;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.Binary64;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.lang.reflect.Array;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Tests the FieldFieldPolynomialSplineFunction implementation.
 *
 */
class FieldPolynomialSplineFunctionTest {

    /** Error tolerance for tests */
    private double tolerance;

    /**
     * Quadratic polynomials used in tests:
     *
     * x^2 + x            [-1, 0)
     * x^2 + x + 2        [0, 1)
     * x^2 + x + 4        [1, 2)
     *
     * Defined so that evaluation using FieldPolynomialSplineFunction evaluation
     * algorithm agrees at knot point boundaries.
     */
    private FieldPolynomialFunction<Binary64>[] polynomials;

    /** Knot points  */
    private Binary64[] knots;

    /** Derivative of test polynomials -- 2x + 1  */
    protected PolynomialFunction dp =
        new PolynomialFunction(new double[] {1d, 2d});


    @Test
    void testConstructor() {
        FieldPolynomialSplineFunction<Binary64> spline =
            new FieldPolynomialSplineFunction<>(knots, polynomials);
        assertArrayEquals(knots, spline.getKnots());
        assertEquals(1d, spline.getPolynomials()[0].getCoefficients()[2].getReal(), 0);
        assertEquals(3, spline.getN());

        try { // too few knots
            new FieldPolynomialSplineFunction<>(new Binary64[] { new Binary64(0) }, polynomials);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try { // too many knots
            new FieldPolynomialSplineFunction<>(new Binary64[] {
                new Binary64(0), new Binary64(1), new Binary64(2),
                new Binary64(3), new Binary64(4)
            }, polynomials);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try { // knots not increasing
            new FieldPolynomialSplineFunction<>(new Binary64[] {
                new Binary64(0), new Binary64(1), new Binary64(3), new Binary64(2)
            }, polynomials);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testValues() {
        FieldPolynomialSplineFunction<Binary64> spline =
            new FieldPolynomialSplineFunction<>(knots, polynomials);
        FieldPolynomialSplineFunction<Binary64> dSpline = spline.polynomialSplineDerivative();

        /**
         * interior points -- spline value at x should equal p(x - knot)
         * where knot is the largest knot point less than or equal to x and p
         * is the polynomial defined over the knot segment to which x belongs.
         */
        double x = -1;
        int index = 0;
        for (int i = 0; i < 10; i++) {
           x += 0.25;
           index = findKnot(knots, x);
           assertEquals(polynomials[index].value(new Binary64(x).subtract(knots[index])).getReal(),
                               spline.value(x).getReal(),
                               tolerance,
                               "spline function evaluation failed for x=" + x);
           assertEquals(dp.value(new Binary64(x).subtract(knots[index])).getReal(),
                               dSpline.value(x).getReal(),
                               tolerance,
                               "spline derivative evaluation failed for x=" + x);
        }

        // knot points -- centering should zero arguments
        for (int i = 0; i < 3; i++) {
            assertEquals(polynomials[i].value(0).getReal(),
                                spline.value(knots[i]).getReal(),
                                tolerance,
                                "spline function evaluation failed for knot=" + knots[i].getReal());
            assertEquals(dp.value(0),
                                dSpline.value(knots[i]).getReal(),
                                tolerance,
                                "spline function evaluation failed for knot=" + knots[i]);
        }

        try { //outside of domain -- under min
            spline.value(-1.5);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try { //outside of domain -- over max
            spline.value(2.5);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testIsValidPoint() {
        final FieldPolynomialSplineFunction<Binary64> spline =
            new FieldPolynomialSplineFunction<>(knots, polynomials);
        final Binary64 xMin = knots[0];
        final Binary64 xMax = knots[knots.length - 1];

        Binary64 x;

        x = xMin;
        assertTrue(spline.isValidPoint(x));
        // Ensure that no exception is thrown.
        spline.value(x);

        x = xMax;
        assertTrue(spline.isValidPoint(x));
        // Ensure that no exception is thrown.
        spline.value(x);

        final Binary64 xRange = xMax.subtract(xMin);
        x = xMin.add(xRange.divide(3.4));
        assertTrue(spline.isValidPoint(x));
        // Ensure that no exception is thrown.
        spline.value(x);

        final Binary64 small = new Binary64(1e-8);
        x = xMin.subtract(small);
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
     protected int findKnot(Binary64[] knots, double x) {
         if (x < knots[0].getReal() || x >= knots[knots.length -1].getReal()) {
             throw new MathIllegalArgumentException(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE,
                                                    x, knots[0], knots[knots.length -1]);
         }
         for (int i = 0; i < knots.length; i++) {
             if (knots[i].getReal() > x) {
                 return i - 1;
             }
         }
         throw new MathIllegalStateException(LocalizedCoreFormats.ILLEGAL_STATE);
     }

     private FieldPolynomialFunction<Binary64> buildD64(double...c) {
         Binary64[] array = new Binary64[c.length];
         for (int i = 0; i < c.length; ++i) {
             array[i] = new Binary64(c[i]);
         }
         return new FieldPolynomialFunction<>(array);
     }

    @BeforeEach
    @SuppressWarnings("unchecked")
    void setUp() {
         tolerance = 1.0e-12;
         polynomials = (FieldPolynomialFunction<Binary64>[]) Array.newInstance(FieldPolynomialFunction.class, 3);
         polynomials[0] = buildD64(0, 1, 1);
         polynomials[1] = buildD64(2, 1, 1);
         polynomials[2] = buildD64(4, 1, 1);
         knots = new Binary64[] {
             new Binary64(-1), new Binary64(0), new Binary64(1), new Binary64(2)
         };
     }

    @Test
    void testValueGradient() {
         // issue #257
         // Given
         final int freeParameters = 1;
         final GradientField field = GradientField.getField(freeParameters);
         final Gradient zero = field.getZero();
         final Gradient time1 = zero.add(1.0);
         final Gradient time2 = zero.add(2.0);
         final Gradient time3 = zero.add(4.0);
         final Gradient x1 = zero.add(4.0);
         final Gradient x2 = Gradient.variable(freeParameters, 0, -2.0);
         final Gradient x3 = zero.add(0.0);
         final Gradient[] times = new Gradient[] {time1, time2, time3};
         final Gradient[] xs = new Gradient[] {x1, x2, x3};
         // When
         final FieldPolynomialSplineFunction<Gradient> spline = new LinearInterpolator().interpolate(times, xs);
         final Gradient actualEvaluation = spline.value(zero.add(3.));
         // Then
         final Gradient expectedEvaluation = Gradient.variable(freeParameters, 0, 0).multiply(0.5).add(-1.0);
         assertEquals(expectedEvaluation, actualEvaluation);
     }

}

