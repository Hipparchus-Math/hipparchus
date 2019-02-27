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

import java.lang.reflect.Array;
import java.util.Arrays;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.Decimal64;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Tests the FieldFieldPolynomialSplineFunction implementation.
 *
 */
public class FieldPolynomialSplineFunctionTest {

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
    private FieldPolynomialFunction<Decimal64>[] polynomials;

    /** Knot points  */
    private Decimal64[] knots;

    /** Derivative of test polynomials -- 2x + 1  */
    protected PolynomialFunction dp =
        new PolynomialFunction(new double[] {1d, 2d});


    @Test
    public void testConstructor() {
        FieldPolynomialSplineFunction<Decimal64> spline =
            new FieldPolynomialSplineFunction<>(knots, polynomials);
        Assert.assertTrue(Arrays.equals(knots, spline.getKnots()));
        Assert.assertEquals(1d, spline.getPolynomials()[0].getCoefficients()[2].getReal(), 0);
        Assert.assertEquals(3, spline.getN());

        try { // too few knots
            new FieldPolynomialSplineFunction<>(new Decimal64[] { new Decimal64(0) }, polynomials);
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try { // too many knots
            new FieldPolynomialSplineFunction<>(new Decimal64[] {
                new Decimal64(0), new Decimal64(1), new Decimal64(2),
                new Decimal64(3), new Decimal64(4)
            }, polynomials);
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try { // knots not increasing
            new FieldPolynomialSplineFunction<>(new Decimal64[] {
                new Decimal64(0), new Decimal64(1), new Decimal64(3), new Decimal64(2)
            }, polynomials);
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testValues() {
        FieldPolynomialSplineFunction<Decimal64> spline =
            new FieldPolynomialSplineFunction<>(knots, polynomials);
        FieldPolynomialSplineFunction<Decimal64> dSpline = spline.polynomialSplineDerivative();

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
           Assert.assertEquals("spline function evaluation failed for x=" + x,
                               polynomials[index].value(new Decimal64(x).subtract(knots[index])).getReal(),
                               spline.value(x).getReal(),
                               tolerance);
           Assert.assertEquals("spline derivative evaluation failed for x=" + x,
                               dp.value(new Decimal64(x).subtract(knots[index])).getReal(),
                               dSpline.value(x).getReal(),
                               tolerance);
        }

        // knot points -- centering should zero arguments
        for (int i = 0; i < 3; i++) {
            Assert.assertEquals("spline function evaluation failed for knot=" + knots[i].getReal(),
                                polynomials[i].value(0).getReal(),
                                spline.value(knots[i]).getReal(),
                                tolerance);
            Assert.assertEquals("spline function evaluation failed for knot=" + knots[i],
                                dp.value(0),
                                dSpline.value(knots[i]).getReal(),
                                tolerance);
        }

        try { //outside of domain -- under min
            spline.value(-1.5);
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try { //outside of domain -- over max
            spline.value(2.5);
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testIsValidPoint() {
        final FieldPolynomialSplineFunction<Decimal64> spline =
            new FieldPolynomialSplineFunction<>(knots, polynomials);
        final Decimal64 xMin = knots[0];
        final Decimal64 xMax = knots[knots.length - 1];

        Decimal64 x;

        x = xMin;
        Assert.assertTrue(spline.isValidPoint(x));
        // Ensure that no exception is thrown.
        spline.value(x);

        x = xMax;
        Assert.assertTrue(spline.isValidPoint(x));
        // Ensure that no exception is thrown.
        spline.value(x);

        final Decimal64 xRange = xMax.subtract(xMin);
        x = xMin.add(xRange.divide(3.4));
        Assert.assertTrue(spline.isValidPoint(x));
        // Ensure that no exception is thrown.
        spline.value(x);

        final Decimal64 small = new Decimal64(1e-8);
        x = xMin.subtract(small);
        Assert.assertFalse(spline.isValidPoint(x));
        // Ensure that an exception would have been thrown.
        try {
            spline.value(x);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException expected) {}
    }

    /**
     *  Do linear search to find largest knot point less than or equal to x.
     *  Implementation does binary search.
     */
     protected int findKnot(Decimal64[] knots, double x) {
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

     private FieldPolynomialFunction<Decimal64> buildD64(double...c) {
         Decimal64[] array = new Decimal64[c.length];
         for (int i = 0; i < c.length; ++i) {
             array[i] = new Decimal64(c[i]);
         }
         return new FieldPolynomialFunction<>(array);
     }

     @Before
     @SuppressWarnings("unchecked")
     public void setUp() {
         tolerance = 1.0e-12;
         polynomials = (FieldPolynomialFunction<Decimal64>[]) Array.newInstance(FieldPolynomialFunction.class, 3);
         polynomials[0] = buildD64(0, 1, 1);
         polynomials[1] = buildD64(2, 1, 1);
         polynomials[2] = buildD64(4, 1, 1);
         knots = new Decimal64[] {
             new Decimal64(-1), new Decimal64(0), new Decimal64(1), new Decimal64(2)
         };
     }

}

