/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements. See the NOTICE file distributed with this
 * work for additional information regarding copyright ownership. The ASF
 * licenses this file to You under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * https://www.apache.org/licenses/LICENSE-2.0 Unless required by applicable law
 * or agreed to in writing, software distributed under the License is
 * distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied. See the License for the specific language
 * governing permissions and limitations under the License.
 */
package org.hipparchus.util;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathRuntimeException;
import org.junit.jupiter.api.Test;

import java.math.RoundingMode;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the {@link Precision} class.
 */
class PrecisionTest {
    @Test
    void testEqualsWithRelativeTolerance() {
        assertTrue(Precision.equalsWithRelativeTolerance(0d, 0d, 0d));
        assertTrue(Precision.equalsWithRelativeTolerance(0d, 1 / Double.NEGATIVE_INFINITY, 0d));

        final double eps = 1e-14;
        assertFalse(Precision.equalsWithRelativeTolerance(1.987654687654968, 1.987654687654988, eps));
        assertTrue(Precision.equalsWithRelativeTolerance(1.987654687654968, 1.987654687654987, eps));
        assertFalse(Precision.equalsWithRelativeTolerance(1.987654687654968, 1.987654687654948, eps));
        assertTrue(Precision.equalsWithRelativeTolerance(1.987654687654968, 1.987654687654949, eps));

        assertFalse(Precision.equalsWithRelativeTolerance(Precision.SAFE_MIN, 0.0, eps));

        assertFalse(Precision.equalsWithRelativeTolerance(1.0000000000001e-300, 1e-300, eps));
        assertTrue(Precision.equalsWithRelativeTolerance(1.00000000000001e-300, 1e-300, eps));

        assertFalse(Precision.equalsWithRelativeTolerance(Double.NEGATIVE_INFINITY, 1.23, eps));
        assertFalse(Precision.equalsWithRelativeTolerance(Double.POSITIVE_INFINITY, 1.23, eps));

        assertTrue(Precision.equalsWithRelativeTolerance(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, eps));
        assertTrue(Precision.equalsWithRelativeTolerance(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, eps));
        assertFalse(Precision.equalsWithRelativeTolerance(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, eps));

        assertFalse(Precision.equalsWithRelativeTolerance(Double.NaN, 1.23, eps));
        assertFalse(Precision.equalsWithRelativeTolerance(Double.NaN, Double.NaN, eps));
    }

    @Test
    void testEqualsIncludingNaN() {
        double[] testArray = {
            Double.NaN,
            Double.POSITIVE_INFINITY,
            Double.NEGATIVE_INFINITY,
            1d,
            0d };
        for (int i = 0; i < testArray.length; i++) {
            for (int j = 0; j < testArray.length; j++) {
                if (i == j) {
                    assertTrue(Precision.equalsIncludingNaN(testArray[i], testArray[j]));
                    assertTrue(Precision.equalsIncludingNaN(testArray[j], testArray[i]));
                } else {
                    assertFalse(Precision.equalsIncludingNaN(testArray[i], testArray[j]));
                    assertFalse(Precision.equalsIncludingNaN(testArray[j], testArray[i]));
                }
            }
        }
    }

    @Test
    void testEqualsWithAllowedDelta() {
        assertTrue(Precision.equals(153.0000, 153.0000, .0625));
        assertTrue(Precision.equals(153.0000, 153.0625, .0625));
        assertTrue(Precision.equals(152.9375, 153.0000, .0625));
        assertFalse(Precision.equals(153.0000, 153.0625, .0624));
        assertFalse(Precision.equals(152.9374, 153.0000, .0625));
        assertFalse(Precision.equals(Double.NaN, Double.NaN, 1.0));
        assertTrue(Precision.equals(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0));
        assertTrue(Precision.equals(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, 1.0));
        assertFalse(Precision.equals(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0));
    }

    @Test
    void testMath475() {
        final double a = 1.7976931348623182E16;
        final double b = FastMath.nextUp(a);

        double diff = FastMath.abs(a - b);
        // Because they are adjacent floating point numbers, "a" and "b" are
        // considered equal even though the allowed error is smaller than
        // their difference.
        assertTrue(Precision.equals(a, b, 0.5 * diff));

        final double c = FastMath.nextUp(b);
        diff = FastMath.abs(a - c);
        // Because "a" and "c" are not adjacent, the tolerance is taken into
        // account for assessing equality.
        assertTrue(Precision.equals(a, c, diff));
        assertFalse(Precision.equals(a, c, (1 - 1e-16) * diff));
    }

    @Test
    void testEqualsIncludingNaNWithAllowedDelta() {
        assertTrue(Precision.equalsIncludingNaN(153.0000, 153.0000, .0625));
        assertTrue(Precision.equalsIncludingNaN(153.0000, 153.0625, .0625));
        assertTrue(Precision.equalsIncludingNaN(152.9375, 153.0000, .0625));
        assertTrue(Precision.equalsIncludingNaN(Double.NaN, Double.NaN, 1.0));
        assertTrue(Precision.equalsIncludingNaN(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0));
        assertTrue(Precision.equalsIncludingNaN(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, 1.0));
        assertFalse(Precision.equalsIncludingNaN(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0));
        assertFalse(Precision.equalsIncludingNaN(153.0000, 153.0625, .0624));
        assertFalse(Precision.equalsIncludingNaN(152.9374, 153.0000, .0625));
    }

    // Tests for floating point equality
    @Test
    void testFloatEqualsWithAllowedUlps() {
        assertTrue(Precision.equals(0.0f, -0.0f),"+0.0f == -0.0f");
        assertTrue(Precision.equals(0.0f, -0.0f, 1),"+0.0f == -0.0f (1 ulp)");
        float oneFloat = 1.0f;
        assertTrue(Precision.equals(oneFloat, Float.intBitsToFloat(1 + Float.floatToIntBits(oneFloat))),"1.0f == 1.0f + 1 ulp");
        assertTrue(Precision.equals(oneFloat, Float.intBitsToFloat(1 + Float.floatToIntBits(oneFloat)), 1),"1.0f == 1.0f + 1 ulp (1 ulp)");
        assertFalse(Precision.equals(oneFloat, Float.intBitsToFloat(2 + Float.floatToIntBits(oneFloat)), 1),"1.0f != 1.0f + 2 ulp (1 ulp)");

        assertTrue(Precision.equals(153.0f, 153.0f, 1));

        // These tests need adjusting for floating point precision
//        Assertions.assertTrue(Precision.equals(153.0f, 153.00000000000003f, 1));
//        Assertions.assertFalse(Precision.equals(153.0f, 153.00000000000006f, 1));
//        Assertions.assertTrue(Precision.equals(153.0f, 152.99999999999997f, 1));
//        Assertions.assertFalse(Precision.equals(153f, 152.99999999999994f, 1));
//
//        Assertions.assertTrue(Precision.equals(-128.0f, -127.99999999999999f, 1));
//        Assertions.assertFalse(Precision.equals(-128.0f, -127.99999999999997f, 1));
//        Assertions.assertTrue(Precision.equals(-128.0f, -128.00000000000003f, 1));
//        Assertions.assertFalse(Precision.equals(-128.0f, -128.00000000000006f, 1));

        assertTrue(Precision.equals(Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY, 1));
        assertTrue(Precision.equals(Double.MAX_VALUE, Float.POSITIVE_INFINITY, 1));

        assertTrue(Precision.equals(Float.NEGATIVE_INFINITY, Float.NEGATIVE_INFINITY, 1));
        assertTrue(Precision.equals(-Float.MAX_VALUE, Float.NEGATIVE_INFINITY, 1));

        assertFalse(Precision.equals(Float.NaN, Float.NaN, 1));
        assertFalse(Precision.equals(Float.NaN, Float.NaN, 0));
        assertFalse(Precision.equals(Float.NaN, 0, 0));
        assertFalse(Precision.equals(Float.NaN, Float.POSITIVE_INFINITY, 0));
        assertFalse(Precision.equals(Float.NaN, Float.NEGATIVE_INFINITY, 0));

        assertFalse(Precision.equals(Float.NEGATIVE_INFINITY, Float.POSITIVE_INFINITY, 100000));
    }

    @Test
    void testEqualsWithAllowedUlps() {
        assertTrue(Precision.equals(0.0, -0.0, 1));

        assertTrue(Precision.equals(1.0, 1 + FastMath.ulp(1d), 1));
        assertFalse(Precision.equals(1.0, 1 + 2 * FastMath.ulp(1d), 1));

        final double nUp1 = FastMath.nextAfter(1d, Double.POSITIVE_INFINITY);
        final double nnUp1 = FastMath.nextAfter(nUp1, Double.POSITIVE_INFINITY);
        assertTrue(Precision.equals(1.0, nUp1, 1));
        assertTrue(Precision.equals(nUp1, nnUp1, 1));
        assertFalse(Precision.equals(1.0, nnUp1, 1));

        assertTrue(Precision.equals(0.0, FastMath.ulp(0d), 1));
        assertTrue(Precision.equals(0.0, -FastMath.ulp(0d), 1));

        assertTrue(Precision.equals(153.0, 153.0, 1));

        assertTrue(Precision.equals(153.0, 153.00000000000003, 1));
        assertFalse(Precision.equals(153.0, 153.00000000000006, 1));
        assertTrue(Precision.equals(153.0, 152.99999999999997, 1));
        assertFalse(Precision.equals(153, 152.99999999999994, 1));

        assertTrue(Precision.equals(-128.0, -127.99999999999999, 1));
        assertFalse(Precision.equals(-128.0, -127.99999999999997, 1));
        assertTrue(Precision.equals(-128.0, -128.00000000000003, 1));
        assertFalse(Precision.equals(-128.0, -128.00000000000006, 1));

        assertTrue(Precision.equals(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, 1));
        assertTrue(Precision.equals(Double.MAX_VALUE, Double.POSITIVE_INFINITY, 1));

        assertTrue(Precision.equals(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        assertTrue(Precision.equals(-Double.MAX_VALUE, Double.NEGATIVE_INFINITY, 1));

        assertFalse(Precision.equals(Double.NaN, Double.NaN, 1));
        assertFalse(Precision.equals(Double.NaN, Double.NaN, 0));
        assertFalse(Precision.equals(Double.NaN, 0, 0));
        assertFalse(Precision.equals(Double.NaN, Double.POSITIVE_INFINITY, 0));
        assertFalse(Precision.equals(Double.NaN, Double.NEGATIVE_INFINITY, 0));

        assertFalse(Precision.equals(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 100000));
    }

    @Test
    void testEqualsIncludingNaNWithAllowedUlps() {
        assertTrue(Precision.equalsIncludingNaN(0.0, -0.0, 1));

        assertTrue(Precision.equalsIncludingNaN(1.0, 1 + FastMath.ulp(1d), 1));
        assertFalse(Precision.equalsIncludingNaN(1.0, 1 + 2 * FastMath.ulp(1d), 1));

        final double nUp1 = FastMath.nextAfter(1d, Double.POSITIVE_INFINITY);
        final double nnUp1 = FastMath.nextAfter(nUp1, Double.POSITIVE_INFINITY);
        assertTrue(Precision.equalsIncludingNaN(1.0, nUp1, 1));
        assertTrue(Precision.equalsIncludingNaN(nUp1, nnUp1, 1));
        assertFalse(Precision.equalsIncludingNaN(1.0, nnUp1, 1));

        assertTrue(Precision.equalsIncludingNaN(0.0, FastMath.ulp(0d), 1));
        assertTrue(Precision.equalsIncludingNaN(0.0, -FastMath.ulp(0d), 1));

        assertTrue(Precision.equalsIncludingNaN(153.0, 153.0, 1));

        assertTrue(Precision.equalsIncludingNaN(153.0, 153.00000000000003, 1));
        assertFalse(Precision.equalsIncludingNaN(153.0, 153.00000000000006, 1));
        assertTrue(Precision.equalsIncludingNaN(153.0, 152.99999999999997, 1));
        assertFalse(Precision.equalsIncludingNaN(153, 152.99999999999994, 1));

        assertTrue(Precision.equalsIncludingNaN(-128.0, -127.99999999999999, 1));
        assertFalse(Precision.equalsIncludingNaN(-128.0, -127.99999999999997, 1));
        assertTrue(Precision.equalsIncludingNaN(-128.0, -128.00000000000003, 1));
        assertFalse(Precision.equalsIncludingNaN(-128.0, -128.00000000000006, 1));

        assertTrue(Precision.equalsIncludingNaN(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, 1));
        assertTrue(Precision.equalsIncludingNaN(Double.MAX_VALUE, Double.POSITIVE_INFINITY, 1));

        assertTrue(Precision.equalsIncludingNaN(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        assertTrue(Precision.equalsIncludingNaN(-Double.MAX_VALUE, Double.NEGATIVE_INFINITY, 1));

        assertTrue(Precision.equalsIncludingNaN(Double.NaN, Double.NaN, 1));

        assertFalse(Precision.equalsIncludingNaN(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 100000));
    }

    @Test
    void testCompareToEpsilon() {
        assertEquals(0, Precision.compareTo(152.33, 152.32, .011));
        assertTrue(Precision.compareTo(152.308, 152.32, .011) < 0);
        assertTrue(Precision.compareTo(152.33, 152.318, .011) > 0);
        assertEquals(0, Precision.compareTo(Double.MIN_VALUE, +0.0, Double.MIN_VALUE));
        assertEquals(0, Precision.compareTo(Double.MIN_VALUE, -0.0, Double.MIN_VALUE));
    }

    @Test
    void testCompareToMaxUlps() {
        double a     = 152.32;
        double delta = FastMath.ulp(a);
        for (int i = 0; i <= 10; ++i) {
            if (i <= 5) {
                assertEquals( 0, Precision.compareTo(a, a + i * delta, 5));
                assertEquals( 0, Precision.compareTo(a, a - i * delta, 5));
            } else {
                assertEquals(-1, Precision.compareTo(a, a + i * delta, 5));
                assertEquals(+1, Precision.compareTo(a, a - i * delta, 5));
            }
        }

        assertEquals( 0, Precision.compareTo(-0.0, 0.0, 0));

        assertEquals(-1, Precision.compareTo(-Double.MIN_VALUE, -0.0, 0));
        assertEquals( 0, Precision.compareTo(-Double.MIN_VALUE, -0.0, 1));
        assertEquals(-1, Precision.compareTo(-Double.MIN_VALUE, +0.0, 0));
        assertEquals( 0, Precision.compareTo(-Double.MIN_VALUE, +0.0, 1));

        assertEquals(+1, Precision.compareTo( Double.MIN_VALUE, -0.0, 0));
        assertEquals( 0, Precision.compareTo( Double.MIN_VALUE, -0.0, 1));
        assertEquals(+1, Precision.compareTo( Double.MIN_VALUE, +0.0, 0));
        assertEquals( 0, Precision.compareTo( Double.MIN_VALUE, +0.0, 1));

        assertEquals(-1, Precision.compareTo(-Double.MIN_VALUE, Double.MIN_VALUE, 0));
        assertEquals(-1, Precision.compareTo(-Double.MIN_VALUE, Double.MIN_VALUE, 1));
        assertEquals( 0, Precision.compareTo(-Double.MIN_VALUE, Double.MIN_VALUE, 2));

        assertEquals( 0, Precision.compareTo(Double.MAX_VALUE, Double.POSITIVE_INFINITY, 1));
        assertEquals(-1, Precision.compareTo(Double.MAX_VALUE, Double.POSITIVE_INFINITY, 0));

        assertEquals(+1, Precision.compareTo(Double.MAX_VALUE, Double.NaN, Integer.MAX_VALUE));
        assertEquals(+1, Precision.compareTo(Double.NaN, Double.MAX_VALUE, Integer.MAX_VALUE));
    }

    @Test
    void testRoundDouble() {
        double x = 1.234567890;
        assertEquals(1.23, Precision.round(x, 2), 0.0);
        assertEquals(1.235, Precision.round(x, 3), 0.0);
        assertEquals(1.2346, Precision.round(x, 4), 0.0);

        // JIRA MATH-151
        assertEquals(39.25, Precision.round(39.245, 2), 0.0);
        assertEquals(39.24, Precision.round(39.245, 2, RoundingMode.DOWN), 0.0);
        double xx = 39.0;
        xx += 245d / 1000d;
        assertEquals(39.25, Precision.round(xx, 2), 0.0);

        // BZ 35904
        assertEquals(30.1d, Precision.round(30.095d, 2), 0.0d);
        assertEquals(30.1d, Precision.round(30.095d, 1), 0.0d);
        assertEquals(33.1d, Precision.round(33.095d, 1), 0.0d);
        assertEquals(33.1d, Precision.round(33.095d, 2), 0.0d);
        assertEquals(50.09d, Precision.round(50.085d, 2), 0.0d);
        assertEquals(50.19d, Precision.round(50.185d, 2), 0.0d);
        assertEquals(50.01d, Precision.round(50.005d, 2), 0.0d);
        assertEquals(30.01d, Precision.round(30.005d, 2), 0.0d);
        assertEquals(30.65d, Precision.round(30.645d, 2), 0.0d);

        assertEquals(1.24, Precision.round(x, 2, RoundingMode.CEILING), 0.0);
        assertEquals(1.235, Precision.round(x, 3, RoundingMode.CEILING), 0.0);
        assertEquals(1.2346, Precision.round(x, 4, RoundingMode.CEILING), 0.0);
        assertEquals(-1.23, Precision.round(-x, 2, RoundingMode.CEILING), 0.0);
        assertEquals(-1.234, Precision.round(-x, 3, RoundingMode.CEILING), 0.0);
        assertEquals(-1.2345, Precision.round(-x, 4, RoundingMode.CEILING), 0.0);

        assertEquals(1.23, Precision.round(x, 2, RoundingMode.DOWN), 0.0);
        assertEquals(1.234, Precision.round(x, 3, RoundingMode.DOWN), 0.0);
        assertEquals(1.2345, Precision.round(x, 4, RoundingMode.DOWN), 0.0);
        assertEquals(-1.23, Precision.round(-x, 2, RoundingMode.DOWN), 0.0);
        assertEquals(-1.234, Precision.round(-x, 3, RoundingMode.DOWN), 0.0);
        assertEquals(-1.2345, Precision.round(-x, 4, RoundingMode.DOWN), 0.0);

        assertEquals(1.23, Precision.round(x, 2, RoundingMode.FLOOR), 0.0);
        assertEquals(1.234, Precision.round(x, 3, RoundingMode.FLOOR), 0.0);
        assertEquals(1.2345, Precision.round(x, 4, RoundingMode.FLOOR), 0.0);
        assertEquals(-1.24, Precision.round(-x, 2, RoundingMode.FLOOR), 0.0);
        assertEquals(-1.235, Precision.round(-x, 3, RoundingMode.FLOOR), 0.0);
        assertEquals(-1.2346, Precision.round(-x, 4, RoundingMode.FLOOR), 0.0);

        assertEquals(1.23, Precision.round(x, 2, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(1.235, Precision.round(x, 3, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(1.2346, Precision.round(x, 4, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(-1.23, Precision.round(-x, 2, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(-1.235, Precision.round(-x, 3, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(-1.2346, Precision.round(-x, 4, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(1.234, Precision.round(1.2345, 3, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(-1.234, Precision.round(-1.2345, 3, RoundingMode.HALF_DOWN), 0.0);

        assertEquals(1.23, Precision.round(x, 2, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(1.235, Precision.round(x, 3, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(1.2346, Precision.round(x, 4, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(-1.23, Precision.round(-x, 2, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(-1.235, Precision.round(-x, 3, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(-1.2346, Precision.round(-x, 4, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(1.234, Precision.round(1.2345, 3, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(-1.234, Precision.round(-1.2345, 3, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(1.236, Precision.round(1.2355, 3, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(-1.236, Precision.round(-1.2355, 3, RoundingMode.HALF_EVEN), 0.0);

        assertEquals(1.23, Precision.round(x, 2, RoundingMode.HALF_UP), 0.0);
        assertEquals(1.235, Precision.round(x, 3, RoundingMode.HALF_UP), 0.0);
        assertEquals(1.2346, Precision.round(x, 4, RoundingMode.HALF_UP), 0.0);
        assertEquals(-1.23, Precision.round(-x, 2, RoundingMode.HALF_UP), 0.0);
        assertEquals(-1.235, Precision.round(-x, 3, RoundingMode.HALF_UP), 0.0);
        assertEquals(-1.2346, Precision.round(-x, 4, RoundingMode.HALF_UP), 0.0);
        assertEquals(1.235, Precision.round(1.2345, 3, RoundingMode.HALF_UP), 0.0);
        assertEquals(-1.235, Precision.round(-1.2345, 3, RoundingMode.HALF_UP), 0.0);

        assertEquals(-1.23, Precision.round(-1.23, 2, RoundingMode.UNNECESSARY), 0.0);
        assertEquals(1.23, Precision.round(1.23, 2, RoundingMode.UNNECESSARY), 0.0);

        try {
            Precision.round(1.234, 2, RoundingMode.UNNECESSARY);
            fail();
        } catch (ArithmeticException ex) {
            // expected
        }

        assertEquals(1.24, Precision.round(x, 2, RoundingMode.UP), 0.0);
        assertEquals(1.235, Precision.round(x, 3, RoundingMode.UP), 0.0);
        assertEquals(1.2346, Precision.round(x, 4, RoundingMode.UP), 0.0);
        assertEquals(-1.24, Precision.round(-x, 2, RoundingMode.UP), 0.0);
        assertEquals(-1.235, Precision.round(-x, 3, RoundingMode.UP), 0.0);
        assertEquals(-1.2346, Precision.round(-x, 4, RoundingMode.UP), 0.0);

        // MATH-151
        assertEquals(39.25, Precision.round(39.245, 2, RoundingMode.HALF_UP), 0.0);

        // special values
        UnitTestUtils.customAssertEquals(Double.NaN, Precision.round(Double.NaN, 2), 0.0);
        assertEquals(0.0, Precision.round(0.0, 2), 0.0);
        assertEquals(Double.POSITIVE_INFINITY, Precision.round(Double.POSITIVE_INFINITY, 2), 0.0);
        assertEquals(Double.NEGATIVE_INFINITY, Precision.round(Double.NEGATIVE_INFINITY, 2), 0.0);
        // comparison of positive and negative zero is not possible -> always equal thus do string comparison
        assertEquals("-0.0", Double.toString(Precision.round(-0.0, 0)));
        assertEquals("-0.0", Double.toString(Precision.round(-1e-10, 0)));
    }

    @Test
    void testRoundFloat() {
        float x = 1.234567890f;
        assertEquals(1.23f, Precision.round(x, 2), 0.0);
        assertEquals(1.235f, Precision.round(x, 3), 0.0);
        assertEquals(1.2346f, Precision.round(x, 4), 0.0);

        // BZ 35904
        assertEquals(30.1f, Precision.round(30.095f, 2), 0.0f);
        assertEquals(30.1f, Precision.round(30.095f, 1), 0.0f);
        assertEquals(50.09f, Precision.round(50.085f, 2), 0.0f);
        assertEquals(50.19f, Precision.round(50.185f, 2), 0.0f);
        assertEquals(50.01f, Precision.round(50.005f, 2), 0.0f);
        assertEquals(30.01f, Precision.round(30.005f, 2), 0.0f);
        assertEquals(30.65f, Precision.round(30.645f, 2), 0.0f);

        assertEquals(1.24f, Precision.round(x, 2, RoundingMode.CEILING), 0.0);
        assertEquals(1.235f, Precision.round(x, 3, RoundingMode.CEILING), 0.0);
        assertEquals(1.2346f, Precision.round(x, 4, RoundingMode.CEILING), 0.0);
        assertEquals(-1.23f, Precision.round(-x, 2, RoundingMode.CEILING), 0.0);
        assertEquals(-1.234f, Precision.round(-x, 3, RoundingMode.CEILING), 0.0);
        assertEquals(-1.2345f, Precision.round(-x, 4, RoundingMode.CEILING), 0.0);

        assertEquals(1.23f, Precision.round(x, 2, RoundingMode.DOWN), 0.0);
        assertEquals(1.234f, Precision.round(x, 3, RoundingMode.DOWN), 0.0);
        assertEquals(1.2345f, Precision.round(x, 4, RoundingMode.DOWN), 0.0);
        assertEquals(-1.23f, Precision.round(-x, 2, RoundingMode.DOWN), 0.0);
        assertEquals(-1.234f, Precision.round(-x, 3, RoundingMode.DOWN), 0.0);
        assertEquals(-1.2345f, Precision.round(-x, 4, RoundingMode.DOWN), 0.0);

        assertEquals(1.23f, Precision.round(x, 2, RoundingMode.FLOOR), 0.0);
        assertEquals(1.234f, Precision.round(x, 3, RoundingMode.FLOOR), 0.0);
        assertEquals(1.2345f, Precision.round(x, 4, RoundingMode.FLOOR), 0.0);
        assertEquals(-1.24f, Precision.round(-x, 2, RoundingMode.FLOOR), 0.0);
        assertEquals(-1.235f, Precision.round(-x, 3, RoundingMode.FLOOR), 0.0);
        assertEquals(-1.2346f, Precision.round(-x, 4, RoundingMode.FLOOR), 0.0);

        assertEquals(1.23f, Precision.round(x, 2, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(1.235f, Precision.round(x, 3, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(1.2346f, Precision.round(x, 4, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(-1.23f, Precision.round(-x, 2, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(-1.235f, Precision.round(-x, 3, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(-1.2346f, Precision.round(-x, 4, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(1.234f, Precision.round(1.2345f, 3, RoundingMode.HALF_DOWN), 0.0);
        assertEquals(-1.234f, Precision.round(-1.2345f, 3, RoundingMode.HALF_DOWN), 0.0);

        assertEquals(1.23f, Precision.round(x, 2, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(1.235f, Precision.round(x, 3, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(1.2346f, Precision.round(x, 4, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(-1.23f, Precision.round(-x, 2, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(-1.235f, Precision.round(-x, 3, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(-1.2346f, Precision.round(-x, 4, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(1.234f, Precision.round(1.2345f, 3, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(-1.234f, Precision.round(-1.2345f, 3, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(1.236f, Precision.round(1.2355f, 3, RoundingMode.HALF_EVEN), 0.0);
        assertEquals(-1.236f, Precision.round(-1.2355f, 3, RoundingMode.HALF_EVEN), 0.0);

        assertEquals(1.23f, Precision.round(x, 2, RoundingMode.HALF_UP), 0.0);
        assertEquals(1.235f, Precision.round(x, 3, RoundingMode.HALF_UP), 0.0);
        assertEquals(1.2346f, Precision.round(x, 4, RoundingMode.HALF_UP), 0.0);
        assertEquals(-1.23f, Precision.round(-x, 2, RoundingMode.HALF_UP), 0.0);
        assertEquals(-1.235f, Precision.round(-x, 3, RoundingMode.HALF_UP), 0.0);
        assertEquals(-1.2346f, Precision.round(-x, 4, RoundingMode.HALF_UP), 0.0);
        assertEquals(1.235f, Precision.round(1.2345f, 3, RoundingMode.HALF_UP), 0.0);
        assertEquals(-1.235f, Precision.round(-1.2345f, 3, RoundingMode.HALF_UP), 0.0);

        assertEquals(-1.23f, Precision.round(-1.23f, 2, RoundingMode.UNNECESSARY), 0.0);
        assertEquals(1.23f, Precision.round(1.23f, 2, RoundingMode.UNNECESSARY), 0.0);

        try {
            Precision.round(1.234f, 2, RoundingMode.UNNECESSARY);
            fail();
        } catch (MathRuntimeException ex) {
            // success
        }

        assertEquals(1.24f, Precision.round(x, 2, RoundingMode.UP), 0.0);
        assertEquals(1.235f, Precision.round(x, 3, RoundingMode.UP), 0.0);
        assertEquals(1.2346f, Precision.round(x, 4, RoundingMode.UP), 0.0);
        assertEquals(-1.24f, Precision.round(-x, 2, RoundingMode.UP), 0.0);
        assertEquals(-1.235f, Precision.round(-x, 3, RoundingMode.UP), 0.0);
        assertEquals(-1.2346f, Precision.round(-x, 4, RoundingMode.UP), 0.0);

        // special values
        UnitTestUtils.customAssertEquals(Float.NaN, Precision.round(Float.NaN, 2), 0.0f);
        assertEquals(0.0f, Precision.round(0.0f, 2), 0.0f);
        assertEquals(Float.POSITIVE_INFINITY, Precision.round(Float.POSITIVE_INFINITY, 2), 0.0f);
        assertEquals(Float.NEGATIVE_INFINITY, Precision.round(Float.NEGATIVE_INFINITY, 2), 0.0f);
        // comparison of positive and negative zero is not possible -> always equal thus do string comparison
        assertEquals("-0.0", Float.toString(Precision.round(-0.0f, 0)));
        assertEquals("-0.0", Float.toString(Precision.round(-1e-10f, 0)));

        // MATH-1070
        assertEquals(0.0f, Precision.round(0f, 2, RoundingMode.UP), 0.0f);
        assertEquals(0.05f, Precision.round(0.05f, 2, RoundingMode.UP), 0.0f);
        assertEquals(0.06f, Precision.round(0.051f, 2, RoundingMode.UP), 0.0f);
        assertEquals(0.06f, Precision.round(0.0505f, 2, RoundingMode.UP), 0.0f);
        assertEquals(0.06f, Precision.round(0.059f, 2, RoundingMode.UP), 0.0f);
    }


    @Test
    void testIssue721() {
        assertEquals(-53,   FastMath.getExponent(Precision.EPSILON));
        assertEquals(-1022, FastMath.getExponent(Precision.SAFE_MIN));
    }


    @Test
    void testRepresentableDelta() {
        int nonRepresentableCount = 0;
        final double x = 100;
        final int numTrials = 10000;
        for (int i = 0; i < numTrials; i++) {
            final double originalDelta = FastMath.random();
            final double delta = Precision.representableDelta(x, originalDelta);
            if (delta != originalDelta) {
                ++nonRepresentableCount;
            }
        }

        assertTrue(nonRepresentableCount / (double) numTrials > 0.9);
    }

    @Test
    void testMath843() {
        final double afterEpsilon = FastMath.nextAfter(Precision.EPSILON,
                                                       Double.POSITIVE_INFINITY);

        // a) 1 + EPSILON is equal to 1.
        assertEquals(1, 1 + Precision.EPSILON);

        // b) 1 + "the number after EPSILON" is not equal to 1.
        assertNotEquals(1, 1 + afterEpsilon, 0.0);
    }

    @Test
    void testMath1127() {
        assertFalse(Precision.equals(2.0, -2.0, 1));
        assertTrue(Precision.equals(0.0, -0.0, 0));
        assertFalse(Precision.equals(2.0f, -2.0f, 1));
        assertTrue(Precision.equals(0.0f, -0.0f, 0));
    }

    @Test
    void testIsMathematicalIntegerDouble() {

        assertFalse(Precision.isMathematicalInteger(Double.NaN));
        assertFalse(Precision.isMathematicalInteger(Double.POSITIVE_INFINITY));
        assertFalse(Precision.isMathematicalInteger(Double.NEGATIVE_INFINITY));
        assertFalse(Precision.isMathematicalInteger(Double.MIN_NORMAL));
        assertFalse(Precision.isMathematicalInteger(Double.MIN_VALUE));

        assertTrue(Precision.isMathematicalInteger(-0.0));
        assertTrue(Precision.isMathematicalInteger(+0.0));

        for (int i = -1000; i < 1000; ++i) {
            final double d = i;
            assertTrue(Precision.isMathematicalInteger(d));
            assertFalse(Precision.isMathematicalInteger(FastMath.nextAfter(d, Double.POSITIVE_INFINITY)));
            assertFalse(Precision.isMathematicalInteger(FastMath.nextAfter(d, Double.NEGATIVE_INFINITY)));
        }

        double minNoFractional = 0x1l << 52;
        assertTrue(Precision.isMathematicalInteger(minNoFractional));
        assertFalse(Precision.isMathematicalInteger(minNoFractional - 0.5));
        assertTrue(Precision.isMathematicalInteger(minNoFractional + 0.5));

    }

    @Test
    void testIsMathematicalIntegerFloat() {

        assertFalse(Precision.isMathematicalInteger(Float.NaN));
        assertFalse(Precision.isMathematicalInteger(Float.POSITIVE_INFINITY));
        assertFalse(Precision.isMathematicalInteger(Float.NEGATIVE_INFINITY));
        assertFalse(Precision.isMathematicalInteger(Float.MIN_NORMAL));
        assertFalse(Precision.isMathematicalInteger(Float.MIN_VALUE));

        assertTrue(Precision.isMathematicalInteger(-0.0f));
        assertTrue(Precision.isMathematicalInteger(+0.0f));

        for (int i = -1000; i < 1000; ++i) {
            final float f = i;
            assertTrue(Precision.isMathematicalInteger(f));
            assertFalse(Precision.isMathematicalInteger(FastMath.nextAfter(f, Float.POSITIVE_INFINITY)));
            assertFalse(Precision.isMathematicalInteger(FastMath.nextAfter(f, Float.NEGATIVE_INFINITY)));
        }

        float minNoFractional = 0x1l << 23;
        assertTrue(Precision.isMathematicalInteger(minNoFractional));
        assertFalse(Precision.isMathematicalInteger(minNoFractional - 0.5f));
        assertTrue(Precision.isMathematicalInteger(minNoFractional + 0.5f));

    }

}
