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
import org.hipparchus.complex.Complex;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.Well1024a;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the {@link MathArrays} class.
 */
public class MathArraysTest {

    private double[] testArray = {0, 1, 2, 3, 4, 5};
    private double[] testWeightsArray = {0.3, 0.2, 1.3, 1.1, 1.0, 1.8};
    private double[] testNegativeWeightsArray = {-0.3, 0.2, -1.3, 1.1, 1.0, 1.8};
    private double[] nullArray = null;
    private double[] singletonArray = {0};

    @Test
    void testScale() {
        final double[] test = new double[] { -2.5, -1, 0, 1, 2.5 };
        final double[] correctTest = test.clone();
        final double[] correctScaled = new double[]{5.25, 2.1, 0, -2.1, -5.25};

        final double[] scaled = MathArrays.scale(-2.1, test);

        // Make sure test has not changed
        for (int i = 0; i < test.length; i++) {
            assertEquals(correctTest[i], test[i], 0);
        }

        // Test scaled values
        for (int i = 0; i < scaled.length; i++) {
            assertEquals(correctScaled[i], scaled[i], 0);
        }
    }

    @Test
    void testScaleInPlace() {
        final double[] test = new double[] { -2.5, -1, 0, 1, 2.5 };
        final double[] correctScaled = new double[]{5.25, 2.1, 0, -2.1, -5.25};
        MathArrays.scaleInPlace(-2.1, test);

        // Make sure test has changed
        for (int i = 0; i < test.length; i++) {
            assertEquals(correctScaled[i], test[i], 0);
        }
    }

    @Test
    void testEbeAddPrecondition() {
        assertThrows(MathIllegalArgumentException.class, () -> MathArrays.ebeAdd(new double[3], new double[4]));
    }

    @Test
    void testEbeSubtractPrecondition() {
        assertThrows(MathIllegalArgumentException.class, () -> MathArrays.ebeSubtract(new double[3], new double[4]));
    }

    @Test
    void testEbeMultiplyPrecondition() {
        assertThrows(MathIllegalArgumentException.class, () -> MathArrays.ebeMultiply(new double[3], new double[4]));
    }

    @Test
    void testEbeDividePrecondition() {
        assertThrows(MathIllegalArgumentException.class, () -> MathArrays.ebeDivide(new double[3], new double[4]));
    }

    @Test
    void testEbeAdd() {
        final double[] a = { 0, 1, 2 };
        final double[] b = { 3, 5, 7 };
        final double[] r = MathArrays.ebeAdd(a, b);

        for (int i = 0; i < a.length; i++) {
            assertEquals(a[i] + b[i], r[i], 0);
        }
    }

    @Test
    void testEbeSubtract() {
        final double[] a = { 0, 1, 2 };
        final double[] b = { 3, 5, 7 };
        final double[] r = MathArrays.ebeSubtract(a, b);

        for (int i = 0; i < a.length; i++) {
            assertEquals(a[i] - b[i], r[i], 0);
        }
    }

    @Test
    void testEbeMultiply() {
        final double[] a = { 0, 1, 2 };
        final double[] b = { 3, 5, 7 };
        final double[] r = MathArrays.ebeMultiply(a, b);

        for (int i = 0; i < a.length; i++) {
            assertEquals(a[i] * b[i], r[i], 0);
        }
    }

    @Test
    void testEbeDivide() {
        final double[] a = { 0, 1, 2 };
        final double[] b = { 3, 5, 7 };
        final double[] r = MathArrays.ebeDivide(a, b);

        for (int i = 0; i < a.length; i++) {
            assertEquals(a[i] / b[i], r[i], 0);
        }
    }

    @Test
    void testL1DistanceDouble() {
        double[] p1 = { 2.5,  0.0 };
        double[] p2 = { -0.5, 4.0 };
        assertTrue(Precision.equals(7.0, MathArrays.distance1(p1, p2), 1));
    }

    @Test
    void testL1DistanceInt() {
        int[] p1 = { 3, 0 };
        int[] p2 = { 0, 4 };
        assertEquals(7, MathArrays.distance1(p1, p2));
    }

    @Test
    void testL2DistanceDouble() {
        double[] p1 = { 2.5,  0.0 };
        double[] p2 = { -0.5, 4.0 };
        assertTrue(Precision.equals(5.0, MathArrays.distance(p1, p2), 1));
    }

    @Test
    void testL2DistanceInt() {
        int[] p1 = { 3, 0 };testArrayEquals();
        int[] p2 = { 0, 4 };
        assertTrue(Precision.equals(5, MathArrays.distance(p1, p2), 1));
    }

    @Test
    void testLInfDistanceDouble() {
        double[] p1 = { 2.5,  0.0 };
        double[] p2 = { -0.5, 4.0 };
        assertTrue(Precision.equals(4.0, MathArrays.distanceInf(p1, p2), 1));
    }

    @Test
    void testLInfDistanceInt() {
        int[] p1 = { 3, 0 };
        int[] p2 = { 0, 4 };
        assertEquals(4, MathArrays.distanceInf(p1, p2));
    }

    @Test
    void testCosAngle2D() {
        double expected;

        final double[] v1 = { 1, 0 };
        expected = 1;
        assertEquals(expected, MathArrays.cosAngle(v1, v1), 0d);

        final double[] v2 = { 0, 1 };
        expected = 0;
        assertEquals(expected, MathArrays.cosAngle(v1, v2), 0d);

        final double[] v3 = { 7, 7 };
        expected = Math.sqrt(2) / 2;
        assertEquals(expected, MathArrays.cosAngle(v1, v3), 1e-15);
        assertEquals(expected, MathArrays.cosAngle(v3, v2), 1e-15);

        final double[] v4 = { -5, 0 };
        expected = -1;
        assertEquals(expected, MathArrays.cosAngle(v1, v4), 0);

        final double[] v5 = { -100, 100 };
        expected = 0;
        assertEquals(expected, MathArrays.cosAngle(v3, v5), 0);
    }

    @Test
    void testCosAngle3D() {
        double expected;

        final double[] v1 = { 1, 1, 0 };
        expected = 1;
        assertEquals(expected, MathArrays.cosAngle(v1, v1), 1e-15);

        final double[] v2 = { 1, 1, 1 };
        expected = Math.sqrt(2) / Math.sqrt(3);
        assertEquals(expected, MathArrays.cosAngle(v1, v2), 1e-15);
    }

    @Test
    void testCosAngleExtreme() {
        double expected;

        final double tiny = 1e-200;
        final double[] v1 = { tiny, tiny };
        final double big = 1e200;
        final double[] v2 = { -big, -big };
        expected = -1;
        assertEquals(expected, MathArrays.cosAngle(v1, v2), 1e-15);

        final double[] v3 = { big, -big };
        expected = 0;
        assertEquals(expected, MathArrays.cosAngle(v1, v3), 1e-15);
    }

    @Test
    void testCheckOrder() {
        MathArrays.checkOrder(new double[] {-15, -5.5, -1, 2, 15},
                             MathArrays.OrderDirection.INCREASING, true);
        MathArrays.checkOrder(new double[] {-15, -5.5, -1, 2, 2},
                             MathArrays.OrderDirection.INCREASING, false);
        MathArrays.checkOrder(new double[] {3, -5.5, -11, -27.5},
                             MathArrays.OrderDirection.DECREASING, true);
        MathArrays.checkOrder(new double[] {3, 0, 0, -5.5, -11, -27.5},
                             MathArrays.OrderDirection.DECREASING, false);

        try {
            MathArrays.checkOrder(new double[] {-15, -5.5, -1, -1, 2, 15},
                                 MathArrays.OrderDirection.INCREASING, true);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        try {
            MathArrays.checkOrder(new double[] {-15, -5.5, -1, -2, 2},
                                 MathArrays.OrderDirection.INCREASING, false);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        try {
            MathArrays.checkOrder(new double[] {3, 3, -5.5, -11, -27.5},
                                 MathArrays.OrderDirection.DECREASING, true);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        try {
            MathArrays.checkOrder(new double[] {3, -1, 0, -5.5, -11, -27.5},
                                 MathArrays.OrderDirection.DECREASING, false);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        try {
            MathArrays.checkOrder(new double[] {3, 0, -5.5, -11, -10},
                                 MathArrays.OrderDirection.DECREASING, false);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
    }

    @Test
    void testIsMonotonic() {
        assertFalse(MathArrays.isMonotonic(new double[] { -15, -5.5, -1, -1, 2, 15 },
                                                  MathArrays.OrderDirection.INCREASING, true));
        assertTrue(MathArrays.isMonotonic(new double[] { -15, -5.5, -1, 0, 2, 15 },
                                                 MathArrays.OrderDirection.INCREASING, true));
        assertFalse(MathArrays.isMonotonic(new double[] { -15, -5.5, -1, -2, 2 },
                                                  MathArrays.OrderDirection.INCREASING, false));
        assertTrue(MathArrays.isMonotonic(new double[] { -15, -5.5, -1, -1, 2 },
                                                 MathArrays.OrderDirection.INCREASING, false));
        assertFalse(MathArrays.isMonotonic(new double[] { 3, 3, -5.5, -11, -27.5 },
                                                  MathArrays.OrderDirection.DECREASING, true));
        assertTrue(MathArrays.isMonotonic(new double[] { 3, 2, -5.5, -11, -27.5 },
                                                 MathArrays.OrderDirection.DECREASING, true));
        assertFalse(MathArrays.isMonotonic(new double[] { 3, -1, 0, -5.5, -11, -27.5 },
                                                  MathArrays.OrderDirection.DECREASING, false));
        assertTrue(MathArrays.isMonotonic(new double[] { 3, 0, 0, -5.5, -11, -27.5 },
                                                 MathArrays.OrderDirection.DECREASING, false));
    }

    @Test
    void testIsMonotonicComparable() {
        assertFalse(MathArrays.isMonotonic(new Double[] { Double.valueOf(-15),
                                                                 Double.valueOf(-5.5),
                                                                 Double.valueOf(-1),
                                                                 Double.valueOf(-1),
                                                                 Double.valueOf(2),
                                                                 Double.valueOf(15) },
                MathArrays.OrderDirection.INCREASING, true));
        assertTrue(MathArrays.isMonotonic(new Double[] { Double.valueOf(-15),
                                                                Double.valueOf(-5.5),
                                                                Double.valueOf(-1),
                                                                Double.valueOf(0),
                                                                Double.valueOf(2),
                                                                Double.valueOf(15) },
                MathArrays.OrderDirection.INCREASING, true));
        assertFalse(MathArrays.isMonotonic(new Double[] { Double.valueOf(-15),
                                                                 Double.valueOf(-5.5),
                                                                 Double.valueOf(-1),
                                                                 Double.valueOf(-2),
                                                                 Double.valueOf(2) },
                MathArrays.OrderDirection.INCREASING, false));
        assertTrue(MathArrays.isMonotonic(new Double[] { Double.valueOf(-15),
                                                                Double.valueOf(-5.5),
                                                                Double.valueOf(-1),
                                                                Double.valueOf(-1),
                                                                Double.valueOf(2) },
                MathArrays.OrderDirection.INCREASING, false));
        assertFalse(MathArrays.isMonotonic(new Double[] { Double.valueOf(3),
                                                                 Double.valueOf(3),
                                                                 Double.valueOf(-5.5),
                                                                 Double.valueOf(-11),
                                                                 Double.valueOf(-27.5) },
                MathArrays.OrderDirection.DECREASING, true));
        assertTrue(MathArrays.isMonotonic(new Double[] { Double.valueOf(3),
                                                                Double.valueOf(2),
                                                                Double.valueOf(-5.5),
                                                                Double.valueOf(-11),
                                                                Double.valueOf(-27.5) },
                MathArrays.OrderDirection.DECREASING, true));
        assertFalse(MathArrays.isMonotonic(new Double[] { Double.valueOf(3),
                                                                 Double.valueOf(-1),
                                                                 Double.valueOf(0),
                                                                 Double.valueOf(-5.5),
                                                                 Double.valueOf(-11),
                                                                 Double.valueOf(-27.5) },
                MathArrays.OrderDirection.DECREASING, false));
        assertTrue(MathArrays.isMonotonic(new Double[] { Double.valueOf(3),
                                                                Double.valueOf(0),
                                                                Double.valueOf(0),
                                                                Double.valueOf(-5.5),
                                                                Double.valueOf(-11),
                                                                Double.valueOf(-27.5) },
                MathArrays.OrderDirection.DECREASING, false));
    }

    @Test
    void testCheckRectangular() {
        final long[][] rect = new long[][] {{0, 1}, {2, 3}};
        final long[][] ragged = new long[][] {{0, 1}, {2}};
        final long[][] nullArray = null;
        final long[][] empty = new long[][] {};
        MathArrays.checkRectangular(rect);
        MathArrays.checkRectangular(empty);
        try {
            MathArrays.checkRectangular(ragged);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
        try {
            MathArrays.checkRectangular(nullArray);
            fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // Expected
        }
    }

    @Test
    void testCheckPositive() {
        final double[] positive = new double[] {1, 2, 3};
        final double[] nonNegative = new double[] {0, 1, 2};
        final double[] nullArray = null;
        final double[] empty = new double[] {};
        MathArrays.checkPositive(positive);
        MathArrays.checkPositive(empty);
        try {
            MathArrays.checkPositive(nullArray);
            fail("Expecting NullPointerException");
        } catch (NullPointerException ex) {
            // Expected
        }
        try {
            MathArrays.checkPositive(nonNegative);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
    }

    @Test
    void testCheckNonNegative() {
        final long[] nonNegative = new long[] {0, 1};
        final long[] hasNegative = new long[] {-1};
        final long[] nullArray = null;
        final long[] empty = new long[] {};
        MathArrays.checkNonNegative(nonNegative);
        MathArrays.checkNonNegative(empty);
        try {
            MathArrays.checkNonNegative(nullArray);
            fail("Expecting NullPointerException");
        } catch (NullPointerException ex) {
            // Expected
        }
        try {
            MathArrays.checkNonNegative(hasNegative);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
    }

    @Test
    void testCheckNonNegative2D() {
        final long[][] nonNegative = new long[][] {{0, 1}, {1, 0}};
        final long[][] hasNegative = new long[][] {{-1}, {0}};
        final long[][] nullArray = null;
        final long[][] empty = new long[][] {};
        MathArrays.checkNonNegative(nonNegative);
        MathArrays.checkNonNegative(empty);
        try {
            MathArrays.checkNonNegative(nullArray);
            fail("Expecting NullPointerException");
        } catch (NullPointerException ex) {
            // Expected
        }
        try {
            MathArrays.checkNonNegative(hasNegative);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
    }

    @Test
    void testCheckNotNaN() {
        final double[] withoutNaN = { Double.NEGATIVE_INFINITY,
                                      -Double.MAX_VALUE,
                                      -1, 0,
                                      Double.MIN_VALUE,
                                      FastMath.ulp(1d),
                                      1, 3, 113, 4769,
                                      Double.MAX_VALUE,
                                      Double.POSITIVE_INFINITY };

        final double[] withNaN = { Double.NEGATIVE_INFINITY,
                                   -Double.MAX_VALUE,
                                   -1, 0,
                                   Double.MIN_VALUE,
                                   FastMath.ulp(1d),
                                   1, 3, 113, 4769,
                                   Double.MAX_VALUE,
                                   Double.POSITIVE_INFINITY,
                                   Double.NaN };


        final double[] nullArray = null;
        final double[] empty = new double[] {};
        MathArrays.checkNotNaN(withoutNaN);
        MathArrays.checkNotNaN(empty);
        try {
            MathArrays.checkNotNaN(nullArray);
            fail("Expecting NullPointerException");
        } catch (NullPointerException ex) {
            // Expected
        }
        try {
            MathArrays.checkNotNaN(withNaN);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // Expected
        }
    }

    @Test
    void testCheckEqualLength1() {
        assertThrows(MathIllegalArgumentException.class, () -> MathArrays.checkEqualLength(new double[]{1, 2, 3},
            new double[]{1, 2, 3, 4}));
    }

    @Test
    void testCheckEqualLength2() {
        final double[] a = new double[] {-1, -12, -23, -34};
        final double[] b = new double[] {56, 67, 78, 89};
        assertTrue(MathArrays.checkEqualLength(a, b, false));
    }

    @Test
    void testSortInPlace() {
        final double[] x1 = {2,   5,  -3, 1,  4};
        final double[] x2 = {4,  25,   9, 1, 16};
        final double[] x3 = {8, 125, -27, 1, 64};

        MathArrays.sortInPlace(x1, x2, x3);

        assertEquals(-3,  x1[0], FastMath.ulp(1d));
        assertEquals(9,   x2[0], FastMath.ulp(1d));
        assertEquals(-27, x3[0], FastMath.ulp(1d));

        assertEquals(1, x1[1], FastMath.ulp(1d));
        assertEquals(1, x2[1], FastMath.ulp(1d));
        assertEquals(1, x3[1], FastMath.ulp(1d));

        assertEquals(2, x1[2], FastMath.ulp(1d));
        assertEquals(4, x2[2], FastMath.ulp(1d));
        assertEquals(8, x3[2], FastMath.ulp(1d));

        assertEquals(4,  x1[3], FastMath.ulp(1d));
        assertEquals(16, x2[3], FastMath.ulp(1d));
        assertEquals(64, x3[3], FastMath.ulp(1d));

        assertEquals(5,   x1[4], FastMath.ulp(1d));
        assertEquals(25,  x2[4], FastMath.ulp(1d));
        assertEquals(125, x3[4], FastMath.ulp(1d));
    }

    @Test
    void testSortInPlaceDecresasingOrder() {
        final double[] x1 = {2,   5,  -3, 1,  4};
        final double[] x2 = {4,  25,   9, 1, 16};
        final double[] x3 = {8, 125, -27, 1, 64};

        MathArrays.sortInPlace(x1,
                               MathArrays.OrderDirection.DECREASING,
                               x2, x3);

        assertEquals(-3,  x1[4], FastMath.ulp(1d));
        assertEquals(9,   x2[4], FastMath.ulp(1d));
        assertEquals(-27, x3[4], FastMath.ulp(1d));

        assertEquals(1, x1[3], FastMath.ulp(1d));
        assertEquals(1, x2[3], FastMath.ulp(1d));
        assertEquals(1, x3[3], FastMath.ulp(1d));

        assertEquals(2, x1[2], FastMath.ulp(1d));
        assertEquals(4, x2[2], FastMath.ulp(1d));
        assertEquals(8, x3[2], FastMath.ulp(1d));

        assertEquals(4,  x1[1], FastMath.ulp(1d));
        assertEquals(16, x2[1], FastMath.ulp(1d));
        assertEquals(64, x3[1], FastMath.ulp(1d));

        assertEquals(5,   x1[0], FastMath.ulp(1d));
        assertEquals(25,  x2[0], FastMath.ulp(1d));
        assertEquals(125, x3[0], FastMath.ulp(1d));
    }

    /** Example in javadoc */
    @Test
    void testSortInPlaceExample() {
        final double[] x = {3, 1, 2};
        final double[] y = {1, 2, 3};
        final double[] z = {0, 5, 7};
        MathArrays.sortInPlace(x, y, z);
        final double[] sx = {1, 2, 3};
        final double[] sy = {2, 3, 1};
        final double[] sz = {5, 7, 0};
        assertArrayEquals(sx, x);
        assertArrayEquals(sy, y);
        assertArrayEquals(sz, z);
    }

    @Test
    void testSortInPlaceFailures() {
        final double[] nullArray = null;
        final double[] one = {1};
        final double[] two = {1, 2};
        final double[] onep = {2};
        try {
            MathArrays.sortInPlace(one, two);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MathArrays.sortInPlace(one, nullArray);
            fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
        try {
            MathArrays.sortInPlace(one, onep, nullArray);
            fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

    // MATH-1005
    @Test
    void testLinearCombinationWithSingleElementArray() {
        final double[] a = { 1.23456789 };
        final double[] b = { 98765432.1 };

        assertEquals(a[0] * b[0], MathArrays.linearCombination(a, b), 0d);
    }

    @Test
    void testLinearCombination1() {
        final double[] a = new double[] {
            -1321008684645961.0 / 268435456.0,
            -5774608829631843.0 / 268435456.0,
            -7645843051051357.0 / 8589934592.0
        };
        final double[] b = new double[] {
            -5712344449280879.0 / 2097152.0,
            -4550117129121957.0 / 2097152.0,
            8846951984510141.0 / 131072.0
        };

        final double abSumInline = MathArrays.linearCombination(a[0], b[0],
                                                                a[1], b[1],
                                                                a[2], b[2]);
        final double abSumArray = MathArrays.linearCombination(a, b);

        assertEquals(abSumInline, abSumArray, 0);
        assertEquals(-1.8551294182586248737720779899, abSumInline, 1.0e-15);

        final double naive = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        assertTrue(FastMath.abs(naive - abSumInline) > 1.5);

    }

    @Test
    void testLinearCombinationSignedZero() {
        assertTrue(FastMath.copySign(1, MathArrays.linearCombination(+0.0, 1.0, +0.0, 1.0)) > 0);
        assertTrue(FastMath.copySign(1, MathArrays.linearCombination(-0.0, 1.0, -0.0, 1.0)) < 0);
        assertTrue(FastMath.copySign(1, MathArrays.linearCombination(+0.0, 1.0, +0.0, 1.0, +0.0, 1.0)) > 0);
        assertTrue(FastMath.copySign(1, MathArrays.linearCombination(-0.0, 1.0, -0.0, 1.0, -0.0, 1.0)) < 0);
        assertTrue(FastMath.copySign(1, MathArrays.linearCombination(+0.0, 1.0, +0.0, 1.0, +0.0, 1.0)) > 0);
        assertTrue(FastMath.copySign(1, MathArrays.linearCombination(-0.0, 1.0, -0.0, 1.0, -0.0, 1.0)) < 0);
        assertTrue(FastMath.copySign(1, MathArrays.linearCombination(+0.0, 1.0, +0.0, 1.0, +0.0, 1.0, +0.0, 1.0)) > 0);
        assertTrue(FastMath.copySign(1, MathArrays.linearCombination(-0.0, 1.0, -0.0, 1.0, -0.0, 1.0, -0.0, 1.0)) < 0);
        assertTrue(FastMath.copySign(1, MathArrays.linearCombination(new double[] {+0.0, +0.0}, new double[] {1.0, 1.0})) > 0);
        assertTrue(FastMath.copySign(1, MathArrays.linearCombination(new double[] {-0.0, -0.0}, new double[] {1.0, 1.0})) < 0);
    }

    @Test
    void testLinearCombination2() {
        // we compare accurate versus naive dot product implementations
        // on regular vectors (i.e. not extreme cases like in the previous test)
        Well1024a random = new Well1024a(553267312521321234l);

        for (int i = 0; i < 10000; ++i) {
            final double ux = 1e17 * random.nextDouble();
            final double uy = 1e17 * random.nextDouble();
            final double uz = 1e17 * random.nextDouble();
            final double vx = 1e17 * random.nextDouble();
            final double vy = 1e17 * random.nextDouble();
            final double vz = 1e17 * random.nextDouble();
            final double sInline = MathArrays.linearCombination(ux, vx,
                                                                uy, vy,
                                                                uz, vz);
            final double sArray = MathArrays.linearCombination(new double[] {ux, uy, uz},
                                                               new double[] {vx, vy, vz});
            assertEquals(sInline, sArray, 0);
        }
    }

    @Test
    void testLinearCombinationHuge() {
        int scale = 971;
        final double[] a = new double[] {
                                         -1321008684645961.0 / 268435456.0,
                                         -5774608829631843.0 / 268435456.0,
                                         -7645843051051357.0 / 8589934592.0
                                     };
        final double[] b = new double[] {
                                         -5712344449280879.0 / 2097152.0,
                                         -4550117129121957.0 / 2097152.0,
                                          8846951984510141.0 / 131072.0
                                     };

        double[] scaledA = new double[a.length];
        double[] scaledB = new double[b.length];
        for (int i = 0; i < scaledA.length; ++i) {
            scaledA[i] = FastMath.scalb(a[i], -scale);
            scaledB[i] = FastMath.scalb(b[i], +scale);
        }
        final double abSumInline = MathArrays.linearCombination(scaledA[0], scaledB[0],
                                                                scaledA[1], scaledB[1],
                                                                scaledA[2], scaledB[2]);
        final double abSumArray = MathArrays.linearCombination(scaledA, scaledB);

        assertEquals(abSumInline, abSumArray, 0);
        assertEquals(-1.8551294182586248737720779899, abSumInline, 1.0e-15);

        final double naive = scaledA[0] * scaledB[0] + scaledA[1] * scaledB[1] + scaledA[2] * scaledB[2];
        assertTrue(FastMath.abs(naive - abSumInline) > 1.5);

    }

    @Test
    void testLinearCombinationInfinite() {
        final double[][] a = new double[][] {
            { 1, 2, 3, 4},
            { 1, Double.POSITIVE_INFINITY, 3, 4},
            { 1, 2, Double.POSITIVE_INFINITY, 4},
            { 1, Double.POSITIVE_INFINITY, 3, Double.NEGATIVE_INFINITY},
            { 1, 2, 3, 4},
            { 1, 2, 3, 4},
            { 1, 2, 3, 4},
            { 1, 2, 3, 4}
        };
        final double[][] b = new double[][] {
            { 1, -2, 3, 4},
            { 1, -2, 3, 4},
            { 1, -2, 3, 4},
            { 1, -2, 3, 4},
            { 1, Double.POSITIVE_INFINITY, 3, 4},
            { 1, -2, Double.POSITIVE_INFINITY, 4},
            { 1, Double.POSITIVE_INFINITY, 3, Double.NEGATIVE_INFINITY},
            { Double.NaN, -2, 3, 4}
        };

        assertEquals(-3,
                            MathArrays.linearCombination(a[0][0], b[0][0],
                                                         a[0][1], b[0][1]),
                            1.0e-10);
        assertEquals(6,
                            MathArrays.linearCombination(a[0][0], b[0][0],
                                                         a[0][1], b[0][1],
                                                         a[0][2], b[0][2]),
                            1.0e-10);
        assertEquals(22,
                            MathArrays.linearCombination(a[0][0], b[0][0],
                                                         a[0][1], b[0][1],
                                                         a[0][2], b[0][2],
                                                         a[0][3], b[0][3]),
                            1.0e-10);
        assertEquals(22, MathArrays.linearCombination(a[0], b[0]), 1.0e-10);

        assertEquals(Double.NEGATIVE_INFINITY,
                            MathArrays.linearCombination(a[1][0], b[1][0],
                                                         a[1][1], b[1][1]),
                            1.0e-10);
        assertEquals(Double.NEGATIVE_INFINITY,
                            MathArrays.linearCombination(a[1][0], b[1][0],
                                                         a[1][1], b[1][1],
                                                         a[1][2], b[1][2]),
                            1.0e-10);
        assertEquals(Double.NEGATIVE_INFINITY,
                            MathArrays.linearCombination(a[1][0], b[1][0],
                                                         a[1][1], b[1][1],
                                                         a[1][2], b[1][2],
                                                         a[1][3], b[1][3]),
                            1.0e-10);
        assertEquals(Double.NEGATIVE_INFINITY, MathArrays.linearCombination(a[1], b[1]), 1.0e-10);

        assertEquals(-3,
                            MathArrays.linearCombination(a[2][0], b[2][0],
                                                         a[2][1], b[2][1]),
                            1.0e-10);
        assertEquals(Double.POSITIVE_INFINITY,
                            MathArrays.linearCombination(a[2][0], b[2][0],
                                                         a[2][1], b[2][1],
                                                         a[2][2], b[2][2]),
                            1.0e-10);
        assertEquals(Double.POSITIVE_INFINITY,
                            MathArrays.linearCombination(a[2][0], b[2][0],
                                                         a[2][1], b[2][1],
                                                         a[2][2], b[2][2],
                                                         a[2][3], b[2][3]),
                            1.0e-10);
        assertEquals(Double.POSITIVE_INFINITY, MathArrays.linearCombination(a[2], b[2]), 1.0e-10);

        assertEquals(Double.NEGATIVE_INFINITY,
                            MathArrays.linearCombination(a[3][0], b[3][0],
                                                         a[3][1], b[3][1]),
                            1.0e-10);
        assertEquals(Double.NEGATIVE_INFINITY,
                            MathArrays.linearCombination(a[3][0], b[3][0],
                                                         a[3][1], b[3][1],
                                                         a[3][2], b[3][2]),
                            1.0e-10);
        assertEquals(Double.NEGATIVE_INFINITY,
                            MathArrays.linearCombination(a[3][0], b[3][0],
                                                         a[3][1], b[3][1],
                                                         a[3][2], b[3][2],
                                                         a[3][3], b[3][3]),
                            1.0e-10);
        assertEquals(Double.NEGATIVE_INFINITY, MathArrays.linearCombination(a[3], b[3]), 1.0e-10);

        assertEquals(Double.POSITIVE_INFINITY,
                            MathArrays.linearCombination(a[4][0], b[4][0],
                                                         a[4][1], b[4][1]),
                            1.0e-10);
        assertEquals(Double.POSITIVE_INFINITY,
                            MathArrays.linearCombination(a[4][0], b[4][0],
                                                         a[4][1], b[4][1],
                                                         a[4][2], b[4][2]),
                            1.0e-10);
        assertEquals(Double.POSITIVE_INFINITY,
                            MathArrays.linearCombination(a[4][0], b[4][0],
                                                         a[4][1], b[4][1],
                                                         a[4][2], b[4][2],
                                                         a[4][3], b[4][3]),
                            1.0e-10);
        assertEquals(Double.POSITIVE_INFINITY, MathArrays.linearCombination(a[4], b[4]), 1.0e-10);

        assertEquals(-3,
                            MathArrays.linearCombination(a[5][0], b[5][0],
                                                         a[5][1], b[5][1]),
                            1.0e-10);
        assertEquals(Double.POSITIVE_INFINITY,
                            MathArrays.linearCombination(a[5][0], b[5][0],
                                                         a[5][1], b[5][1],
                                                         a[5][2], b[5][2]),
                            1.0e-10);
        assertEquals(Double.POSITIVE_INFINITY,
                            MathArrays.linearCombination(a[5][0], b[5][0],
                                                         a[5][1], b[5][1],
                                                         a[5][2], b[5][2],
                                                         a[5][3], b[5][3]),
                            1.0e-10);
        assertEquals(Double.POSITIVE_INFINITY, MathArrays.linearCombination(a[5], b[5]), 1.0e-10);

        assertEquals(Double.POSITIVE_INFINITY,
                            MathArrays.linearCombination(a[6][0], b[6][0],
                                                         a[6][1], b[6][1]),
                            1.0e-10);
        assertEquals(Double.POSITIVE_INFINITY,
                            MathArrays.linearCombination(a[6][0], b[6][0],
                                                         a[6][1], b[6][1],
                                                         a[6][2], b[6][2]),
                            1.0e-10);
        assertTrue(Double.isNaN(MathArrays.linearCombination(a[6][0], b[6][0],
                                                                    a[6][1], b[6][1],
                                                                    a[6][2], b[6][2],
                                                                    a[6][3], b[6][3])));
        assertTrue(Double.isNaN(MathArrays.linearCombination(a[6], b[6])));

        assertTrue(Double.isNaN(MathArrays.linearCombination(a[7][0], b[7][0],
                                                                    a[7][1], b[7][1])));
        assertTrue(Double.isNaN(MathArrays.linearCombination(a[7][0], b[7][0],
                                                                    a[7][1], b[7][1],
                                                                    a[7][2], b[7][2])));
        assertTrue(Double.isNaN(MathArrays.linearCombination(a[7][0], b[7][0],
                                                                    a[7][1], b[7][1],
                                                                    a[7][2], b[7][2],
                                                                    a[7][3], b[7][3])));
        assertTrue(Double.isNaN(MathArrays.linearCombination(a[7], b[7])));
    }

    @Test
    void testArrayEquals() {
        assertFalse(MathArrays.equals(new double[] { 1d }, null));
        assertFalse(MathArrays.equals(null, new double[] { 1d }));
        assertTrue(MathArrays.equals((double[]) null, (double[]) null));

        assertFalse(MathArrays.equals(new double[] { 1d }, new double[0]));
        assertTrue(MathArrays.equals(new double[] { 1d }, new double[] { 1d }));
        assertTrue(MathArrays.equals(new double[] { Double.POSITIVE_INFINITY,
                                                           Double.NEGATIVE_INFINITY, 1d, 0d },
                                            new double[] { Double.POSITIVE_INFINITY,
                                                           Double.NEGATIVE_INFINITY, 1d, 0d }));
        assertFalse(MathArrays.equals(new double[] { Double.NaN },
                                             new double[] { Double.NaN }));
        assertFalse(MathArrays.equals(new double[] { Double.POSITIVE_INFINITY },
                                             new double[] { Double.NEGATIVE_INFINITY }));
        assertFalse(MathArrays.equals(new double[] { 1d },
                                             new double[] { FastMath.nextAfter(FastMath.nextAfter(1d, 2d), 2d) }));

    }

    @Test
    void testArrayEqualsIncludingNaN() {
        assertFalse(MathArrays.equalsIncludingNaN(new double[] { 1d }, null));
        assertFalse(MathArrays.equalsIncludingNaN(null, new double[] { 1d }));
        assertTrue(MathArrays.equalsIncludingNaN((double[]) null, (double[]) null));

        assertFalse(MathArrays.equalsIncludingNaN(new double[] { 1d }, new double[0]));
        assertTrue(MathArrays.equalsIncludingNaN(new double[] { 1d }, new double[] { 1d }));
        assertTrue(MathArrays.equalsIncludingNaN(new double[] { Double.NaN, Double.POSITIVE_INFINITY,
                                                                       Double.NEGATIVE_INFINITY, 1d, 0d },
                                                        new double[] { Double.NaN, Double.POSITIVE_INFINITY,
                                                                       Double.NEGATIVE_INFINITY, 1d, 0d }));
        assertFalse(MathArrays.equalsIncludingNaN(new double[] { Double.POSITIVE_INFINITY },
                                                         new double[] { Double.NEGATIVE_INFINITY }));
        assertFalse(MathArrays.equalsIncludingNaN(new double[] { 1d },
                                                         new double[] { FastMath.nextAfter(FastMath.nextAfter(1d, 2d), 2d) }));
    }

    @Test
    void testLongArrayEquals() {
        assertFalse(MathArrays.equals(new long[] { 1L }, null));
        assertFalse(MathArrays.equals(null, new long[] { 1L }));
        assertTrue(MathArrays.equals((long[]) null, (long[]) null));

        assertFalse(MathArrays.equals(new long[] { 1L }, new long[0]));
        assertTrue(MathArrays.equals(new long[] { 1L }, new long[] { 1L }));
        assertTrue(MathArrays.equals(new long[] { 0x7ff0000000000000L,
                                                           0xfff0000000000000L, 1L, 0L },
                                            new long[] { 0x7ff0000000000000L,
                                                           0xfff0000000000000L, 1L, 0L }));
        assertFalse(MathArrays.equals(new long[] { 0x7ff0000000000000L },
                                             new long[] { 0xfff0000000000000L }));

    }

    @Test
    void testIntArrayEquals() {
        assertFalse(MathArrays.equals(new int[] { 1 }, null));
        assertFalse(MathArrays.equals(null, new int[] { 1 }));
        assertTrue(MathArrays.equals((int[]) null, (int[]) null));

        assertFalse(MathArrays.equals(new int[] { 1 }, new int[0]));
        assertTrue(MathArrays.equals(new int[] { 1 }, new int[] { 1 }));
        assertTrue(MathArrays.equals(new int[] { Integer.MAX_VALUE,
                                                        Integer.MIN_VALUE, 1, 0 },
                                            new int[] { Integer.MAX_VALUE,
                                                        Integer.MIN_VALUE, 1, 0 }));
        assertFalse(MathArrays.equals(new int[] { Integer.MAX_VALUE },
                                             new int[] { Integer.MIN_VALUE }));

    }

    @Test
    void testByteArrayEquals() {
        assertFalse(MathArrays.equals(new byte[] { 1 }, null));
        assertFalse(MathArrays.equals(null, new byte[] { 1 }));
        assertTrue(MathArrays.equals((byte[]) null, (byte[]) null));

        assertFalse(MathArrays.equals(new byte[] { 1 }, new byte[0]));
        assertTrue(MathArrays.equals(new byte[] { 1 }, new byte[] { 1 }));
        assertTrue(MathArrays.equals(new byte[] { Byte.MAX_VALUE,
                                                         Byte.MIN_VALUE, 1, 0 },
                                            new byte[] { Byte.MAX_VALUE,
                                                         Byte.MIN_VALUE, 1, 0 }));
        assertFalse(MathArrays.equals(new byte[] { Byte.MAX_VALUE },
                                             new byte[] { Byte.MIN_VALUE }));

    }

    @Test
    void testShortArrayEquals() {
        assertFalse(MathArrays.equals(new short[] { 1 }, null));
        assertFalse(MathArrays.equals(null, new short[] { 1 }));
        assertTrue(MathArrays.equals((short[]) null, (short[]) null));

        assertFalse(MathArrays.equals(new short[] { 1 }, new short[0]));
        assertTrue(MathArrays.equals(new short[] { 1 }, new short[] { 1 }));
        assertTrue(MathArrays.equals(new short[] { Short.MAX_VALUE,
                                                          Short.MIN_VALUE, 1, 0 },
                                            new short[] { Short.MAX_VALUE,
                                                          Short.MIN_VALUE, 1, 0 }));
        assertFalse(MathArrays.equals(new short[] { Short.MAX_VALUE },
                                             new short[] { Short.MIN_VALUE }));

    }

    @Test
    void testFloatArrayEquals() {
        assertFalse(MathArrays.equals(new float[] { 1 }, null));
        assertFalse(MathArrays.equals(null, new float[] { 1 }));
        assertTrue(MathArrays.equals((float[]) null, (float[]) null));

        assertFalse(MathArrays.equals(new float[] { 1 }, new float[0]));
        assertTrue(MathArrays.equals(new float[] { 1 }, new float[] { 1 }));
        assertTrue(MathArrays.equals(new float[] { Float.MAX_VALUE,
                                                          Float.MIN_VALUE, 1, 0 },
                                            new float[] { Float.MAX_VALUE,
                                                          Float.MIN_VALUE, 1, 0 }));
        assertFalse(MathArrays.equals(new float[] { Float.MAX_VALUE },
                                             new float[] { Float.MIN_VALUE }));

    }

    @Test
    void testDoubleArrayEquals() {
        assertFalse(MathArrays.equals(new double[] { 1 }, null));
        assertFalse(MathArrays.equals(null, new double[] { 1 }));
        assertTrue(MathArrays.equals((double[]) null, (double[]) null));

        assertFalse(MathArrays.equals(new double[] { 1 }, new double[0]));
        assertTrue(MathArrays.equals(new double[] { 1 }, new double[] { 1 }));
        assertTrue(MathArrays.equals(new double[] { Double.MAX_VALUE,
                                                          Double.MIN_VALUE, 1, 0 },
                                            new double[] { Double.MAX_VALUE,
                                                          Double.MIN_VALUE, 1, 0 }));
        assertFalse(MathArrays.equals(new double[] { Double.MAX_VALUE },
                                             new double[] { Double.MIN_VALUE }));

    }

    @Test
    void testFieldArrayEquals() {
        assertFalse(MathArrays.equals(new Complex[] { Complex.I }, null));
        assertFalse(MathArrays.equals(null, new Complex[] { Complex.I }));
        assertTrue(MathArrays.equals((Complex[]) null, (Complex[]) null));

        assertFalse(MathArrays.equals(new Complex[] { Complex.I }, new Complex[0]));
        assertTrue(MathArrays.equals(new Complex[] { Complex.I }, new Complex[] { Complex.I }));
        assertTrue(MathArrays.equals(new Complex[] { new Complex(Double.MAX_VALUE, 0),
                                                     new Complex(Double.MIN_VALUE, 0),
                                                     Complex.ONE, Complex.ZERO },
                                     new Complex[] { new Complex(Double.MAX_VALUE, 0),
                                                     new Complex(Double.MIN_VALUE, 0),
                                                     Complex.ONE, Complex.ZERO }));
        assertFalse(MathArrays.equals(new Complex[] { new Complex(Double.MAX_VALUE, 0) },
                                             new Complex[] { new Complex(Double.MIN_VALUE, 0) }));

    }

    @Test
    void testNormalizeArray() {
        double[] testValues1 = new double[] {1, 1, 2};
        UnitTestUtils.customAssertEquals(new double[] { .25, .25, .5},
                                         MathArrays.normalizeArray(testValues1, 1),
                                         Double.MIN_VALUE);

        double[] testValues2 = new double[] {-1, -1, 1};
        UnitTestUtils.customAssertEquals(new double[] { 1, 1, -1},
                                         MathArrays.normalizeArray(testValues2, 1),
                                         Double.MIN_VALUE);

        // Ignore NaNs
        double[] testValues3 = new double[] {-1, -1, Double.NaN, 1, Double.NaN};
        UnitTestUtils.customAssertEquals(new double[] { 1, 1, Double.NaN, -1, Double.NaN},
                                         MathArrays.normalizeArray(testValues3, 1),
                                         Double.MIN_VALUE);

        // Zero sum -> MathRuntimeException
        double[] zeroSum = new double[] {-1, 1};
        try {
            MathArrays.normalizeArray(zeroSum, 1);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        // Infinite elements -> MathRuntimeException
        double[] hasInf = new double[] {1, 2, 1, Double.NEGATIVE_INFINITY};
        try {
            MathArrays.normalizeArray(hasInf, 1);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {}

        // Infinite target -> MathIllegalArgumentException
        try {
            MathArrays.normalizeArray(testValues1, Double.POSITIVE_INFINITY);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {}

        // NaN target -> MathIllegalArgumentException
        try {
            MathArrays.normalizeArray(testValues1, Double.NaN);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {}
    }

    @Test
    void testConvolve() {
        /* Test Case (obtained via SciPy)
         * x=[1.2,-1.8,1.4]
         * h=[1,0.8,0.5,0.3]
         * convolve(x,h) -> array([ 1.2 , -0.84,  0.56,  0.58,  0.16,  0.42])
         */
        double[] x1 = { 1.2, -1.8, 1.4 };
        double[] h1 = { 1, 0.8, 0.5, 0.3 };
        double[] y1 = { 1.2, -0.84, 0.56, 0.58, 0.16, 0.42 };
        double tolerance = 1e-13;

        double[] yActual = MathArrays.convolve(x1, h1);
        assertArrayEquals(y1, yActual, tolerance);

        double[] x2 = { 1, 2, 3 };
        double[] h2 = { 0, 1, 0.5 };
        double[] y2 = { 0, 1, 2.5, 4, 1.5 };

        yActual = MathArrays.convolve(x2, h2);
        assertArrayEquals(y2, yActual, tolerance);

        try {
            MathArrays.convolve(new double[]{1, 2}, null);
            fail("an exception should have been thrown");
        } catch (NullArgumentException e) {
            // expected behavior
        }

        try {
            MathArrays.convolve(null, new double[]{1, 2});
            fail("an exception should have been thrown");
        } catch (NullArgumentException e) {
            // expected behavior
        }

        try {
            MathArrays.convolve(new double[]{1, 2}, new double[]{});
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // expected behavior
        }

        try {
            MathArrays.convolve(new double[]{}, new double[]{1, 2});
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // expected behavior
        }

        try {
            MathArrays.convolve(new double[]{}, new double[]{});
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // expected behavior
        }
    }

    @Test
    void testShuffleTail() {
        final int[] orig = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        final int[] list = orig.clone();
        final int start = 4;
        MathArrays.shuffle(list, start, MathArrays.Position.TAIL, new Well1024a(7654321L));

        // Ensure that all entries below index "start" did not move.
        for (int i = 0; i < start; i++) {
            assertEquals(orig[i], list[i]);
        }

        // Ensure that at least one entry has moved.
        boolean ok = false;
        for (int i = start; i < orig.length - 1; i++) {
            if (orig[i] != list[i]) {
                ok = true;
                break;
            }
        }
        assertTrue(ok);
    }

    @Test
    void testShuffleHead() {
        final int[] orig = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        final int[] list = orig.clone();
        final int start = 4;
        MathArrays.shuffle(list, start, MathArrays.Position.HEAD, new Well1024a(1234567L));

        // Ensure that all entries above index "start" did not move.
        for (int i = start + 1; i < orig.length; i++) {
            assertEquals(orig[i], list[i]);
        }

        // Ensure that at least one entry has moved.
        boolean ok = false;
        for (int i = 0; i <= start; i++) {
            if (orig[i] != list[i]) {
                ok = true;
                break;
            }
        }
        assertTrue(ok);
    }

    @Test
    void testNatural() {
        final int n = 4;
        final int[] expected = {0, 1, 2, 3};

        final int[] natural = MathArrays.natural(n);
        for (int i = 0; i < n; i++) {
            assertEquals(expected[i], natural[i]);
        }
    }

    @Test
    void testNaturalZero() {
        final int[] natural = MathArrays.natural(0);
        assertEquals(0, natural.length);
    }

    @Test
    void testSequence() {
        final int size = 4;
        final int start = 5;
        final int stride = 2;
        final int[] expected = {5, 7, 9, 11};

        final int[] seq = MathArrays.sequence(size, start, stride);
        for (int i = 0; i < size; i++) {
            assertEquals(expected[i], seq[i]);
        }
    }

    @Test
    void testSequenceZero() {
        final int[] seq = MathArrays.sequence(0, 12345, 6789);
        assertEquals(0, seq.length);
    }

    @Test
    void testVerifyValuesPositive() {
        for (int j = 0; j < 6; j++) {
            for (int i = 1; i < (7 - j); i++) {
                assertTrue(MathArrays.verifyValues(testArray, 0, i));
            }
        }
        assertTrue(MathArrays.verifyValues(singletonArray, 0, 1));
        assertTrue(MathArrays.verifyValues(singletonArray, 0, 0, true));
    }

    @Test
    void testVerifyValuesNegative() {
        assertFalse(MathArrays.verifyValues(singletonArray, 0, 0));
        assertFalse(MathArrays.verifyValues(testArray, 0, 0));
        try {
            MathArrays.verifyValues(singletonArray, 2, 1);  // start past end
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MathArrays.verifyValues(testArray, 0, 7);  // end past end
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MathArrays.verifyValues(testArray, -1, 1);  // start negative
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MathArrays.verifyValues(testArray, 0, -1);  // length negative
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MathArrays.verifyValues(nullArray, 0, 1);  // null array
            fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
        try {
            MathArrays.verifyValues(testArray, nullArray, 0, 1);  // null weights array
            fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
        try {
            MathArrays.verifyValues(singletonArray, testWeightsArray, 0, 1);  // weights.length != value.length
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MathArrays.verifyValues(testArray, testNegativeWeightsArray, 0, 6);  // can't have negative weights
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    void testConcatenate() {
        final double[] u = new double[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        final double[] x = new double[] {0, 1, 2};
        final double[] y = new double[] {3, 4, 5, 6, 7, 8};
        final double[] z = new double[] {9};
        assertArrayEquals(u, MathArrays.concatenate(x, y, z), 0);
    }

    @Test
    void testConcatentateSingle() {
        final double[] x = new double[] {0, 1, 2};
        assertArrayEquals(x, MathArrays.concatenate(x), 0);
    }

    public void testConcatenateEmptyArguments() {
        final double[] x = new double[] {0, 1, 2};
        final double[] y = new double[] {3};
        final double[] z = new double[] {};
        final double[] u = new double[] {0, 1, 2, 3};
        assertArrayEquals(u,  MathArrays.concatenate(x, z, y), 0);
        assertArrayEquals(u,  MathArrays.concatenate(x, y, z), 0);
        assertArrayEquals(u,  MathArrays.concatenate(z, x, y), 0);
        assertEquals(0,  MathArrays.concatenate(z, z, z).length);
    }

    @Test
    void testConcatenateNullArguments() {
        assertThrows(NullPointerException.class, () -> {
            final double[] x = new double[]{0, 1, 2};
            MathArrays.concatenate(x, null);
        });
    }

    @Test
    void testUnique() {
        final double[] x = {0, 9, 3, 0, 11, 7, 3, 5, -1, -2};
        final double[] values = {11, 9, 7, 5, 3, 0, -1, -2};
        assertArrayEquals(values, MathArrays.unique(x), 0);
    }

    @Test
    void testUniqueInfiniteValues() {
        final double [] x = {0, Double.NEGATIVE_INFINITY, 3, Double.NEGATIVE_INFINITY,
            3, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        final double[] u = {Double.POSITIVE_INFINITY, 3, 0, Double.NEGATIVE_INFINITY};
        assertArrayEquals(u , MathArrays.unique(x), 0);
    }

    @Test
    void testUniqueNaNValues() {
        final double[] x = new double[] {10, 2, Double.NaN, Double.NaN, Double.NaN,
            Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY};
        final double[] u = MathArrays.unique(x);
        assertEquals(5, u.length);
        assertTrue(Double.isNaN(u[0]));
        assertEquals(Double.POSITIVE_INFINITY, u[1], 0);
        assertEquals(10, u[2], 0);
        assertEquals(2, u[3], 0);
        assertEquals(Double.NEGATIVE_INFINITY, u[4], 0);
    }

    @Test
    void testUniqueNullArgument() {
        assertThrows(NullPointerException.class, () -> MathArrays.unique(null));
    }

    @Test
    void testBuildArray1() {
        Binary64Field field = Binary64Field.getInstance();
        Binary64[] array = MathArrays.buildArray(field, 3);
        assertEquals(3, array.length);
        for (Binary64 d : array) {
            assertSame(field.getZero(), d);
        }
    }

    @Test
    void testBuildArray2AllIndices() {
        Binary64Field field = Binary64Field.getInstance();
        Binary64[][] array = MathArrays.buildArray(field, 3, 2);
        assertEquals(3, array.length);
        for (Binary64[] a1 : array) {
            assertEquals(2, a1.length);
            for (Binary64 d : a1) {
                assertSame(field.getZero(), d);
            }
        }
    }

    @Test
    void testBuildArray2MissingLastIndex() {
        Binary64Field field = Binary64Field.getInstance();
        Binary64[][] array = MathArrays.buildArray(field, 3, -1);
        assertEquals(3, array.length);
        for (Binary64[] a1 : array) {
            assertNull(a1);
        }
    }

    @Test
    void testBuildArray3AllIndices() {
        Binary64Field field = Binary64Field.getInstance();
        Binary64[][][] array = MathArrays.buildArray(field, 3, 2, 4);
        assertEquals(3, array.length);
        for (Binary64[][] a1 : array) {
            assertEquals(2, a1.length);
            for (Binary64[] a2 : a1) {
                assertEquals(4, a2.length);
                for (Binary64 d : a2) {
                    assertSame(field.getZero(), d);
                }
            }
        }
    }

    @Test
    void testBuildArray3MissingLastIndex() {
        Binary64Field field = Binary64Field.getInstance();
        Binary64[][][] array = MathArrays.buildArray(field, 3, 2, -1);
        assertEquals(3, array.length);
        for (Binary64[][] a1 : array) {
            assertEquals(2, a1.length);
            for (Binary64[] a2 : a1) {
                assertNull(a2);
            }
        }
    }

}
