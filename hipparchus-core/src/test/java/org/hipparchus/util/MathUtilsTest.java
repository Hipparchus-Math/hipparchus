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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.FieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.util.MathUtils.FieldSumAndResidual;
import org.hipparchus.util.MathUtils.SumAndResidual;
import org.junit.jupiter.api.Test;

import java.util.regex.Pattern;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the MathUtils class.
 */
final class MathUtilsTest {
    @Test
    void testEqualsDouble() {
        final double x = 1234.5678;
        assertTrue(MathUtils.equals(x, x));
        assertFalse(MathUtils.equals(x, -x));

        // Special cases (cf. semantics of JDK's "Double").
        // 1. NaN
        assertTrue(MathUtils.equals(Double.NaN, Double.NaN));
        // 2. Negative zero
        final double mZero = -0d;
        final double zero = 0d;
        assertTrue(MathUtils.equals(zero, zero));
        assertTrue(MathUtils.equals(mZero, mZero));
        assertFalse(MathUtils.equals(mZero, zero));
    }

    @Test
    void testHash() {
        double[] testArray = {
            Double.NaN,
            Double.POSITIVE_INFINITY,
            Double.NEGATIVE_INFINITY,
            1d,
            0d,
            1E-14,
            (1 + 1E-14),
            Double.MIN_VALUE,
            Double.MAX_VALUE };
        for (int i = 0; i < testArray.length; i++) {
            for (int j = 0; j < testArray.length; j++) {
                if (i == j) {
                    assertEquals(MathUtils.hash(testArray[i]), MathUtils.hash(testArray[j]));
                    assertEquals(MathUtils.hash(testArray[j]), MathUtils.hash(testArray[i]));
                } else {
                    assertTrue(MathUtils.hash(testArray[i]) != MathUtils.hash(testArray[j]));
                    assertTrue(MathUtils.hash(testArray[j]) != MathUtils.hash(testArray[i]));
                }
            }
        }
    }

    @Test
    void testArrayHash() {
        assertEquals(0, MathUtils.hash(null));
        assertEquals(MathUtils.hash(new double[] {
                                      Double.NaN, Double.POSITIVE_INFINITY,
                                      Double.NEGATIVE_INFINITY, 1d, 0d
                                    }),
                     MathUtils.hash(new double[] {
                                      Double.NaN, Double.POSITIVE_INFINITY,
                                      Double.NEGATIVE_INFINITY, 1d, 0d
                                    }));
        assertNotEquals(MathUtils.hash(new double[]{1d}), MathUtils.hash(new double[]{FastMath.nextAfter(1d, 2d)}));
        assertNotEquals(MathUtils.hash(new double[]{1d}), MathUtils.hash(new double[]{1d, 1d}));
    }

    /**
     * Make sure that permuted arrays do not hash to the same value.
     */
    @Test
    void testPermutedArrayHash() {
        double[] original = new double[10];
        double[] permuted = new double[10];
        RandomDataGenerator random = new RandomDataGenerator(100);

        // Generate 10 distinct random values
        for (int i = 0; i < 10; i++) {
            original[i] = random.nextUniform(i + 0.5, i + 0.75);
        }

        // Generate a random permutation, making sure it is not the identity
        boolean isIdentity = true;
        do {
            int[] permutation = random.nextPermutation(10, 10);
            for (int i = 0; i < 10; i++) {
                if (i != permutation[i]) {
                    isIdentity = false;
                }
                permuted[i] = original[permutation[i]];
            }
        } while (isIdentity);

        // Verify that permuted array has different hash
        assertNotEquals(MathUtils.hash(original), MathUtils.hash(permuted));
    }

    @Test
    void testIndicatorByte() {
        assertEquals((byte)1, MathUtils.copySign((byte)1, (byte)2));
        assertEquals((byte)1, MathUtils.copySign((byte)1, (byte)0));
        assertEquals((byte)(-1), MathUtils.copySign((byte)1, (byte)(-2)));
    }

    @Test
    void testIndicatorInt() {
        assertEquals(1, MathUtils.copySign(1, 2));
        assertEquals(1, MathUtils.copySign(1, 0));
        assertEquals((-1), MathUtils.copySign(1, -2));
    }

    @Test
    void testIndicatorLong() {
        assertEquals(1L, MathUtils.copySign(1L, 2L));
        assertEquals(1L, MathUtils.copySign(1L, 0L));
        assertEquals(-1L, MathUtils.copySign(1L, -2L));
    }

    @Test
    void testIndicatorShort() {
        assertEquals((short)1, MathUtils.copySign((short)1, (short)2));
        assertEquals((short)1, MathUtils.copySign((short)1, (short)0));
        assertEquals((short)(-1), MathUtils.copySign((short)1, (short)(-2)));
    }

    @Test
    void testNormalizeAngle() {
        for (double a = -15.0; a <= 15.0; a += 0.1) {
            for (double b = -15.0; b <= 15.0; b += 0.2) {
                double c = MathUtils.normalizeAngle(a, b);
                assertTrue((b - FastMath.PI) <= c);
                assertTrue(c <= (b + FastMath.PI));
                double twoK = FastMath.rint((a - c) / FastMath.PI);
                assertEquals(c, a - twoK * FastMath.PI, 1.0e-14);
            }
        }
    }

    @Test
    void testFieldNormalizeAngle() {
        doTestFieldNormalizeAngle(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestFieldNormalizeAngle(final Field<T> field) {
        final T zero = field.getZero();
        for (double a = -15.0; a <= 15.0; a += 0.1) {
            for (double b = -15.0; b <= 15.0; b += 0.2) {
                T c = MathUtils.normalizeAngle(zero.add(a), zero.add(b));
                double cR = c.getReal();
                assertTrue((b - FastMath.PI) <= cR);
                assertTrue(cR <= (b + FastMath.PI));
                double twoK = FastMath.rint((a - cR) / FastMath.PI);
                assertEquals(cR, a - twoK * FastMath.PI, 1.0e-14);
            }
        }
    }

    @Test
    void testReduce() {
        final double period = -12.222;
        final double offset = 13;

        final double delta = 1.5;

        double orig = offset + 122456789 * period + delta;
        double expected = delta;
        assertEquals(expected,
                            MathUtils.reduce(orig, period, offset),
                            1e-7);
        assertEquals(expected,
                            MathUtils.reduce(orig, -period, offset),
                            1e-7);

        orig = offset - 123356789 * period - delta;
        expected = FastMath.abs(period) - delta;
        assertEquals(expected,
                            MathUtils.reduce(orig, period, offset),
                            1e-6);
        assertEquals(expected,
                            MathUtils.reduce(orig, -period, offset),
                            1e-6);

        orig = offset - 123446789 * period + delta;
        expected = delta;
        assertEquals(expected,
                            MathUtils.reduce(orig, period, offset),
                            1e-6);
        assertEquals(expected,
                            MathUtils.reduce(orig, -period, offset),
                            1e-6);

        assertTrue(Double.isNaN(MathUtils.reduce(orig, Double.NaN, offset)));
        assertTrue(Double.isNaN(MathUtils.reduce(Double.NaN, period, offset)));
        assertTrue(Double.isNaN(MathUtils.reduce(orig, period, Double.NaN)));
        assertTrue(Double.isNaN(MathUtils.reduce(orig, period,
                Double.POSITIVE_INFINITY)));
        assertTrue(Double.isNaN(MathUtils.reduce(Double.POSITIVE_INFINITY,
                period, offset)));
        assertTrue(Double.isNaN(MathUtils.reduce(orig,
                Double.POSITIVE_INFINITY, offset)));
        assertTrue(Double.isNaN(MathUtils.reduce(orig,
                Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)));
        assertTrue(Double.isNaN(MathUtils.reduce(Double.POSITIVE_INFINITY,
                period, Double.POSITIVE_INFINITY)));
        assertTrue(Double.isNaN(MathUtils.reduce(Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY, offset)));
        assertTrue(Double.isNaN(MathUtils.reduce(Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY,  Double.POSITIVE_INFINITY)));
    }

    @Test
    void testReduceComparedWithNormalizeAngle() {
        final double tol = Math.ulp(1d);
        final double period = 2 * Math.PI;
        for (double a = -15; a <= 15; a += 0.5) {
            for (double center = -15; center <= 15; center += 1) {
                final double nA = MathUtils.normalizeAngle(a, center);
                final double offset = center - Math.PI;
                final double r = MathUtils.reduce(a, period, offset);
                assertEquals(nA, r + offset, tol);
            }
        }
    }

    @Test
    void testSignByte() {
        final byte one = (byte) 1;
        assertEquals((byte) 1, MathUtils.copySign(one, (byte) 2));
        assertEquals((byte) (-1), MathUtils.copySign(one, (byte) (-2)));
    }

    @Test
    void testSignInt() {
        final int one = 1;
        assertEquals(1, MathUtils.copySign(one, 2));
        assertEquals((-1), MathUtils.copySign(one, -2));
    }

    @Test
    void testSignLong() {
        final long one = 1L;
        assertEquals(1L, MathUtils.copySign(one, 2L));
        assertEquals(-1L, MathUtils.copySign(one, -2L));
    }

    @Test
    void testSignShort() {
        final short one = (short) 1;
        assertEquals((short) 1, MathUtils.copySign(one, (short) 2));
        assertEquals((short) (-1), MathUtils.copySign(one, (short) (-2)));
    }

    @Test
    void testCheckFinite() {
        try {
            MathUtils.checkFinite(Double.POSITIVE_INFINITY);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            assertEquals(LocalizedCoreFormats.NOT_FINITE_NUMBER, e.getSpecifier());
        }
        try {
            MathUtils.checkFinite(Double.NEGATIVE_INFINITY);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            assertEquals(LocalizedCoreFormats.NOT_FINITE_NUMBER, e.getSpecifier());
        }
        try {
            MathUtils.checkFinite(Double.NaN);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            assertEquals(LocalizedCoreFormats.NOT_FINITE_NUMBER, e.getSpecifier());
        }

        try {
            MathUtils.checkFinite(new double[] {0, -1, Double.POSITIVE_INFINITY, -2, 3});
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            assertEquals(LocalizedCoreFormats.NOT_FINITE_NUMBER, e.getSpecifier());
        }
        try {
            MathUtils.checkFinite(new double[] {1, Double.NEGATIVE_INFINITY, -2, 3});
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            assertEquals(LocalizedCoreFormats.NOT_FINITE_NUMBER, e.getSpecifier());
        }
        try {
            MathUtils.checkFinite(new double[] {4, 3, -1, Double.NaN, -2, 1});
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            assertEquals(LocalizedCoreFormats.NOT_FINITE_NUMBER, e.getSpecifier());
        }
    }

    @Test
    void testCheckNotNull1() {
        try {
            Object obj = null;
            MathUtils.checkNotNull(obj);
        } catch (NullArgumentException e) {
            // Expected.
        }
    }

    @Test
    void testCheckNotNull2() {
        try {
            double[] array = null;
            MathUtils.checkNotNull(array, LocalizedCoreFormats.INPUT_ARRAY);
        } catch (NullArgumentException e) {
            // Expected.
        }
    }

    @Test
    void testCopySignByte() {
        byte a = MathUtils.copySign(Byte.MIN_VALUE, (byte) -1);
        assertEquals(Byte.MIN_VALUE, a);

        final byte minValuePlusOne = Byte.MIN_VALUE + (byte) 1;
        a = MathUtils.copySign(minValuePlusOne, (byte) 1);
        assertEquals(Byte.MAX_VALUE, a);

        a = MathUtils.copySign(Byte.MAX_VALUE, (byte) -1);
        assertEquals(minValuePlusOne, a);

        final byte one = 1;
        byte val = -2;
        a = MathUtils.copySign(val, one);
        assertEquals(-val, a);

        final byte minusOne = -one;
        val = 2;
        a = MathUtils.copySign(val, minusOne);
        assertEquals(-val, a);

        val = 0;
        a = MathUtils.copySign(val, minusOne);
        assertEquals(val, a);

        val = 0;
        a = MathUtils.copySign(val, one);
        assertEquals(val, a);
    }

    @Test
    void testCopySignByte2() {
        assertThrows(MathRuntimeException.class, () -> MathUtils.copySign(Byte.MIN_VALUE, (byte) 1));
    }

    @Test
    public void testHipparchusVersion() {
        final Pattern pattern = Pattern.compile("unknown|[0-9.]*(?:-SNAPSHOT)?");
        assertTrue(pattern.matcher(MathUtils.getHipparchusVersion()).matches());
    }

    /**
     * Tests {@link MathUtils#twoSum(double, double)}.
     */
    @Test
    void testTwoSum() {
        //      sum = 0.30000000000000004
        // residual = -2.7755575615628914E-17
        final double a1 = 0.1;
        final double b1 = 0.2;
        final SumAndResidual result1 = MathUtils.twoSum(a1, b1);
        assertEquals(a1 + b1, result1.getSum(), 0.);
        assertEquals(a1 + b1, result1.getSum() + result1.getResidual(), 0.);
        assertNotEquals(0., result1.getResidual(), 0.);

        //      sum = -1580.3205849419005
        // residual = -1.1368683772161603E-13
        final double a2 = -615.7212034581913;
        final double b2 = -964.5993814837093;
        final SumAndResidual result2 = MathUtils.twoSum(a2, b2);
        assertEquals(a2 + b2, result2.getSum(), 0.);
        assertEquals(a2 + b2, result2.getSum() + result2.getResidual(), 0.);
        assertNotEquals(0., result2.getResidual(), 0.);

        //      sum = 251.8625825973395
        // residual = 1.4210854715202004E-14
        final double a3 = 60.348375484313706;
        final double b3 = 191.5142071130258;
        final SumAndResidual result3 = MathUtils.twoSum(a3, b3);
        assertEquals(a3 + b3, result3.getSum(), 0.);
        assertEquals(a3 + b3, result3.getSum() + result3.getResidual(), 0.);
        assertNotEquals(0., result3.getResidual(), 0.);

        //      sum = 622.8319023175123
        // residual = -4.3315557304163255E-14
        final double a4 = 622.8314146170453;
        final double b4 = 0.0004877004669900762;
        final SumAndResidual result4 = MathUtils.twoSum(a4, b4);
        assertEquals(a4 + b4, result4.getSum(), 0.);
        assertEquals(a4 + b4, result4.getSum() + result4.getResidual(), 0.);
        assertNotEquals(0., result4.getResidual(), 0.);
    }

    /**
     * Tests {@link MathUtils#twoSum(FieldElement, FieldElement)}.
     */
    @Test
    void testTwoSumField() {
        final Tuple a = new Tuple(0.1, -615.7212034581913, 60.348375484313706, 622.8314146170453);
        final Tuple b = new Tuple(0.2, -964.5993814837093, 191.5142071130258, 0.0004877004669900762);
        final FieldSumAndResidual<Tuple> result = MathUtils.twoSum(a, b);
        for (int i = 0; i < a.getDimension(); ++i) {
            assertEquals(a.getComponent(i) + b.getComponent(i), result.getSum().getComponent(i), 0.);
            assertEquals(a.getComponent(i) + b.getComponent(i),
                    result.getSum().getComponent(i) + result.getResidual().getComponent(i), 0.);
            assertNotEquals(0., result.getResidual().getComponent(i), 0.);
        }
    }

}
