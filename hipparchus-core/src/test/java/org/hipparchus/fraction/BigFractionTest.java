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
package org.hipparchus.fraction;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;


class BigFractionTest {

    private void customCustomAssertFraction(int expectedNumerator, int expectedDenominator, BigFraction actual) {
        assertEquals(expectedNumerator, actual.getNumeratorAsInt());
        assertEquals(expectedDenominator, actual.getDenominatorAsInt());
    }

    private void customCustomAssertFraction(long expectedNumerator, long expectedDenominator, BigFraction actual) {
        assertEquals(expectedNumerator, actual.getNumeratorAsLong());
        assertEquals(expectedDenominator, actual.getDenominatorAsLong());
    }

    @Test
    void testConstructor() {
        customCustomAssertFraction(0, 1, new BigFraction(0, 1));
        customCustomAssertFraction(0, 1, new BigFraction(0l, 2l));
        customCustomAssertFraction(0, 1, new BigFraction(0, -1));
        customCustomAssertFraction(1, 2, new BigFraction(1, 2));
        customCustomAssertFraction(1, 2, new BigFraction(2, 4));
        customCustomAssertFraction(-1, 2, new BigFraction(-1, 2));
        customCustomAssertFraction(-1, 2, new BigFraction(1, -2));
        customCustomAssertFraction(-1, 2, new BigFraction(-2, 4));
        customCustomAssertFraction(-1, 2, new BigFraction(2, -4));
        customCustomAssertFraction(11, 1, new BigFraction(11));
        customCustomAssertFraction(11, 1, new BigFraction(11l));
        customCustomAssertFraction(11, 1, new BigFraction(new BigInteger("11")));

        customCustomAssertFraction(0, 1, new BigFraction(0.00000000000001, 1.0e-5, 100));
        customCustomAssertFraction(2, 5, new BigFraction(0.40000000000001, 1.0e-5, 100));
        customCustomAssertFraction(15, 1, new BigFraction(15.0000000000001, 1.0e-5, 100));

        assertEquals(0.00000000000001, new BigFraction(0.00000000000001).doubleValue(), 0.0);
        assertEquals(0.40000000000001, new BigFraction(0.40000000000001).doubleValue(), 0.0);
        assertEquals(15.0000000000001, new BigFraction(15.0000000000001).doubleValue(), 0.0);
        customCustomAssertFraction(3602879701896487l, 9007199254740992l, new BigFraction(0.40000000000001));
        customCustomAssertFraction(1055531162664967l, 70368744177664l, new BigFraction(15.0000000000001));
        try {
            new BigFraction(null, BigInteger.ONE);
            fail("Expecting NullArgumentException");
        } catch (NullArgumentException npe) {
            // expected
        }
        try {
            new BigFraction(BigInteger.ONE, null);
            fail("Expecting NullArgumentException");
        } catch (NullArgumentException npe) {
            // expected
        }
        try {
            new BigFraction(BigInteger.ONE, BigInteger.ZERO);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException npe) {
            // expected
        }
    }

    @Test
    void testIsInteger() {
        assertTrue(new BigFraction(12, 12).isInteger());
        assertTrue(new BigFraction(14, 7).isInteger());
        assertFalse(new BigFraction(12, 11).isInteger());
    }

    @Test
    void testGoldenRatio() {
        assertThrows(MathIllegalStateException.class, () -> {
            // the golden ratio is notoriously a difficult number for continuous fraction
            new BigFraction((1 + FastMath.sqrt(5)) / 2, 1.0e-12, 25);
        });
    }

    // MATH-179
    @Test
    void testDoubleConstructor() throws MathIllegalStateException {
        customCustomAssertFraction(1, 2, new BigFraction((double) 1 / (double) 2, 1.0e-5, 100));
        customCustomAssertFraction(1, 3, new BigFraction((double) 1 / (double) 3, 1.0e-5, 100));
        customCustomAssertFraction(2, 3, new BigFraction((double) 2 / (double) 3, 1.0e-5, 100));
        customCustomAssertFraction(1, 4, new BigFraction((double) 1 / (double) 4, 1.0e-5, 100));
        customCustomAssertFraction(3, 4, new BigFraction((double) 3 / (double) 4, 1.0e-5, 100));
        customCustomAssertFraction(1, 5, new BigFraction((double) 1 / (double) 5, 1.0e-5, 100));
        customCustomAssertFraction(2, 5, new BigFraction((double) 2 / (double) 5, 1.0e-5, 100));
        customCustomAssertFraction(3, 5, new BigFraction((double) 3 / (double) 5, 1.0e-5, 100));
        customCustomAssertFraction(4, 5, new BigFraction((double) 4 / (double) 5, 1.0e-5, 100));
        customCustomAssertFraction(1, 6, new BigFraction((double) 1 / (double) 6, 1.0e-5, 100));
        customCustomAssertFraction(5, 6, new BigFraction((double) 5 / (double) 6, 1.0e-5, 100));
        customCustomAssertFraction(1, 7, new BigFraction((double) 1 / (double) 7, 1.0e-5, 100));
        customCustomAssertFraction(2, 7, new BigFraction((double) 2 / (double) 7, 1.0e-5, 100));
        customCustomAssertFraction(3, 7, new BigFraction((double) 3 / (double) 7, 1.0e-5, 100));
        customCustomAssertFraction(4, 7, new BigFraction((double) 4 / (double) 7, 1.0e-5, 100));
        customCustomAssertFraction(5, 7, new BigFraction((double) 5 / (double) 7, 1.0e-5, 100));
        customCustomAssertFraction(6, 7, new BigFraction((double) 6 / (double) 7, 1.0e-5, 100));
        customCustomAssertFraction(1, 8, new BigFraction((double) 1 / (double) 8, 1.0e-5, 100));
        customCustomAssertFraction(3, 8, new BigFraction((double) 3 / (double) 8, 1.0e-5, 100));
        customCustomAssertFraction(5, 8, new BigFraction((double) 5 / (double) 8, 1.0e-5, 100));
        customCustomAssertFraction(7, 8, new BigFraction((double) 7 / (double) 8, 1.0e-5, 100));
        customCustomAssertFraction(1, 9, new BigFraction((double) 1 / (double) 9, 1.0e-5, 100));
        customCustomAssertFraction(2, 9, new BigFraction((double) 2 / (double) 9, 1.0e-5, 100));
        customCustomAssertFraction(4, 9, new BigFraction((double) 4 / (double) 9, 1.0e-5, 100));
        customCustomAssertFraction(5, 9, new BigFraction((double) 5 / (double) 9, 1.0e-5, 100));
        customCustomAssertFraction(7, 9, new BigFraction((double) 7 / (double) 9, 1.0e-5, 100));
        customCustomAssertFraction(8, 9, new BigFraction((double) 8 / (double) 9, 1.0e-5, 100));
        customCustomAssertFraction(1, 10, new BigFraction((double) 1 / (double) 10, 1.0e-5, 100));
        customCustomAssertFraction(3, 10, new BigFraction((double) 3 / (double) 10, 1.0e-5, 100));
        customCustomAssertFraction(7, 10, new BigFraction((double) 7 / (double) 10, 1.0e-5, 100));
        customCustomAssertFraction(9, 10, new BigFraction((double) 9 / (double) 10, 1.0e-5, 100));
        customCustomAssertFraction(1, 11, new BigFraction((double) 1 / (double) 11, 1.0e-5, 100));
        customCustomAssertFraction(2, 11, new BigFraction((double) 2 / (double) 11, 1.0e-5, 100));
        customCustomAssertFraction(3, 11, new BigFraction((double) 3 / (double) 11, 1.0e-5, 100));
        customCustomAssertFraction(4, 11, new BigFraction((double) 4 / (double) 11, 1.0e-5, 100));
        customCustomAssertFraction(5, 11, new BigFraction((double) 5 / (double) 11, 1.0e-5, 100));
        customCustomAssertFraction(6, 11, new BigFraction((double) 6 / (double) 11, 1.0e-5, 100));
        customCustomAssertFraction(7, 11, new BigFraction((double) 7 / (double) 11, 1.0e-5, 100));
        customCustomAssertFraction(8, 11, new BigFraction((double) 8 / (double) 11, 1.0e-5, 100));
        customCustomAssertFraction(9, 11, new BigFraction((double) 9 / (double) 11, 1.0e-5, 100));
        customCustomAssertFraction(10, 11, new BigFraction((double) 10 / (double) 11, 1.0e-5, 100));
    }

    // MATH-181
    @Test
    void testDigitLimitConstructor() {
        customCustomAssertFraction(2, 5, new BigFraction(0.4, 9));
        customCustomAssertFraction(2, 5, new BigFraction(0.4, 99));
        customCustomAssertFraction(2, 5, new BigFraction(0.4, 999));

        customCustomAssertFraction(3, 5, new BigFraction(0.6152, 9));
        customCustomAssertFraction(8, 13, new BigFraction(0.6152, 99));
        customCustomAssertFraction(510, 829, new BigFraction(0.6152, 999));
        customCustomAssertFraction(769, 1250, new BigFraction(0.6152, 9999));

        // MATH-996
        customCustomAssertFraction(1, 2, new BigFraction(0.5000000001, 10));
    }

    // MATH-1029
    @Test
    void testPositiveValueOverflow() {
        customCustomAssertFraction((long) 1e10, 1, new BigFraction(1e10, 1000));
    }

    // MATH-1029
    @Test
    void testNegativeValueOverflow() {
        customCustomAssertFraction((long) -1e10, 1, new BigFraction(-1e10, 1000));
    }

    @Test
    void testEpsilonLimitConstructor() {
        customCustomAssertFraction(2, 5, new BigFraction(0.4, 1.0e-5, 100));

        customCustomAssertFraction(3, 5, new BigFraction(0.6152, 0.02, 100));
        customCustomAssertFraction(8, 13, new BigFraction(0.6152, 1.0e-3, 100));
        customCustomAssertFraction(251, 408, new BigFraction(0.6152, 1.0e-4, 100));
        customCustomAssertFraction(251, 408, new BigFraction(0.6152, 1.0e-5, 100));
        customCustomAssertFraction(510, 829, new BigFraction(0.6152, 1.0e-6, 100));
        customCustomAssertFraction(769, 1250, new BigFraction(0.6152, 1.0e-7, 100));
    }

    @Test
    void testCompareTo() {
        BigFraction first = new BigFraction(1, 2);
        BigFraction second = new BigFraction(1, 3);
        BigFraction third = new BigFraction(1, 2);

        assertEquals(0, first.compareTo(first));
        assertEquals(0, first.compareTo(third));
        assertEquals(1, first.compareTo(second));
        assertEquals(-1, second.compareTo(first));

        // these two values are different approximations of PI
        // the first  one is approximately PI - 3.07e-18
        // the second one is approximately PI + 1.936e-17
        BigFraction pi1 = new BigFraction(1068966896, 340262731);
        BigFraction pi2 = new BigFraction( 411557987, 131002976);
        assertEquals(-1, pi1.compareTo(pi2));
        assertEquals( 1, pi2.compareTo(pi1));
        assertEquals(0.0, pi1.doubleValue() - pi2.doubleValue(), 1.0e-20);

    }

    @Test
    void testDoubleValue() {
        BigFraction first = new BigFraction(1, 2);
        BigFraction second = new BigFraction(1, 3);

        assertEquals(0.5, first.doubleValue(), 0.0);
        assertEquals(1.0 / 3.0, second.doubleValue(), 0.0);
    }

    // MATH-744
    @Test
    void testDoubleValueForLargeNumeratorAndDenominator() {
        final BigInteger pow400 = BigInteger.TEN.pow(400);
        final BigInteger pow401 = BigInteger.TEN.pow(401);
        final BigInteger two = new BigInteger("2");
        final BigFraction large = new BigFraction(pow401.add(BigInteger.ONE),
                                                  pow400.multiply(two));

        assertEquals(5, large.doubleValue(), 1e-15);
    }

    // MATH-744
    @Test
    void testFloatValueForLargeNumeratorAndDenominator() {
        final BigInteger pow400 = BigInteger.TEN.pow(400);
        final BigInteger pow401 = BigInteger.TEN.pow(401);
        final BigInteger two = new BigInteger("2");
        final BigFraction large = new BigFraction(pow401.add(BigInteger.ONE),
                                                  pow400.multiply(two));

        assertEquals(5, large.floatValue(), 1e-15);
    }

    @Test
    void testFloatValue() {
        BigFraction first = new BigFraction(1, 2);
        BigFraction second = new BigFraction(1, 3);

        assertEquals(0.5f, first.floatValue(), 0.0f);
        assertEquals((float) (1.0 / 3.0), second.floatValue(), 0.0f);
    }

    @Test
    void testIntValue() {
        BigFraction first = new BigFraction(1, 2);
        BigFraction second = new BigFraction(3, 2);

        assertEquals(0, first.intValue());
        assertEquals(1, second.intValue());
    }

    @Test
    void testLongValue() {
        BigFraction first = new BigFraction(1, 2);
        BigFraction second = new BigFraction(3, 2);

        assertEquals(0L, first.longValue());
        assertEquals(1L, second.longValue());
    }

    @Test
    void testConstructorDouble() {
        customCustomAssertFraction(1, 2, new BigFraction(0.5));
        customCustomAssertFraction(6004799503160661l, 18014398509481984l, new BigFraction(1.0 / 3.0));
        customCustomAssertFraction(6124895493223875l, 36028797018963968l, new BigFraction(17.0 / 100.0));
        customCustomAssertFraction(1784551352345559l, 562949953421312l, new BigFraction(317.0 / 100.0));
        customCustomAssertFraction(-1, 2, new BigFraction(-0.5));
        customCustomAssertFraction(-6004799503160661l, 18014398509481984l, new BigFraction(-1.0 / 3.0));
        customCustomAssertFraction(-6124895493223875l, 36028797018963968l, new BigFraction(17.0 / -100.0));
        customCustomAssertFraction(-1784551352345559l, 562949953421312l, new BigFraction(-317.0 / 100.0));
        for (double v : new double[] { Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}) {
            try {
                new BigFraction(v);
                fail("Expecting MathIllegalArgumentException");
            } catch (MathIllegalArgumentException iae) {
                // expected
            }
        }
        assertEquals(1l, new BigFraction(Double.MAX_VALUE).getDenominatorAsLong());
        assertEquals(1l, new BigFraction(Double.longBitsToDouble(0x0010000000000000L)).getNumeratorAsLong());
        assertEquals(1l, new BigFraction(Double.MIN_VALUE).getNumeratorAsLong());
    }

    @Test
    void testAbs() {
        BigFraction a = new BigFraction(10, 21);
        BigFraction b = new BigFraction(-10, 21);
        BigFraction c = new BigFraction(10, -21);

        customCustomAssertFraction(10, 21, a.abs());
        customCustomAssertFraction(10, 21, b.abs());
        customCustomAssertFraction(10, 21, c.abs());
    }

    @Test
    void testSignum() {
        assertEquals(-1, new BigFraction(4, -5).signum());
        assertEquals(-1, new BigFraction(-4, 5).signum());
        assertEquals( 0, new BigFraction(0).signum());
        assertEquals(+1, new BigFraction(-4, -5).signum());
        assertEquals(+1, new BigFraction(4, 5).signum());
    }

    @Test
    void testReciprocal() {
        BigFraction f = null;

        f = new BigFraction(50, 75);
        f = f.reciprocal();
        assertEquals(3, f.getNumeratorAsInt());
        assertEquals(2, f.getDenominatorAsInt());

        f = new BigFraction(4, 3);
        f = f.reciprocal();
        assertEquals(3, f.getNumeratorAsInt());
        assertEquals(4, f.getDenominatorAsInt());

        f = new BigFraction(-15, 47);
        f = f.reciprocal();
        assertEquals(-47, f.getNumeratorAsInt());
        assertEquals(15, f.getDenominatorAsInt());

        f = new BigFraction(0, 3);
        try {
            f = f.reciprocal();
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
        }

        // large values
        f = new BigFraction(Integer.MAX_VALUE, 1);
        f = f.reciprocal();
        assertEquals(1, f.getNumeratorAsInt());
        assertEquals(Integer.MAX_VALUE, f.getDenominatorAsInt());
    }

    @Test
    void testNegate() {
        BigFraction f = null;

        f = new BigFraction(50, 75);
        f = f.negate();
        assertEquals(-2, f.getNumeratorAsInt());
        assertEquals(3, f.getDenominatorAsInt());

        f = new BigFraction(-50, 75);
        f = f.negate();
        assertEquals(2, f.getNumeratorAsInt());
        assertEquals(3, f.getDenominatorAsInt());

        // large values
        f = new BigFraction(Integer.MAX_VALUE - 1, Integer.MAX_VALUE);
        f = f.negate();
        assertEquals(Integer.MIN_VALUE + 2, f.getNumeratorAsInt());
        assertEquals(Integer.MAX_VALUE, f.getDenominatorAsInt());

    }

    @Test
    void testAdd() {
        BigFraction a = new BigFraction(1, 2);
        BigFraction b = new BigFraction(2, 3);

        customCustomAssertFraction(1, 1, a.add(a));
        customCustomAssertFraction(7, 6, a.add(b));
        customCustomAssertFraction(7, 6, b.add(a));
        customCustomAssertFraction(4, 3, b.add(b));

        BigFraction f1 = new BigFraction(Integer.MAX_VALUE - 1, 1);
        BigFraction f2 = BigFraction.ONE;
        BigFraction f = f1.add(f2);
        assertEquals(Integer.MAX_VALUE, f.getNumeratorAsInt());
        assertEquals(1, f.getDenominatorAsInt());

        f1 = new BigFraction(-1, 13 * 13 * 2 * 2);
        f2 = new BigFraction(-2, 13 * 17 * 2);
        f = f1.add(f2);
        assertEquals(13 * 13 * 17 * 2 * 2, f.getDenominatorAsInt());
        assertEquals(-17 - 2 * 13 * 2, f.getNumeratorAsInt());

        try {
            f.add((BigFraction) null);
            fail("expecting NullArgumentException");
        } catch (NullArgumentException ex) {
        }

        // if this fraction is added naively, it will overflow.
        // check that it doesn't.
        f1 = new BigFraction(1, 32768 * 3);
        f2 = new BigFraction(1, 59049);
        f = f1.add(f2);
        assertEquals(52451, f.getNumeratorAsInt());
        assertEquals(1934917632, f.getDenominatorAsInt());

        f1 = new BigFraction(Integer.MIN_VALUE, 3);
        f2 = new BigFraction(1, 3);
        f = f1.add(f2);
        assertEquals(Integer.MIN_VALUE + 1, f.getNumeratorAsInt());
        assertEquals(3, f.getDenominatorAsInt());

        f1 = new BigFraction(Integer.MAX_VALUE - 1, 1);
        f = f1.add(BigInteger.ONE);
        assertEquals(Integer.MAX_VALUE, f.getNumeratorAsInt());
        assertEquals(1, f.getDenominatorAsInt());

        f = f.add(BigInteger.ZERO);
        assertEquals(Integer.MAX_VALUE, f.getNumeratorAsInt());
        assertEquals(1, f.getDenominatorAsInt());

        f1 = new BigFraction(Integer.MAX_VALUE - 1, 1);
        f = f1.add(1);
        assertEquals(Integer.MAX_VALUE, f.getNumeratorAsInt());
        assertEquals(1, f.getDenominatorAsInt());

        f = f.add(0);
        assertEquals(Integer.MAX_VALUE, f.getNumeratorAsInt());
        assertEquals(1, f.getDenominatorAsInt());

        f1 = new BigFraction(Integer.MAX_VALUE - 1, 1);
        f = f1.add(1l);
        assertEquals(Integer.MAX_VALUE, f.getNumeratorAsInt());
        assertEquals(1, f.getDenominatorAsInt());

        f = f.add(0l);
        assertEquals(Integer.MAX_VALUE, f.getNumeratorAsInt());
        assertEquals(1, f.getDenominatorAsInt());

    }

    @Test
    void testDivide() {
        BigFraction a = new BigFraction(1, 2);
        BigFraction b = new BigFraction(2, 3);

        customCustomAssertFraction(1, 1, a.divide(a));
        customCustomAssertFraction(3, 4, a.divide(b));
        customCustomAssertFraction(4, 3, b.divide(a));
        customCustomAssertFraction(1, 1, b.divide(b));

        BigFraction f1 = new BigFraction(3, 5);
        BigFraction f2 = BigFraction.ZERO;
        try {
            f1.divide(f2);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
        }

        f1 = new BigFraction(0, 5);
        f2 = new BigFraction(2, 7);
        BigFraction f = f1.divide(f2);
        assertSame(BigFraction.ZERO, f);

        f1 = new BigFraction(2, 7);
        f2 = BigFraction.ONE;
        f = f1.divide(f2);
        assertEquals(2, f.getNumeratorAsInt());
        assertEquals(7, f.getDenominatorAsInt());

        f1 = new BigFraction(1, Integer.MAX_VALUE);
        f = f1.divide(f1);
        assertEquals(1, f.getNumeratorAsInt());
        assertEquals(1, f.getDenominatorAsInt());

        f1 = new BigFraction(Integer.MIN_VALUE, Integer.MAX_VALUE);
        f2 = new BigFraction(1, Integer.MAX_VALUE);
        f = f1.divide(f2);
        assertEquals(Integer.MIN_VALUE, f.getNumeratorAsInt());
        assertEquals(1, f.getDenominatorAsInt());

        try {
            f.divide((BigFraction) null);
            fail("expecting NullArgumentException");
        } catch (NullArgumentException ex) {
        }

        f1 = new BigFraction(Integer.MIN_VALUE, Integer.MAX_VALUE);
        f = f1.divide(BigInteger.valueOf(Integer.MIN_VALUE));
        assertEquals(Integer.MAX_VALUE, f.getDenominatorAsInt());
        assertEquals(1, f.getNumeratorAsInt());

        f1 = new BigFraction(Integer.MIN_VALUE, Integer.MAX_VALUE);
        f = f1.divide(Integer.MIN_VALUE);
        assertEquals(Integer.MAX_VALUE, f.getDenominatorAsInt());
        assertEquals(1, f.getNumeratorAsInt());

        f1 = new BigFraction(Integer.MIN_VALUE, Integer.MAX_VALUE);
        f = f1.divide((long) Integer.MIN_VALUE);
        assertEquals(Integer.MAX_VALUE, f.getDenominatorAsInt());
        assertEquals(1, f.getNumeratorAsInt());

    }

    @Test
    void testMultiply() {
        BigFraction a = new BigFraction(1, 2);
        BigFraction b = new BigFraction(2, 3);

        customCustomAssertFraction(1, 4, a.multiply(a));
        customCustomAssertFraction(1, 3, a.multiply(b));
        customCustomAssertFraction(1, 3, b.multiply(a));
        customCustomAssertFraction(4, 9, b.multiply(b));

        BigFraction f1 = new BigFraction(Integer.MAX_VALUE, 1);
        BigFraction f2 = new BigFraction(Integer.MIN_VALUE, Integer.MAX_VALUE);
        BigFraction f = f1.multiply(f2);
        assertEquals(Integer.MIN_VALUE, f.getNumeratorAsInt());
        assertEquals(1, f.getDenominatorAsInt());

        f = f2.multiply(Integer.MAX_VALUE);
        assertEquals(Integer.MIN_VALUE, f.getNumeratorAsInt());
        assertEquals(1, f.getDenominatorAsInt());

        f = f2.multiply((long) Integer.MAX_VALUE);
        assertEquals(Integer.MIN_VALUE, f.getNumeratorAsInt());
        assertEquals(1, f.getDenominatorAsInt());

        try {
            f.multiply((BigFraction) null);
            fail("expecting NullArgumentException");
        } catch (NullArgumentException ex) {
        }

    }

    @Test
    void testSubtract() {
        BigFraction a = new BigFraction(1, 2);
        BigFraction b = new BigFraction(2, 3);

        customCustomAssertFraction(0, 1, a.subtract(a));
        customCustomAssertFraction(-1, 6, a.subtract(b));
        customCustomAssertFraction(1, 6, b.subtract(a));
        customCustomAssertFraction(0, 1, b.subtract(b));

        BigFraction f = new BigFraction(1, 1);
        try {
            f.subtract((BigFraction) null);
            fail("expecting NullArgumentException");
        } catch (NullArgumentException ex) {
        }

        // if this fraction is subtracted naively, it will overflow.
        // check that it doesn't.
        BigFraction f1 = new BigFraction(1, 32768 * 3);
        BigFraction f2 = new BigFraction(1, 59049);
        f = f1.subtract(f2);
        assertEquals(-13085, f.getNumeratorAsInt());
        assertEquals(1934917632, f.getDenominatorAsInt());

        f1 = new BigFraction(Integer.MIN_VALUE, 3);
        f2 = new BigFraction(1, 3).negate();
        f = f1.subtract(f2);
        assertEquals(Integer.MIN_VALUE + 1, f.getNumeratorAsInt());
        assertEquals(3, f.getDenominatorAsInt());

        f1 = new BigFraction(Integer.MAX_VALUE, 1);
        f2 = BigFraction.ONE;
        f = f1.subtract(f2);
        assertEquals(Integer.MAX_VALUE - 1, f.getNumeratorAsInt());
        assertEquals(1, f.getDenominatorAsInt());

    }

    @Test
    void testBigDecimalValue() {
        assertEquals(new BigDecimal(0.5), new BigFraction(1, 2).bigDecimalValue());
        assertEquals(new BigDecimal("0.0003"), new BigFraction(3, 10000).bigDecimalValue());
        assertEquals(new BigDecimal("0"), new BigFraction(1, 3).bigDecimalValue(RoundingMode.DOWN));
        assertEquals(new BigDecimal("0.333"), new BigFraction(1, 3).bigDecimalValue(3, RoundingMode.DOWN));
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testEqualsAndHashCode() {
        BigFraction zero = new BigFraction(0, 1);
        BigFraction nullFraction = null;
        assertEquals(zero, zero);
        assertNotEquals(zero, nullFraction);
        assertNotEquals(zero, Double.valueOf(0));
        BigFraction zero2 = new BigFraction(0, 2);
        assertEquals(zero, zero2);
        assertEquals(zero.hashCode(), zero2.hashCode());
        BigFraction one = new BigFraction(1, 1);
        assertFalse((one.equals(zero) || zero.equals(one)));
        assertEquals(BigFraction.ONE, one);
    }

    @Test
    void testGCD() {
      BigFraction first = new BigFraction(1, 3);
      BigFraction second = new BigFraction(2, 5);
      BigFraction third = new BigFraction(3, 7);
      BigFraction gcd1 = first.gcd(second);
        assertEquals(gcd1, BigFraction.getReducedFraction(1, 15));
      BigFraction gcd2 = gcd1.gcd(third);
        assertEquals(gcd2, BigFraction.getReducedFraction(1, 105));

      // example from https://math.stackexchange.com/a/151089
      BigFraction x = new BigFraction(3, 7);
      BigFraction y = new BigFraction(12, 22);
      BigFraction gcd = x.gcd(y);
        assertEquals(gcd, BigFraction.getReducedFraction(3, 77));

      x = new BigFraction(13, 6);
      y = new BigFraction(3, 4);
      gcd = x.gcd(y);
        assertEquals(gcd, BigFraction.getReducedFraction(1, 12));

    }

    @Test
    void testLCM() {
      BigFraction first = new BigFraction(1, 3);
      BigFraction second = new BigFraction(2, 5);
      BigFraction third = new BigFraction(3, 7);
      BigFraction lcm1 = first.lcm(second);
        assertEquals(lcm1, BigFraction.getReducedFraction(2, 1));
      BigFraction lcm2 = lcm1.lcm(third);
        assertEquals(lcm2, BigFraction.getReducedFraction(6, 1));
    }

    @Test
    void testGetReducedFraction() {
        BigFraction threeFourths = new BigFraction(3, 4);
        assertEquals(threeFourths, BigFraction.getReducedFraction(6, 8));
        assertEquals(BigFraction.ZERO, BigFraction.getReducedFraction(0, -1));
        try {
            BigFraction.getReducedFraction(1, 0);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        assertEquals(-1, BigFraction.getReducedFraction(2, Integer.MIN_VALUE).getNumeratorAsInt());
        assertEquals(-1, BigFraction.getReducedFraction(1, -1).getNumeratorAsInt());
    }

    @Test
    void testPercentage() {
        assertEquals(50.0, new BigFraction(1, 2).percentageValue(), 1.0e-15);
    }

    @Test
    void testPow() {
        assertEquals(new BigFraction(8192, 1594323), new BigFraction(2, 3).pow(13));
        assertEquals(new BigFraction(8192, 1594323), new BigFraction(2, 3).pow(13l));
        assertEquals(new BigFraction(8192, 1594323), new BigFraction(2, 3).pow(BigInteger.valueOf(13l)));
        assertEquals(BigFraction.ONE, new BigFraction(2, 3).pow(0));
        assertEquals(BigFraction.ONE, new BigFraction(2, 3).pow(0l));
        assertEquals(BigFraction.ONE, new BigFraction(2, 3).pow(BigInteger.valueOf(0l)));
        assertEquals(new BigFraction(1594323, 8192), new BigFraction(2, 3).pow(-13));
        assertEquals(new BigFraction(1594323, 8192), new BigFraction(2, 3).pow(-13l));
        assertEquals(new BigFraction(1594323, 8192), new BigFraction(2, 3).pow(BigInteger.valueOf(-13l)));
    }

    @Test
    void testMath340() {
        BigFraction fractionA = new BigFraction(0.00131);
        BigFraction fractionB = new BigFraction(.37).reciprocal();
        BigFraction errorResult = fractionA.multiply(fractionB);
        BigFraction correctResult = new BigFraction(fractionA.getNumerator().multiply(fractionB.getNumerator()),
                                                    fractionA.getDenominator().multiply(fractionB.getDenominator()));
        assertEquals(correctResult, errorResult);
    }

    @Test
    void testNormalizedEquals() {
        assertEquals(new BigFraction(237, -3871), new BigFraction(-51l, 833l));
    }

    @Test
    void testSerial() {
        BigFraction[] fractions = {
            new BigFraction(3, 4), BigFraction.ONE, BigFraction.ZERO,
            new BigFraction(17), new BigFraction(FastMath.PI, 1000),
            new BigFraction(-5, 2)
        };
        for (BigFraction fraction : fractions) {
            assertEquals(fraction, UnitTestUtils.serializeAndRecover(fraction));
        }
    }

    @Test
    void testConvergents() {
        // OEIS A002485, Numerators of convergents to Pi (https://oeis.org/A002485)
        // 0, 1, 3, 22, 333, 355, 103993, 104348, 208341, 312689, 833719, 1146408, 4272943, 5419351, 80143857, 165707065, 245850922
        // OEIS A002486, Apart from two leading terms (which are present by convention), denominators of convergents to Pi (https://oeis.org/A002486)
        // 1, 0, 1,  7, 106, 113,  33102,  33215,  66317,  99532, 265381,  364913, 1360120, 1725033, 25510582,  52746197, 78256779
        List<BigFraction> convergents = BigFraction.convergents(FastMath.PI, 20).collect(Collectors.toList());
        assertEquals(new BigFraction(       3,        1), convergents.get( 0));
        assertEquals(new BigFraction(      22,        7), convergents.get( 1));
        assertEquals(new BigFraction(     333,      106), convergents.get( 2));
        assertEquals(new BigFraction(     355,      113), convergents.get( 3));
        assertEquals(new BigFraction(  103993,    33102), convergents.get( 4));
        assertEquals(new BigFraction(  104348,    33215), convergents.get( 5));
        assertEquals(new BigFraction(  208341,    66317), convergents.get( 6));
        assertEquals(new BigFraction(  312689,    99532), convergents.get( 7));
        assertEquals(new BigFraction(  833719,   265381), convergents.get( 8));
        assertEquals(new BigFraction( 1146408,   364913), convergents.get( 9));
        assertEquals(new BigFraction( 4272943,  1360120), convergents.get(10));
        assertEquals(new BigFraction( 5419351,  1725033), convergents.get(11));
        assertEquals(new BigFraction(80143857, 25510582), convergents.get(12));
        assertEquals(13, convergents.size());
    }

    @Test
    void testLimitedConvergents() {
        double value = FastMath.PI;
        assertEquals(new BigFraction(  208341,    66317),
                BigFraction.convergent(value, 7, (p, q) -> Precision.equals(p / (double) q, value, 1)).getKey());
    }

    @Test
    void testTruncatedConvergents() {
        final double value = FastMath.PI;
        assertEquals(new BigFraction(   355,   113),
                BigFraction.convergent(value, 20, (p, q) -> FastMath.abs(p / (double) q - value) < 1.0e-6).getKey());
        assertEquals(new BigFraction(312689, 99532),
                BigFraction.convergent(value, 20, (p, q) -> FastMath.abs(p / (double) q - value) < 1.0e-10).getKey());
    }

    @Test
    void testOutOfRange() {
        BigFraction f = new BigFraction(new BigInteger("1175443811202636889584648110261699215671929253339678037082566183036829784156100300341131818417591797406644569806405529752410539491566888996766640542430075310377605462098357361563685103574645710283612852841417362211504458393792053529953230572830415970785545248189857341548686469982966457542855773477057255734051"),
                                        new BigInteger("32626522339992622633551470546282737778505821290344832738793182277348616222987431136114480634269341408071340993046760559082031250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"));
        // reference value compured using Emacs-calc with 20 digits
        assertEquals(36.027247984128935385, f.doubleValue(), 1.0e-15);
    }

}
