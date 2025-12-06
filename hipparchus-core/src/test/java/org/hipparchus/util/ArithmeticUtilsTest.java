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
package org.hipparchus.util;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.random.RandomDataGenerator;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the {@link ArithmeticUtils} class.
 */
class ArithmeticUtilsTest {

    @Test
    void testAddAndCheck() {
        int big = Integer.MAX_VALUE;
        int bigNeg = Integer.MIN_VALUE;
        assertEquals(big, ArithmeticUtils.addAndCheck(big, 0));
        try {
            ArithmeticUtils.addAndCheck(big, 1);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
        }
        try {
            ArithmeticUtils.addAndCheck(bigNeg, -1);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
        }
    }

    @Test
    void testAddAndCheckLong() {
        long max = Long.MAX_VALUE;
        long min = Long.MIN_VALUE;
        assertEquals(max, ArithmeticUtils.addAndCheck(max, 0L));
        assertEquals(min, ArithmeticUtils.addAndCheck(min, 0L));
        assertEquals(max, ArithmeticUtils.addAndCheck(0L, max));
        assertEquals(min, ArithmeticUtils.addAndCheck(0L, min));
        assertEquals(1, ArithmeticUtils.addAndCheck(-1L, 2L));
        assertEquals(1, ArithmeticUtils.addAndCheck(2L, -1L));
        assertEquals(-3, ArithmeticUtils.addAndCheck(-2L, -1L));
        assertEquals(min, ArithmeticUtils.addAndCheck(min + 1, -1L));
        assertEquals(-1, ArithmeticUtils.addAndCheck(min, max));
        testAddAndCheckLongFailure(max, 1L);
        testAddAndCheckLongFailure(min, -1L);
        testAddAndCheckLongFailure(1L, max);
        testAddAndCheckLongFailure(-1L, min);
        testAddAndCheckLongFailure(max, max);
        testAddAndCheckLongFailure(min, min);
    }

    @Test
    void testGcd() {
        int a = 30;
        int b = 50;
        int c = 77;

        assertEquals(0, ArithmeticUtils.gcd(0, 0));

        assertEquals(b, ArithmeticUtils.gcd(0, b));
        assertEquals(a, ArithmeticUtils.gcd(a, 0));
        assertEquals(b, ArithmeticUtils.gcd(0, -b));
        assertEquals(a, ArithmeticUtils.gcd(-a, 0));

        assertEquals(10, ArithmeticUtils.gcd(a, b));
        assertEquals(10, ArithmeticUtils.gcd(-a, b));
        assertEquals(10, ArithmeticUtils.gcd(a, -b));
        assertEquals(10, ArithmeticUtils.gcd(-a, -b));

        assertEquals(1, ArithmeticUtils.gcd(a, c));
        assertEquals(1, ArithmeticUtils.gcd(-a, c));
        assertEquals(1, ArithmeticUtils.gcd(a, -c));
        assertEquals(1, ArithmeticUtils.gcd(-a, -c));

        assertEquals(3 * (1<<15), ArithmeticUtils.gcd(3 * (1<<20), 9 * (1<<15)));

        assertEquals(Integer.MAX_VALUE, ArithmeticUtils.gcd(Integer.MAX_VALUE, 0));
        assertEquals(Integer.MAX_VALUE, ArithmeticUtils.gcd(-Integer.MAX_VALUE, 0));
        assertEquals(1<<30, ArithmeticUtils.gcd(1<<30, -Integer.MIN_VALUE));
        try {
            // gcd(Integer.MIN_VALUE, 0) > Integer.MAX_VALUE
            ArithmeticUtils.gcd(Integer.MIN_VALUE, 0);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException expected) {
            // expected
        }
        try {
            // gcd(0, Integer.MIN_VALUE) > Integer.MAX_VALUE
            ArithmeticUtils.gcd(0, Integer.MIN_VALUE);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException expected) {
            // expected
        }
        try {
            // gcd(Integer.MIN_VALUE, Integer.MIN_VALUE) > Integer.MAX_VALUE
            ArithmeticUtils.gcd(Integer.MIN_VALUE, Integer.MIN_VALUE);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException expected) {
            // expected
        }
    }

    @Test
    void testGcdConsistency() {
        int[] primeList = {19, 23, 53, 67, 73, 79, 101, 103, 111, 131};
        ArrayList<Integer> primes = new ArrayList<Integer>();
        for (int i = 0; i < primeList.length; i++) {
            primes.add(Integer.valueOf(primeList[i]));
        }
        RandomDataGenerator randomData = new RandomDataGenerator();
        for (int i = 0; i < 20; i++) {
            Object[] sample = randomData.nextSample(primes, 4);
            int p1 = ((Integer) sample[0]).intValue();
            int p2 = ((Integer) sample[1]).intValue();
            int p3 = ((Integer) sample[2]).intValue();
            int p4 = ((Integer) sample[3]).intValue();
            int i1 = p1 * p2 * p3;
            int i2 = p1 * p2 * p4;
            int gcd = p1 * p2;
            assertEquals(gcd, ArithmeticUtils.gcd(i1, i2));
            long l1 = i1;
            long l2 = i2;
            assertEquals(gcd, ArithmeticUtils.gcd(l1, l2));
        }
    }

    @Test
    void testGcdLong(){
        long a = 30;
        long b = 50;
        long c = 77;

        assertEquals(0, ArithmeticUtils.gcd(0L, 0));

        assertEquals(b, ArithmeticUtils.gcd(0, b));
        assertEquals(a, ArithmeticUtils.gcd(a, 0));
        assertEquals(b, ArithmeticUtils.gcd(0, -b));
        assertEquals(a, ArithmeticUtils.gcd(-a, 0));

        assertEquals(10, ArithmeticUtils.gcd(a, b));
        assertEquals(10, ArithmeticUtils.gcd(-a, b));
        assertEquals(10, ArithmeticUtils.gcd(a, -b));
        assertEquals(10, ArithmeticUtils.gcd(-a, -b));

        assertEquals(1, ArithmeticUtils.gcd(a, c));
        assertEquals(1, ArithmeticUtils.gcd(-a, c));
        assertEquals(1, ArithmeticUtils.gcd(a, -c));
        assertEquals(1, ArithmeticUtils.gcd(-a, -c));

        assertEquals(3L * (1L<<45), ArithmeticUtils.gcd(3L * (1L<<50), 9L * (1L<<45)));

        assertEquals(1L<<45, ArithmeticUtils.gcd(1L<<45, Long.MIN_VALUE));

        assertEquals(Long.MAX_VALUE, ArithmeticUtils.gcd(Long.MAX_VALUE, 0L));
        assertEquals(Long.MAX_VALUE, ArithmeticUtils.gcd(-Long.MAX_VALUE, 0L));
        assertEquals(1, ArithmeticUtils.gcd(60247241209L, 153092023L));
        try {
            // gcd(Long.MIN_VALUE, 0) > Long.MAX_VALUE
            ArithmeticUtils.gcd(Long.MIN_VALUE, 0);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException expected) {
            // expected
        }
        try {
            // gcd(0, Long.MIN_VALUE) > Long.MAX_VALUE
            ArithmeticUtils.gcd(0, Long.MIN_VALUE);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException expected) {
            // expected
        }
        try {
            // gcd(Long.MIN_VALUE, Long.MIN_VALUE) > Long.MAX_VALUE
            ArithmeticUtils.gcd(Long.MIN_VALUE, Long.MIN_VALUE);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException expected) {
            // expected
        }
    }


    @Test
    void testLcm() {
        int a = 30;
        int b = 50;
        int c = 77;

        assertEquals(0, ArithmeticUtils.lcm(0, b));
        assertEquals(0, ArithmeticUtils.lcm(a, 0));
        assertEquals(b, ArithmeticUtils.lcm(1, b));
        assertEquals(a, ArithmeticUtils.lcm(a, 1));
        assertEquals(150, ArithmeticUtils.lcm(a, b));
        assertEquals(150, ArithmeticUtils.lcm(-a, b));
        assertEquals(150, ArithmeticUtils.lcm(a, -b));
        assertEquals(150, ArithmeticUtils.lcm(-a, -b));
        assertEquals(2310, ArithmeticUtils.lcm(a, c));

        // Assert that no intermediate value overflows:
        // The naive implementation of lcm(a,b) would be (a*b)/gcd(a,b)
        assertEquals((1<<20)*15, ArithmeticUtils.lcm((1<<20)*3, (1<<20)*5));

        // Special case
        assertEquals(0, ArithmeticUtils.lcm(0, 0));

        try {
            // lcm == abs(MIN_VALUE) cannot be represented as a nonnegative int
            ArithmeticUtils.lcm(Integer.MIN_VALUE, 1);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException expected) {
            // expected
        }

        try {
            // lcm == abs(MIN_VALUE) cannot be represented as a nonnegative int
            ArithmeticUtils.lcm(Integer.MIN_VALUE, 1<<20);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException expected) {
            // expected
        }

        try {
            ArithmeticUtils.lcm(Integer.MAX_VALUE, Integer.MAX_VALUE - 1);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException expected) {
            // expected
        }
    }

    @Test
    void testLcmLong() {
        long a = 30;
        long b = 50;
        long c = 77;

        assertEquals(0, ArithmeticUtils.lcm(0, b));
        assertEquals(0, ArithmeticUtils.lcm(a, 0));
        assertEquals(b, ArithmeticUtils.lcm(1, b));
        assertEquals(a, ArithmeticUtils.lcm(a, 1));
        assertEquals(150, ArithmeticUtils.lcm(a, b));
        assertEquals(150, ArithmeticUtils.lcm(-a, b));
        assertEquals(150, ArithmeticUtils.lcm(a, -b));
        assertEquals(150, ArithmeticUtils.lcm(-a, -b));
        assertEquals(2310, ArithmeticUtils.lcm(a, c));

        assertEquals(Long.MAX_VALUE, ArithmeticUtils.lcm(60247241209L, 153092023L));

        // Assert that no intermediate value overflows:
        // The naive implementation of lcm(a,b) would be (a*b)/gcd(a,b)
        assertEquals((1L<<50)*15, ArithmeticUtils.lcm((1L<<45)*3, (1L<<50)*5));

        // Special case
        assertEquals(0L, ArithmeticUtils.lcm(0L, 0L));

        try {
            // lcm == abs(MIN_VALUE) cannot be represented as a nonnegative int
            ArithmeticUtils.lcm(Long.MIN_VALUE, 1);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException expected) {
            // expected
        }

        try {
            // lcm == abs(MIN_VALUE) cannot be represented as a nonnegative int
            ArithmeticUtils.lcm(Long.MIN_VALUE, 1<<20);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException expected) {
            // expected
        }

        assertEquals((long) Integer.MAX_VALUE * (Integer.MAX_VALUE - 1),
            ArithmeticUtils.lcm((long)Integer.MAX_VALUE, Integer.MAX_VALUE - 1));
        try {
            ArithmeticUtils.lcm(Long.MAX_VALUE, Long.MAX_VALUE - 1);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException expected) {
            // expected
        }
    }

    @Test
    void testMulAndCheck() {
        int big = Integer.MAX_VALUE;
        int bigNeg = Integer.MIN_VALUE;
        assertEquals(big, ArithmeticUtils.mulAndCheck(big, 1));
        try {
            ArithmeticUtils.mulAndCheck(big, 2);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
        }
        try {
            ArithmeticUtils.mulAndCheck(bigNeg, 2);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
        }
    }

    @Test
    void testMulAndCheckLong() {
        long max = Long.MAX_VALUE;
        long min = Long.MIN_VALUE;
        assertEquals(max, ArithmeticUtils.mulAndCheck(max, 1L));
        assertEquals(min, ArithmeticUtils.mulAndCheck(min, 1L));
        assertEquals(0L, ArithmeticUtils.mulAndCheck(max, 0L));
        assertEquals(0L, ArithmeticUtils.mulAndCheck(min, 0L));
        assertEquals(max, ArithmeticUtils.mulAndCheck(1L, max));
        assertEquals(min, ArithmeticUtils.mulAndCheck(1L, min));
        assertEquals(0L, ArithmeticUtils.mulAndCheck(0L, max));
        assertEquals(0L, ArithmeticUtils.mulAndCheck(0L, min));
        assertEquals(1L, ArithmeticUtils.mulAndCheck(-1L, -1L));
        assertEquals(min, ArithmeticUtils.mulAndCheck(min / 2, 2));
        testMulAndCheckLongFailure(max, 2L);
        testMulAndCheckLongFailure(2L, max);
        testMulAndCheckLongFailure(min, 2L);
        testMulAndCheckLongFailure(2L, min);
        testMulAndCheckLongFailure(min, -1L);
        testMulAndCheckLongFailure(-1L, min);
    }

    @Test
    void testSubAndCheck() {
        int big = Integer.MAX_VALUE;
        int bigNeg = Integer.MIN_VALUE;
        assertEquals(big, ArithmeticUtils.subAndCheck(big, 0));
        assertEquals(bigNeg + 1, ArithmeticUtils.subAndCheck(bigNeg, -1));
        assertEquals(-1, ArithmeticUtils.subAndCheck(bigNeg, -big));
        try {
            ArithmeticUtils.subAndCheck(big, -1);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
        }
        try {
            ArithmeticUtils.subAndCheck(bigNeg, 1);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
        }
    }

    @Test
    void testSubAndCheckErrorMessage() {
        int big = Integer.MAX_VALUE;
        try {
            ArithmeticUtils.subAndCheck(big, -1);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
            assertTrue(ex.getMessage().length() > 1);
        }
    }

    @Test
    void testSubAndCheckLong() {
        long max = Long.MAX_VALUE;
        long min = Long.MIN_VALUE;
        assertEquals(max, ArithmeticUtils.subAndCheck(max, 0));
        assertEquals(min, ArithmeticUtils.subAndCheck(min, 0));
        assertEquals(-max, ArithmeticUtils.subAndCheck(0, max));
        assertEquals(min + 1, ArithmeticUtils.subAndCheck(min, -1));
        // min == -1-max
        assertEquals(-1, ArithmeticUtils.subAndCheck(-max - 1, -max));
        assertEquals(max, ArithmeticUtils.subAndCheck(-1, -1 - max));
        testSubAndCheckLongFailure(0L, min);
        testSubAndCheckLongFailure(max, -1L);
        testSubAndCheckLongFailure(min, 1L);
    }

    @Test
    void testPow() {

        assertEquals(1801088541, ArithmeticUtils.pow(21, 7));
        assertEquals(1, ArithmeticUtils.pow(21, 0));
        try {
            ArithmeticUtils.pow(21, -7);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected behavior
        }

        assertEquals(1801088541, ArithmeticUtils.pow(21, 7));
        assertEquals(1, ArithmeticUtils.pow(21, 0));
        try {
            ArithmeticUtils.pow(21, -7);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected behavior
        }

        assertEquals(1801088541L, ArithmeticUtils.pow(21L, 7));
        assertEquals(1L, ArithmeticUtils.pow(21L, 0));
        try {
            ArithmeticUtils.pow(21L, -7);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected behavior
        }

        BigInteger twentyOne = BigInteger.valueOf(21L);
        assertEquals(BigInteger.valueOf(1801088541L), ArithmeticUtils.pow(twentyOne, 7));
        assertEquals(BigInteger.ONE, ArithmeticUtils.pow(twentyOne, 0));
        try {
            ArithmeticUtils.pow(twentyOne, -7);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected behavior
        }

        assertEquals(BigInteger.valueOf(1801088541L), ArithmeticUtils.pow(twentyOne,
                                                                                 7L));
        assertEquals(BigInteger.ONE, ArithmeticUtils.pow(twentyOne, 0L));
        try {
            ArithmeticUtils.pow(twentyOne, -7L);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected behavior
        }

        assertEquals(BigInteger.valueOf(1801088541L), ArithmeticUtils.pow(twentyOne, BigInteger.valueOf(
                        7L)));
        assertEquals(BigInteger.ONE, ArithmeticUtils.pow(twentyOne, BigInteger.ZERO));
        try {
            ArithmeticUtils.pow(twentyOne, BigInteger.valueOf(-7L));
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected behavior
        }

        BigInteger bigOne =
            new BigInteger("1543786922199448028351389769265814882661837148" +
                           "4763915343722775611762713982220306372888519211" +
                           "560905579993523402015636025177602059044911261");
        assertEquals(bigOne, ArithmeticUtils.pow(twentyOne, 103));
        assertEquals(bigOne, ArithmeticUtils.pow(twentyOne, 103L));
        assertEquals(bigOne, ArithmeticUtils.pow(twentyOne, BigInteger.valueOf(
                        103L)));

    }

    @Test
    void testPowIntOverflow() {
        assertThrows(MathRuntimeException.class, () -> ArithmeticUtils.pow(21, 8));
    }

    @Test
    void testPowInt() {
        final int base = 21;

        assertEquals(85766121L,
                            ArithmeticUtils.pow(base, 6));
        assertEquals(1801088541L,
                            ArithmeticUtils.pow(base, 7));
    }

    @Test
    void testPowNegativeIntOverflow() {
        assertThrows(MathRuntimeException.class, () -> ArithmeticUtils.pow(-21, 8));
    }

    @Test
    void testPowNegativeInt() {
        final int base = -21;

        assertEquals(85766121,
                            ArithmeticUtils.pow(base, 6));
        assertEquals(-1801088541,
                            ArithmeticUtils.pow(base, 7));
    }

    @Test
    void testPowMinusOneInt() {
        final int base = -1;
        for (int i = 0; i < 100; i++) {
            final int pow = ArithmeticUtils.pow(base, i);
            assertEquals(i % 2 == 0 ? 1 : -1, pow, "i: " + i);
        }
    }

    @Test
    void testPowOneInt() {
        final int base = 1;
        for (int i = 0; i < 100; i++) {
            final int pow = ArithmeticUtils.pow(base, i);
            assertEquals(1, pow, "i: " + i);
        }
    }

    @Test
    void testPowLongOverflow() {
        assertThrows(MathRuntimeException.class, () -> ArithmeticUtils.pow(21, 15));
    }

    @Test
    void testPowLong() {
        final long base = 21;

        assertEquals(154472377739119461L,
                            ArithmeticUtils.pow(base, 13));
        assertEquals(3243919932521508681L,
                            ArithmeticUtils.pow(base, 14));
    }

    @Test
    void testPowNegativeLongOverflow() {
        assertThrows(MathRuntimeException.class, () -> ArithmeticUtils.pow(-21L, 15));
    }

    @Test
    void testPowNegativeLong() {
        final long base = -21;

        assertEquals(-154472377739119461L,
                            ArithmeticUtils.pow(base, 13));
        assertEquals(3243919932521508681L,
                            ArithmeticUtils.pow(base, 14));
    }

    @Test
    void testPowMinusOneLong() {
        final long base = -1;
        for (int i = 0; i < 100; i++) {
            final long pow = ArithmeticUtils.pow(base, i);
            assertEquals(i % 2 == 0 ? 1 : -1, pow, "i: " + i);
        }
    }

    @Test
    void testPowOneLong() {
        final long base = 1;
        for (int i = 0; i < 100; i++) {
            final long pow = ArithmeticUtils.pow(base, i);
            assertEquals(1, pow, "i: " + i);
        }
    }

    @Test
    void testIsPowerOfTwo() {
        final int n = 1025;
        final boolean[] expected = new boolean[n];
        Arrays.fill(expected, false);
        for (int i = 1; i < expected.length; i *= 2) {
            expected[i] = true;
        }
        for (int i = 0; i < expected.length; i++) {
            final boolean actual = ArithmeticUtils.isPowerOfTwo(i);
            assertEquals(actual, expected[i], Integer.toString(i));
        }
    }

    private void testAddAndCheckLongFailure(long a, long b) {
        try {
            ArithmeticUtils.addAndCheck(a, b);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
            // success
        }
    }

    private void testMulAndCheckLongFailure(long a, long b) {
        try {
            ArithmeticUtils.mulAndCheck(a, b);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
            // success
        }
    }

    private void testSubAndCheckLongFailure(long a, long b) {
        try {
            ArithmeticUtils.subAndCheck(a, b);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
            // success
        }
    }

    /**
     * Testing helper method.
     * @return an array of int numbers containing corner cases:<ul>
     * <li>values near the beginning of int range,</li>
     * <li>values near the end of int range,</li>
     * <li>values near zero</li>
     * <li>and some randomly distributed values.</li>
     * </ul>
     */
    private static int[] getIntSpecialCases() {
        int[] ints = new int[100];
        int i = 0;
        ints[i++] = Integer.MAX_VALUE;
        ints[i++] = Integer.MAX_VALUE - 1;
        ints[i++] = 100;
        ints[i++] = 101;
        ints[i++] = 102;
        ints[i++] = 300;
        ints[i++] = 567;
        for (int j = 0; j < 20; j++) {
            ints[i++] = j;
        }
        for (int j = i - 1; j >= 0; j--) {
            ints[i++] = ints[j] > 0 ? -ints[j] : Integer.MIN_VALUE;
        }
        java.util.Random r = new java.util.Random(System.nanoTime());
        while (i < ints.length) {
            ints[i++] = r.nextInt();
        }
        return ints;
    }

    /**
     * Testing helper method.
     * @return an array of long numbers containing corner cases:<ul>
     * <li>values near the beginning of long range,</li>
     * <li>values near the end of long range,</li>
     * <li>values near the beginning of int range,</li>
     * <li>values near the end of int range,</li>
     * <li>values near zero</li>
     * <li>and some randomly distributed values.</li>
     * </ul>
     */
    private static long[] getLongSpecialCases() {
        long[] longs = new long[100];
        int i = 0;
        longs[i++] = Long.MAX_VALUE;
        longs[i++] = Long.MAX_VALUE - 1L;
        longs[i++] = (long) Integer.MAX_VALUE + 1L;
        longs[i++] = Integer.MAX_VALUE;
        longs[i++] = Integer.MAX_VALUE - 1;
        longs[i++] = 100L;
        longs[i++] = 101L;
        longs[i++] = 102L;
        longs[i++] = 300L;
        longs[i++] = 567L;
        for (int j = 0; j < 20; j++) {
            longs[i++] = j;
        }
        for (int j = i - 1; j >= 0; j--) {
            longs[i++] = longs[j] > 0L ? -longs[j] : Long.MIN_VALUE;
        }
        java.util.Random r = new java.util.Random(System.nanoTime());
        while (i < longs.length) {
            longs[i++] = r.nextLong();
        }
        return longs;
    }

    private static long toUnsignedLong(int number) {
        return number < 0 ? 0x100000000L + (long)number : (long)number;
    }

    private static int remainderUnsignedExpected(int dividend, int divisor) {
        return (int)remainderUnsignedExpected(toUnsignedLong(dividend), toUnsignedLong(divisor));
    }

    private static int divideUnsignedExpected(int dividend, int divisor) {
        return (int)divideUnsignedExpected(toUnsignedLong(dividend), toUnsignedLong(divisor));
    }

    private static BigInteger toUnsignedBigInteger(long number) {
        return number < 0L ? BigInteger.ONE.shiftLeft(64).add(BigInteger.valueOf(number)) : BigInteger.valueOf(number);
    }

    private static long remainderUnsignedExpected(long dividend, long divisor) {
        return toUnsignedBigInteger(dividend).remainder(toUnsignedBigInteger(divisor)).longValue();
    }

    private static long divideUnsignedExpected(long dividend, long divisor) {
        return toUnsignedBigInteger(dividend).divide(toUnsignedBigInteger(divisor)).longValue();
    }

    @Test
    @Timeout(value = 5000L, unit = TimeUnit.MILLISECONDS)
    void testRemainderUnsignedInt() {
        assertEquals(36, ArithmeticUtils.remainderUnsigned(-2147479015, 63));
        assertEquals(6, ArithmeticUtils.remainderUnsigned(-2147479015, 25));
    }

    @Test
    @Timeout(value = 5000L, unit = TimeUnit.MILLISECONDS)
    void testRemainderUnsignedIntSpecialCases() {
        int[] ints = getIntSpecialCases();
        for (int dividend : ints) {
            for (int divisor : ints) {
                if (divisor == 0) {
                    try {
                        ArithmeticUtils.remainderUnsigned(dividend, divisor);
                        fail("Should have failed with ArithmeticException: division by zero");
                    } catch (ArithmeticException e) {
                        // Success.
                    }
                } else {
                    assertEquals(remainderUnsignedExpected(dividend, divisor), ArithmeticUtils.remainderUnsigned(dividend, divisor));
                }
            }
        }
    }

    @Test
    @Timeout(value = 5000L, unit = TimeUnit.MILLISECONDS)
    void testRemainderUnsignedLong() {
        assertEquals(48L, ArithmeticUtils.remainderUnsigned(-2147479015L, 63L));
    }

    //(timeout=5000L)
    @Test
    void testRemainderUnsignedLongSpecialCases() {
        long[] longs = getLongSpecialCases();
        for (long dividend : longs) {
            for (long divisor : longs) {
                if (divisor == 0L) {
                    try {
                        ArithmeticUtils.remainderUnsigned(dividend, divisor);
                        fail("Should have failed with ArithmeticException: division by zero");
                    } catch (ArithmeticException e) {
                        // Success.
                    }
                } else {
                    assertEquals(remainderUnsignedExpected(dividend, divisor), ArithmeticUtils.remainderUnsigned(dividend, divisor));
                }
            }
        }
    }

    @Test
    @Timeout(value = 5000L, unit = TimeUnit.MILLISECONDS)
    void testDivideUnsignedInt() {
        assertEquals(34087115, ArithmeticUtils.divideUnsigned(-2147479015, 63));
        assertEquals(85899531, ArithmeticUtils.divideUnsigned(-2147479015, 25));
        assertEquals(2147483646, ArithmeticUtils.divideUnsigned(-3, 2));
        assertEquals(330382098, ArithmeticUtils.divideUnsigned(-16, 13));
        assertEquals(306783377, ArithmeticUtils.divideUnsigned(-16, 14));
        assertEquals(2, ArithmeticUtils.divideUnsigned(-1, 2147483647));
        assertEquals(2, ArithmeticUtils.divideUnsigned(-2, 2147483647));
        assertEquals(1, ArithmeticUtils.divideUnsigned(-3, 2147483647));
        assertEquals(1, ArithmeticUtils.divideUnsigned(-16, 2147483647));
        assertEquals(1, ArithmeticUtils.divideUnsigned(-16, 2147483646));
    }

    @Test
    @Timeout(value = 5000L, unit = TimeUnit.MILLISECONDS)
    void testDivideUnsignedIntSpecialCases() {
        int[] ints = getIntSpecialCases();
        for (int dividend : ints) {
            for (int divisor : ints) {
                if (divisor == 0) {
                    try {
                        ArithmeticUtils.divideUnsigned(dividend, divisor);
                        fail("Should have failed with ArithmeticException: division by zero");
                    } catch (ArithmeticException e) {
                        // Success.
                    }
                } else {
                    assertEquals(divideUnsignedExpected(dividend, divisor), ArithmeticUtils.divideUnsigned(dividend, divisor));
                }
            }
        }
    }

    @Test
    @Timeout(value = 5000L, unit = TimeUnit.MILLISECONDS)
    void testDivideUnsignedLong() {
        assertEquals(292805461453366231L, ArithmeticUtils.divideUnsigned(-2147479015L, 63L));
    }

    @Test
    @Timeout(value = 5000L, unit = TimeUnit.MILLISECONDS)
    void testDivideUnsignedLongSpecialCases() {
        long[] longs = getLongSpecialCases();
        for (long dividend : longs) {
            for (long divisor : longs) {
                if (divisor == 0L) {
                    try {
                        ArithmeticUtils.divideUnsigned(dividend, divisor);
                        fail("Should have failed with ArithmeticException: division by zero");
                    } catch (ArithmeticException e) {
                        // Success.
                    }
                } else {
                    assertEquals(divideUnsignedExpected(dividend, divisor), ArithmeticUtils.divideUnsigned(dividend, divisor));
                }
            }
        }
    }
}
