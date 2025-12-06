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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the {@link CombinatoricsUtils} class.
 */
class CombinatoricsUtilsTest {

    /** cached binomial coefficients */
    private static final List<Map<Integer, Long>> binomialCache = new ArrayList<Map<Integer, Long>>();

    /** Verify that b(0,0) = 1 */
    @Test
    void test0Choose0() {
        assertEquals(1d, CombinatoricsUtils.binomialCoefficientDouble(0, 0), 0);
        assertEquals(0d, CombinatoricsUtils.binomialCoefficientLog(0, 0), 0);
        assertEquals(1, CombinatoricsUtils.binomialCoefficient(0, 0));
    }

    @Test
    void testBinomialCoefficient() {
        long[] bcoef5 = {
            1,
            5,
            10,
            10,
            5,
            1 };
        long[] bcoef6 = {
            1,
            6,
            15,
            20,
            15,
            6,
            1 };
        for (int i = 0; i < 6; i++) {
            assertEquals(bcoef5[i], CombinatoricsUtils.binomialCoefficient(5, i), "5 choose " + i);
        }
        for (int i = 0; i < 7; i++) {
            assertEquals(bcoef6[i], CombinatoricsUtils.binomialCoefficient(6, i), "6 choose " + i);
        }

        for (int n = 1; n < 10; n++) {
            for (int k = 0; k <= n; k++) {
                assertEquals(binomialCoefficient(n, k), CombinatoricsUtils.binomialCoefficient(n, k), n + " choose " + k);
                assertEquals(binomialCoefficient(n, k), CombinatoricsUtils.binomialCoefficientDouble(n, k), Double.MIN_VALUE, n + " choose " + k);
                assertEquals(FastMath.log(binomialCoefficient(n, k)), CombinatoricsUtils.binomialCoefficientLog(n, k), 10E-12, n + " choose " + k);
            }
        }

        int[] n = { 34, 66, 100, 1500, 1500 };
        int[] k = { 17, 33, 10, 1500 - 4, 4 };
        for (int i = 0; i < n.length; i++) {
            long expected = binomialCoefficient(n[i], k[i]);
            assertEquals(expected,
                CombinatoricsUtils.binomialCoefficient(n[i], k[i]),
                n[i] + " choose " + k[i]);
            assertEquals(expected,
                CombinatoricsUtils.binomialCoefficientDouble(n[i], k[i]), 0.0, n[i] + " choose " + k[i]);
            assertEquals(FastMath.log(expected),
                CombinatoricsUtils.binomialCoefficientLog(n[i], k[i]), 0.0, "log(" + n[i] + " choose " + k[i] + ")");
        }
    }

    @Test
    void testBinomialCoefficientFail() {
        try {
            CombinatoricsUtils.binomialCoefficient(4, 5);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }

        try {
            CombinatoricsUtils.binomialCoefficientDouble(4, 5);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }

        try {
            CombinatoricsUtils.binomialCoefficientLog(4, 5);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }

        try {
            CombinatoricsUtils.binomialCoefficient(-1, -2);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            CombinatoricsUtils.binomialCoefficientDouble(-1, -2);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            CombinatoricsUtils.binomialCoefficientLog(-1, -2);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }

        try {
            CombinatoricsUtils.binomialCoefficient(67, 30);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
            // ignored
        }
        try {
            CombinatoricsUtils.binomialCoefficient(67, 34);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
            // ignored
        }
        double x = CombinatoricsUtils.binomialCoefficientDouble(1030, 515);
        assertTrue(Double
            .isInfinite(x), "expecting infinite binomial coefficient");
    }

    /**
     * Tests correctness for large n and sharpness of upper bound in API doc
     * JIRA: MATH-241
     */
    @Test
    void testBinomialCoefficientLarge() throws Exception {
        // This tests all legal and illegal values for n <= 200.
        for (int n = 0; n <= 200; n++) {
            for (int k = 0; k <= n; k++) {
                long ourResult = -1;
                long exactResult = -1;
                boolean shouldThrow = false;
                boolean didThrow = false;
                try {
                    ourResult = CombinatoricsUtils.binomialCoefficient(n, k);
                } catch (MathRuntimeException ex) {
                    didThrow = true;
                }
                try {
                    exactResult = binomialCoefficient(n, k);
                } catch (MathRuntimeException ex) {
                    shouldThrow = true;
                }
                assertEquals(exactResult, ourResult, n + " choose " + k);
                assertEquals(shouldThrow, didThrow, n + " choose " + k);
                assertTrue((n > 66 || !didThrow), n + " choose " + k);

                if (!shouldThrow && exactResult > 1) {
                    assertEquals(1.,
                        CombinatoricsUtils.binomialCoefficientDouble(n, k) / exactResult, 1e-10, n + " choose " + k);
                    assertEquals(1,
                        CombinatoricsUtils.binomialCoefficientLog(n, k) / FastMath.log(exactResult), 1e-10, n + " choose " + k);
                }
            }
        }

        long ourResult = CombinatoricsUtils.binomialCoefficient(300, 3);
        long exactResult = binomialCoefficient(300, 3);
        assertEquals(exactResult, ourResult);

        ourResult = CombinatoricsUtils.binomialCoefficient(700, 697);
        exactResult = binomialCoefficient(700, 697);
        assertEquals(exactResult, ourResult);

        // This one should throw
        try {
            CombinatoricsUtils.binomialCoefficient(700, 300);
            fail("Expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
            // Expected
        }

        int n = 10000;
        ourResult = CombinatoricsUtils.binomialCoefficient(n, 3);
        exactResult = binomialCoefficient(n, 3);
        assertEquals(exactResult, ourResult);
        assertEquals(1, CombinatoricsUtils.binomialCoefficientDouble(n, 3) / exactResult, 1e-10);
        assertEquals(1, CombinatoricsUtils.binomialCoefficientLog(n, 3) / FastMath.log(exactResult), 1e-10);

    }

    @Test
    void testFactorial() {
        for (int i = 1; i < 21; i++) {
            assertEquals(factorial(i), CombinatoricsUtils.factorial(i), i + "! ");
            assertEquals(factorial(i), CombinatoricsUtils.factorialDouble(i), Double.MIN_VALUE, i + "! ");
            assertEquals(FastMath.log(factorial(i)), CombinatoricsUtils.factorialLog(i), 10E-12, i + "! ");
        }

        assertEquals(1, CombinatoricsUtils.factorial(0), "0");
        assertEquals(1.0d, CombinatoricsUtils.factorialDouble(0), 1E-14, "0");
        assertEquals(0.0d, CombinatoricsUtils.factorialLog(0), 1E-14, "0");
    }

    @Test
    void testFactorialFail() {
        try {
            CombinatoricsUtils.factorial(-1);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            CombinatoricsUtils.factorialDouble(-1);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            CombinatoricsUtils.factorialLog(-1);
            fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            CombinatoricsUtils.factorial(21);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
            // ignored
        }
        assertTrue(Double.isInfinite(CombinatoricsUtils.factorialDouble(171)), "expecting infinite factorial value");
    }

    @Test
    void testStirlingS2() {

        assertEquals(1, CombinatoricsUtils.stirlingS2(0, 0));

        for (int n = 1; n < 30; ++n) {
            assertEquals(0, CombinatoricsUtils.stirlingS2(n, 0));
            assertEquals(1, CombinatoricsUtils.stirlingS2(n, 1));
            if (n > 2) {
                assertEquals((1l << (n - 1)) - 1l, CombinatoricsUtils.stirlingS2(n, 2));
                assertEquals(CombinatoricsUtils.binomialCoefficient(n, 2),
                                    CombinatoricsUtils.stirlingS2(n, n - 1));
            }
            assertEquals(1, CombinatoricsUtils.stirlingS2(n, n));
        }
        assertEquals(536870911l, CombinatoricsUtils.stirlingS2(30, 2));
        assertEquals(576460752303423487l, CombinatoricsUtils.stirlingS2(60, 2));

        assertEquals(   25, CombinatoricsUtils.stirlingS2( 5, 3));
        assertEquals(   90, CombinatoricsUtils.stirlingS2( 6, 3));
        assertEquals(   65, CombinatoricsUtils.stirlingS2( 6, 4));
        assertEquals(  301, CombinatoricsUtils.stirlingS2( 7, 3));
        assertEquals(  350, CombinatoricsUtils.stirlingS2( 7, 4));
        assertEquals(  140, CombinatoricsUtils.stirlingS2( 7, 5));
        assertEquals(  966, CombinatoricsUtils.stirlingS2( 8, 3));
        assertEquals( 1701, CombinatoricsUtils.stirlingS2( 8, 4));
        assertEquals( 1050, CombinatoricsUtils.stirlingS2( 8, 5));
        assertEquals(  266, CombinatoricsUtils.stirlingS2( 8, 6));
        assertEquals( 3025, CombinatoricsUtils.stirlingS2( 9, 3));
        assertEquals( 7770, CombinatoricsUtils.stirlingS2( 9, 4));
        assertEquals( 6951, CombinatoricsUtils.stirlingS2( 9, 5));
        assertEquals( 2646, CombinatoricsUtils.stirlingS2( 9, 6));
        assertEquals(  462, CombinatoricsUtils.stirlingS2( 9, 7));
        assertEquals( 9330, CombinatoricsUtils.stirlingS2(10, 3));
        assertEquals(34105, CombinatoricsUtils.stirlingS2(10, 4));
        assertEquals(42525, CombinatoricsUtils.stirlingS2(10, 5));
        assertEquals(22827, CombinatoricsUtils.stirlingS2(10, 6));
        assertEquals( 5880, CombinatoricsUtils.stirlingS2(10, 7));
        assertEquals(  750, CombinatoricsUtils.stirlingS2(10, 8));

    }

    @Test
    void testStirlingS2NegativeN() {
        assertThrows(MathIllegalArgumentException.class, () -> CombinatoricsUtils.stirlingS2(3, -1));
    }

    @Test
    void testStirlingS2LargeK() {
        assertThrows(MathIllegalArgumentException.class, () -> CombinatoricsUtils.stirlingS2(3, 4));
    }

    @Test
    void testStirlingS2Overflow() {
        assertThrows(MathRuntimeException.class, () -> CombinatoricsUtils.stirlingS2(26, 9));
    }

    @Test
    void testCheckBinomial1() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            // n < 0
            CombinatoricsUtils.checkBinomial(-1, -2);
        });
    }

    @Test
    void testCheckBinomial2() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            // k > n
            CombinatoricsUtils.checkBinomial(4, 5);
        });
    }

    @Test
    void testCheckBinomial3() {
        // OK (no exception thrown)
        CombinatoricsUtils.checkBinomial(5, 4);
    }

    @Test
    void testBellNumber() {
        // OEIS A000110: http://oeis.org/A000110
        assertEquals(             1l, CombinatoricsUtils.bellNumber( 0));
        assertEquals(             1l, CombinatoricsUtils.bellNumber( 1));
        assertEquals(             2l, CombinatoricsUtils.bellNumber( 2));
        assertEquals(             5l, CombinatoricsUtils.bellNumber( 3));
        assertEquals(            15l, CombinatoricsUtils.bellNumber( 4));
        assertEquals(            52l, CombinatoricsUtils.bellNumber( 5));
        assertEquals(           203l, CombinatoricsUtils.bellNumber( 6));
        assertEquals(           877l, CombinatoricsUtils.bellNumber( 7));
        assertEquals(          4140l, CombinatoricsUtils.bellNumber( 8));
        assertEquals(         21147l, CombinatoricsUtils.bellNumber( 9));
        assertEquals(        115975l, CombinatoricsUtils.bellNumber(10));
        assertEquals(        678570l, CombinatoricsUtils.bellNumber(11));
        assertEquals(       4213597l, CombinatoricsUtils.bellNumber(12));
        assertEquals(      27644437l, CombinatoricsUtils.bellNumber(13));
        assertEquals(     190899322l, CombinatoricsUtils.bellNumber(14));
        assertEquals(    1382958545l, CombinatoricsUtils.bellNumber(15));
        assertEquals(   10480142147l, CombinatoricsUtils.bellNumber(16));
        assertEquals(   82864869804l, CombinatoricsUtils.bellNumber(17));
        assertEquals(  682076806159l, CombinatoricsUtils.bellNumber(18));
        assertEquals( 5832742205057l, CombinatoricsUtils.bellNumber(19));
        assertEquals(51724158235372l, CombinatoricsUtils.bellNumber(20));
    }

    @Test
    void testBellNegative() {
        try {
            CombinatoricsUtils.bellNumber(-1);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.NUMBER_TOO_SMALL, miae.getSpecifier());
            assertEquals(-1, ((Integer) miae.getParts()[0]).intValue());
            assertEquals( 0, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    void testBellLarge() {
        try {
            CombinatoricsUtils.bellNumber(26);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.NUMBER_TOO_LARGE, miae.getSpecifier());
            assertEquals(26, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(25, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    void testPartitions0() {
        List<Integer> emptyList = Collections.emptyList();
        final List<List<Integer>[]> partitions = CombinatoricsUtils.partitions(emptyList).
                                                 collect(Collectors.toList());
        assertEquals(1, partitions.size());
        assertEquals(1, partitions.get(0).length);
        assertEquals(0, partitions.get(0)[0].size());
    }

    @Test
    void testPartitions1() {
        final List<List<Integer>[]> partitions = CombinatoricsUtils.partitions(Arrays.asList(1)).
                                                 collect(Collectors.toList());
        assertEquals(1, partitions.size());
        assertEquals(1, partitions.get(0).length);
        assertEquals(1, partitions.get(0)[0].size());
    }

    @Test
    void testPartitions4() {
        final List<List<Integer>[]> partitions = CombinatoricsUtils.partitions(Arrays.asList(1, 2, 3, 4)).
                                                 collect(Collectors.toList());
        assertEquals(15, partitions.size());

        assertEquals(1, partitions.get(0).length);
        assertArrayEquals(new Integer[] { 1, 2, 3, 4}, partitions.get(0)[0].toArray());

        assertEquals(2, partitions.get(1).length);
        assertArrayEquals(new Integer[] { 1, 2, 3 }, partitions.get(1)[0].toArray());
        assertArrayEquals(new Integer[] { 4 },       partitions.get(1)[1].toArray());

        assertEquals(2, partitions.get(2).length);
        assertArrayEquals(new Integer[] { 1, 2, 4 }, partitions.get(2)[0].toArray());
        assertArrayEquals(new Integer[] { 3 },       partitions.get(2)[1].toArray());

        assertEquals(2, partitions.get(3).length);
        assertArrayEquals(new Integer[] { 1, 2 },    partitions.get(3)[0].toArray());
        assertArrayEquals(new Integer[] { 3, 4 },    partitions.get(3)[1].toArray());

        assertEquals(3, partitions.get(4).length);
        assertArrayEquals(new Integer[] { 1, 2 },    partitions.get(4)[0].toArray());
        assertArrayEquals(new Integer[] { 3 },       partitions.get(4)[1].toArray());
        assertArrayEquals(new Integer[] { 4 },       partitions.get(4)[2].toArray());

        assertEquals(2, partitions.get(5).length);
        assertArrayEquals(new Integer[] { 1, 3, 4 }, partitions.get(5)[0].toArray());
        assertArrayEquals(new Integer[] { 2 },       partitions.get(5)[1].toArray());

        assertEquals(2, partitions.get(6).length);
        assertArrayEquals(new Integer[] { 1, 3 },    partitions.get(6)[0].toArray());
        assertArrayEquals(new Integer[] { 2, 4 },    partitions.get(6)[1].toArray());

        assertEquals(3, partitions.get(7).length);
        assertArrayEquals(new Integer[] { 1, 3 },    partitions.get(7)[0].toArray());
        assertArrayEquals(new Integer[] { 2 },       partitions.get(7)[1].toArray());
        assertArrayEquals(new Integer[] { 4 },       partitions.get(7)[2].toArray());

        assertEquals(2, partitions.get(8).length);
        assertArrayEquals(new Integer[] { 1, 4 },    partitions.get(8)[0].toArray());
        assertArrayEquals(new Integer[] { 2, 3 },    partitions.get(8)[1].toArray());

        assertEquals(2, partitions.get(9).length);
        assertArrayEquals(new Integer[] { 1 },       partitions.get(9)[0].toArray());
        assertArrayEquals(new Integer[] { 2, 3, 4 }, partitions.get(9)[1].toArray());

        assertEquals(3, partitions.get(10).length);
        assertArrayEquals(new Integer[] { 1 },       partitions.get(10)[0].toArray());
        assertArrayEquals(new Integer[] { 2, 3 },    partitions.get(10)[1].toArray());
        assertArrayEquals(new Integer[] { 4 },       partitions.get(10)[2].toArray());

        assertEquals(3, partitions.get(11).length);
        assertArrayEquals(new Integer[] { 1, 4 },    partitions.get(11)[0].toArray());
        assertArrayEquals(new Integer[] { 2 },       partitions.get(11)[1].toArray());
        assertArrayEquals(new Integer[] { 3 },       partitions.get(11)[2].toArray());

        assertEquals(3, partitions.get(12).length);
        assertArrayEquals(new Integer[] { 1 },       partitions.get(12)[0].toArray());
        assertArrayEquals(new Integer[] { 2, 4 },    partitions.get(12)[1].toArray());
        assertArrayEquals(new Integer[] { 3 },       partitions.get(12)[2].toArray());

        assertEquals(3, partitions.get(13).length);
        assertArrayEquals(new Integer[] { 1 },       partitions.get(13)[0].toArray());
        assertArrayEquals(new Integer[] { 2 },       partitions.get(13)[1].toArray());
        assertArrayEquals(new Integer[] { 3, 4},     partitions.get(13)[2].toArray());

        assertEquals(4, partitions.get(14).length);
        assertArrayEquals(new Integer[] { 1 },       partitions.get(14)[0].toArray());
        assertArrayEquals(new Integer[] { 2 },       partitions.get(14)[1].toArray());
        assertArrayEquals(new Integer[] { 3 },       partitions.get(14)[2].toArray());
        assertArrayEquals(new Integer[] { 4 },       partitions.get(14)[3].toArray());

    }

    @Test
    void testPartitions42() {
        final List<List<Integer>[]> partitions = CombinatoricsUtils.partitions(Arrays.asList(1, 2, 3, 4)).
                                                 filter(a -> a.length == 2).
                                                 collect(Collectors.toList());
        assertEquals(7, partitions.size());

        assertEquals(2, partitions.get(0).length);
        assertArrayEquals(new Integer[] { 1, 2, 3 }, partitions.get(0)[0].toArray());
        assertArrayEquals(new Integer[] { 4 },       partitions.get(0)[1].toArray());

        assertEquals(2, partitions.get(1).length);
        assertArrayEquals(new Integer[] { 1, 2, 4 }, partitions.get(1)[0].toArray());
        assertArrayEquals(new Integer[] { 3 },       partitions.get(1)[1].toArray());

        assertEquals(2, partitions.get(2).length);
        assertArrayEquals(new Integer[] { 1, 2 },    partitions.get(2)[0].toArray());
        assertArrayEquals(new Integer[] { 3, 4 },    partitions.get(2)[1].toArray());

        assertEquals(2, partitions.get(3).length);
        assertArrayEquals(new Integer[] { 1, 3, 4 }, partitions.get(3)[0].toArray());
        assertArrayEquals(new Integer[] { 2 },       partitions.get(3)[1].toArray());

        assertEquals(2, partitions.get(4).length);
        assertArrayEquals(new Integer[] { 1, 3 },    partitions.get(4)[0].toArray());
        assertArrayEquals(new Integer[] { 2, 4 },    partitions.get(4)[1].toArray());

        assertEquals(2, partitions.get(5).length);
        assertArrayEquals(new Integer[] { 1, 4 },    partitions.get(5)[0].toArray());
        assertArrayEquals(new Integer[] { 2, 3 },    partitions.get(5)[1].toArray());

        assertEquals(2, partitions.get(6).length);
        assertArrayEquals(new Integer[] { 1 },       partitions.get(6)[0].toArray());
        assertArrayEquals(new Integer[] { 2, 3, 4 }, partitions.get(6)[1].toArray());

    }

    @Test
    void testPartitionsCount() {
        for (int i = 0; i < 12; ++i) {
            List<Integer> list = IntStream.range(0, i).collect(ArrayList::new, ArrayList::add, ArrayList::addAll);
            long partitionsCount = CombinatoricsUtils.partitions(list).count();
            assertEquals(CombinatoricsUtils.bellNumber(i), partitionsCount);
        }
    }

    @Test
    void testExhaustedPartitionsCount() {

        PartitionsIterator<Integer> iterator = new PartitionsIterator<>(Arrays.asList(1, 2, 3));

        assertTrue(iterator.hasNext());
        assertEquals(1, iterator.next().length);
        assertTrue(iterator.hasNext());
        assertEquals(2, iterator.next().length);
        assertTrue(iterator.hasNext());
        assertEquals(2, iterator.next().length);
        assertTrue(iterator.hasNext());
        assertEquals(2, iterator.next().length);
        assertTrue(iterator.hasNext());
        assertEquals(3, iterator.next().length);

        assertFalse(iterator.hasNext());
        try {
            iterator.next();
            fail("an exception should have been thrown");
        } catch (NoSuchElementException e) {
            // expected
        }
    }

    @Test
    void testPermutations0() {
        List<Integer> emptyList = Collections.emptyList();
        final List<List<Integer>> permutations = CombinatoricsUtils.permutations(emptyList).
                                                 collect(Collectors.toList());
        assertEquals(1, permutations.size());
        assertEquals(0, permutations.get(0).size());
    }

    @Test
    void testPermutations1() {
        final List<List<Integer>> permutations = CombinatoricsUtils.permutations(Arrays.asList(1)).
                                                 collect(Collectors.toList());
        assertEquals(1, permutations.size());
        assertEquals(1, permutations.get(0).size());
    }

    @Test
    void testPermutations3() {
        final List<List<Integer>> permutations = CombinatoricsUtils.permutations(Arrays.asList(1, 2, 3)).
                                                 collect(Collectors.toList());
        assertEquals(6, permutations.size());
        assertArrayEquals(new Integer[] { 1, 2, 3 }, permutations.get(0).toArray(new Integer[0]));
        assertArrayEquals(new Integer[] { 1, 3, 2 }, permutations.get(1).toArray(new Integer[0]));
        assertArrayEquals(new Integer[] { 3, 1, 2 }, permutations.get(2).toArray(new Integer[0]));
        assertArrayEquals(new Integer[] { 3, 2, 1 }, permutations.get(3).toArray(new Integer[0]));
        assertArrayEquals(new Integer[] { 2, 3, 1 }, permutations.get(4).toArray(new Integer[0]));
        assertArrayEquals(new Integer[] { 2, 1, 3 }, permutations.get(5).toArray(new Integer[0]));
    }

    @Test
    void testPermutationsCount() {
        for (int i = 0; i < 10; ++i) {
            List<Integer> list = IntStream.range(0, i).collect(ArrayList::new, ArrayList::add, ArrayList::addAll);
            long permutationsCount = CombinatoricsUtils.permutations(list).count();
            assertEquals(CombinatoricsUtils.factorial(i), permutationsCount);
        }
    }

    @Test
    void testExhaustedPermutationsCount() {

        PermutationsIterator<Integer> iterator = new PermutationsIterator<>(Arrays.asList(1, 2, 3));

        assertTrue(iterator.hasNext());
        assertEquals(3, iterator.next().size());
        assertTrue(iterator.hasNext());
        assertEquals(3, iterator.next().size());
        assertTrue(iterator.hasNext());
        assertEquals(3, iterator.next().size());
        assertTrue(iterator.hasNext());
        assertEquals(3, iterator.next().size());
        assertTrue(iterator.hasNext());
        assertEquals(3, iterator.next().size());
        assertTrue(iterator.hasNext());
        assertEquals(3, iterator.next().size());

        assertFalse(iterator.hasNext());
        try {
            iterator.next();
            fail("an exception should have been thrown");
        } catch (NoSuchElementException e) {
            // expected
        }
    }

    /**
     * Exact (caching) recursive implementation to test against
     */
    private long binomialCoefficient(int n, int k) throws MathRuntimeException {
        if (binomialCache.size() > n) {
            Long cachedResult = binomialCache.get(n).get(Integer.valueOf(k));
            if (cachedResult != null) {
                return cachedResult.longValue();
            }
        }
        long result = -1;
        if ((n == k) || (k == 0)) {
            result = 1;
        } else if ((k == 1) || (k == n - 1)) {
            result = n;
        } else {
            // Reduce stack depth for larger values of n
            if (k < n - 100) {
                binomialCoefficient(n - 100, k);
            }
            if (k > 100) {
                binomialCoefficient(n - 100, k - 100);
            }
            result = ArithmeticUtils.addAndCheck(binomialCoefficient(n - 1, k - 1),
                binomialCoefficient(n - 1, k));
        }
        if (result == -1) {
            throw new MathRuntimeException(LocalizedCoreFormats.ARITHMETIC_EXCEPTION);
        }
        for (int i = binomialCache.size(); i < n + 1; i++) {
            binomialCache.add(new HashMap<Integer, Long>());
        }
        binomialCache.get(n).put(Integer.valueOf(k), Long.valueOf(result));
        return result;
    }

    /**
     * Exact direct multiplication implementation to test against
     */
    private long factorial(int n) {
        long result = 1;
        for (int i = 2; i <= n; i++) {
            result *= i;
        }
        return result;
    }
}
