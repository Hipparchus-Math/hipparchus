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
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Comparator;
import java.util.Iterator;
import java.util.NoSuchElementException;

import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Tests for the {@link Combinations} class.
 */
public class CombinationsTest {
    @Test
    public void testAccessor1() {
        final int n = 5;
        final int k = 3;
        Assertions.assertEquals(n, new Combinations(n, k).getN());
    }
    @Test
    public void testAccessor2() {
        final int n = 5;
        final int k = 3;
        Assertions.assertEquals(k, new Combinations(n, k).getK());
    }

    @Test
    public void testLexicographicIterator() {
        checkLexicographicIterator(new Combinations(5, 3));
        checkLexicographicIterator(new Combinations(6, 4));
        checkLexicographicIterator(new Combinations(8, 2));
        checkLexicographicIterator(new Combinations(6, 1));
        checkLexicographicIterator(new Combinations(3, 3));
        checkLexicographicIterator(new Combinations(1, 1));
        checkLexicographicIterator(new Combinations(1, 0));
        checkLexicographicIterator(new Combinations(0, 0));
        checkLexicographicIterator(new Combinations(4, 2));
        checkLexicographicIterator(new Combinations(123, 2));
    }

    @Test
    public void testLexicographicComparatorWrongIterate1() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final int n = 5;
            final int k = 3;
            final Comparator<int[]> comp = new Combinations(n, k).comparator();
            comp.compare(new int[]{1}, new int[]{0, 1, 2});
        });
    }

    @Test
    public void testLexicographicComparatorWrongIterate2() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final int n = 5;
            final int k = 3;
            final Comparator<int[]> comp = new Combinations(n, k).comparator();
            comp.compare(new int[]{0, 1, 2}, new int[]{0, 1, 2, 3});
        });
    }

    @Test
    public void testLexicographicComparatorWrongIterate3() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final int n = 5;
            final int k = 3;
            final Comparator<int[]> comp = new Combinations(n, k).comparator();
            comp.compare(new int[]{1, 2, 5}, new int[]{0, 1, 2});
        });
    }

    @Test
    public void testLexicographicComparatorWrongIterate4() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final int n = 5;
            final int k = 3;
            final Comparator<int[]> comp = new Combinations(n, k).comparator();
            comp.compare(new int[]{1, 2, 4}, new int[]{-1, 1, 2});
        });
    }

    @Test
    public void testLexicographicComparator() {
        final int n = 5;
        final int k = 3;
        final Comparator<int[]> comp = new Combinations(n, k).comparator();
        Assertions.assertEquals(1, comp.compare(new int[] {1, 2, 4},
                                            new int[] {1, 2, 3}));
        Assertions.assertEquals(-1, comp.compare(new int[] {0, 1, 4},
                                             new int[] {0, 2, 4}));
        Assertions.assertEquals(0, comp.compare(new int[] {1, 3, 4},
                                            new int[] {1, 3, 4}));
    }

    /**
     * Check that iterates can be passed unsorted.
     */
    @Test
    public void testLexicographicComparatorUnsorted() {
        final int n = 5;
        final int k = 3;
        final Comparator<int[]> comp = new Combinations(n, k).comparator();
        Assertions.assertEquals(1, comp.compare(new int[] {1, 4, 2},
                                            new int[] {1, 3, 2}));
        Assertions.assertEquals(-1, comp.compare(new int[] {0, 4, 1},
                                             new int[] {0, 4, 2}));
        Assertions.assertEquals(0, comp.compare(new int[] {1, 4, 3},
                                            new int[] {1, 3, 4}));
    }

    @Test
    public void testEmptyCombination() {
        final Iterator<int[]> iter = new Combinations(12345, 0).iterator();
        Assertions.assertTrue(iter.hasNext());
        final int[] c = iter.next();
        Assertions.assertEquals(0, c.length);
        Assertions.assertFalse(iter.hasNext());
    }

    @Test
    public void testFullSetCombination() {
        final int n = 67;
        final Iterator<int[]> iter = new Combinations(n, n).iterator();
        Assertions.assertTrue(iter.hasNext());
        final int[] c = iter.next();
        Assertions.assertEquals(n, c.length);

        for (int i = 0; i < n; i++) {
            Assertions.assertEquals(i, c[i]);
        }

        Assertions.assertFalse(iter.hasNext());
    }

    /**
     * Verifies that the iterator generates a lexicographically
     * increasing sequence of b(n,k) arrays, each having length k
     * and each array itself increasing.
     *
     * @param c Combinations.
     */
    private void checkLexicographicIterator(Combinations c) {
        final Comparator<int[]> comp = c.comparator();
        final int n = c.getN();
        final int k = c.getK();

        int[] lastIterate = null;

        long numIterates = 0;
        for (int[] iterate : c) {
            Assertions.assertEquals(k, iterate.length);

            // Check that the sequence of iterates is ordered.
            if (lastIterate != null) {
                Assertions.assertEquals(1, comp.compare(iterate, lastIterate));
            }

            // Check that each iterate is ordered.
            for (int i = 1; i < iterate.length; i++) {
                Assertions.assertTrue(iterate[i] > iterate[i - 1]);
            }

            lastIterate = iterate;
            ++numIterates;
        }

        // Check the number of iterates.
        Assertions.assertEquals(CombinatoricsUtils.binomialCoefficient(n, k),
                            numIterates);
    }

    @Test
    public void testCombinationsIteratorFail() {
        try {
            new Combinations(4, 5).iterator();
            Assertions.fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }

        try {
            new Combinations(-1, -2).iterator();
            Assertions.fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testLexicographicIteratorUnreachable() {
        // this tests things that could really never happen,
        // as the conditions are tested in the enclosing class before
        // the lexicographic iterator is built
        try {
            Class<?> lexicographicIteratorClass = null;
            for (Class<?> c : Combinations.class.getDeclaredClasses()) {
                if (c.getCanonicalName().endsWith(".LexicographicIterator")) {
                    lexicographicIteratorClass = c;
                }
            }
            Constructor<?> ctr = lexicographicIteratorClass.getDeclaredConstructor(Integer.TYPE, Integer.TYPE);
            Method hasNext = lexicographicIteratorClass.getDeclaredMethod("hasNext");
            Method next = lexicographicIteratorClass.getDeclaredMethod("next");
            Method remove = lexicographicIteratorClass.getDeclaredMethod("remove");
            Assertions.assertFalse((Boolean) hasNext.invoke(ctr.newInstance(3, 0)));
            Assertions.assertFalse((Boolean) hasNext.invoke(ctr.newInstance(3, 3)));
            Assertions.assertTrue((Boolean) hasNext.invoke(ctr.newInstance(3, 2)));
            try {
                next.invoke(ctr.newInstance(3, 0));
                Assertions.fail("an exception should have been thrown");
            } catch (InvocationTargetException ite) {
                Assertions.assertTrue(ite.getCause() instanceof NoSuchElementException);
            }
            try {
                remove.invoke(ctr.newInstance(3, 2));
                Assertions.fail("an exception should have been thrown");
            } catch (InvocationTargetException ite) {
                Assertions.assertTrue(ite.getCause() instanceof UnsupportedOperationException);
            }
        } catch (NoSuchMethodException | IllegalAccessException | IllegalArgumentException |
                 InvocationTargetException | InstantiationException e) {
            Assertions.fail(e.getLocalizedMessage());
        }
    }

    @Test
    public void testSingletonIteratorUnreachable() {
        // this tests things that could really never happen,
        try {
            Class<?> singletonIteratorClass = null;
            for (Class<?> c : Combinations.class.getDeclaredClasses()) {
                if (c.getCanonicalName().endsWith(".SingletonIterator")) {
                    singletonIteratorClass = c;
                }
            }
            Constructor<?> ctr = singletonIteratorClass.getDeclaredConstructor(Integer.TYPE);
            Method hasNext = singletonIteratorClass.getDeclaredMethod("hasNext");
            Method next = singletonIteratorClass.getDeclaredMethod("next");
            Method remove = singletonIteratorClass.getDeclaredMethod("remove");
            Object iterator = ctr.newInstance(3);
            Assertions.assertTrue((Boolean) hasNext.invoke(iterator));
            int[] ret = (int[]) next.invoke(iterator);
            Assertions.assertEquals(3, ret.length);
            Assertions.assertEquals(0, ret[0]);
            Assertions.assertEquals(1, ret[1]);
            Assertions.assertEquals(2, ret[2]);
            Assertions.assertFalse((Boolean) hasNext.invoke(iterator));
            try {
                next.invoke(iterator);
                Assertions.fail("an exception should have been thrown");
            } catch (InvocationTargetException ite) {
                Assertions.assertTrue(ite.getCause() instanceof NoSuchElementException);
            }
            try {
                remove.invoke(iterator);
                Assertions.fail("an exception should have been thrown");
            } catch (InvocationTargetException ite) {
                Assertions.assertTrue(ite.getCause() instanceof UnsupportedOperationException);
            }
        } catch (NoSuchMethodException | IllegalAccessException | IllegalArgumentException |
                 InvocationTargetException | InstantiationException e) {
            Assertions.fail(e.getLocalizedMessage());
        }
    }

}
