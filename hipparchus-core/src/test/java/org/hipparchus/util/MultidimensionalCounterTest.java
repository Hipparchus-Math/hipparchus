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

import java.util.NoSuchElementException;

import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Test cases for the {@link MultidimensionalCounter} class.
 */
public class MultidimensionalCounterTest {
    @Test
    public void testPreconditions() {
        MultidimensionalCounter c;

        try {
            c = new MultidimensionalCounter(0, 1);
            Assertions.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException e) {
            // Expected.
        }
        try {
            c = new MultidimensionalCounter(2, 0);
            Assertions.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException e) {
            // Expected.
        }
        try {
            c = new MultidimensionalCounter(-1, 1);
            Assertions.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException e) {
            // Expected.
        }

        c = new MultidimensionalCounter(2, 3);
        try {
            c.getCount(1, 1, 1);
            Assertions.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException e) {
            // Expected.
        }
        try {
            c.getCount(3, 1);
            Assertions.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException e) {
            // Expected.
        }
        try {
            c.getCount(0, -1);
            Assertions.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException e) {
            // Expected.
        }
        try {
            c.getCounts(-1);
            Assertions.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException e) {
            // Expected.
        }
        try {
            c.getCounts(6);
            Assertions.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException e) {
            // Expected.
        }
    }

    @Test
    public void testIteratorPreconditions() {
        MultidimensionalCounter.Iterator iter = (new MultidimensionalCounter(2, 3)).iterator();
        try {
            iter.getCount(-1);
            Assertions.fail("IndexOutOfBoundsException expected");
        } catch (IndexOutOfBoundsException e) {
            // Expected.
        }
        try {
            iter.getCount(2);
            Assertions.fail("IndexOutOfBoundsException expected");
        } catch (IndexOutOfBoundsException e) {
            // Expected.
        }
    }

    @Test
    public void testIterator() {
        final int dim1 = 3;
        final int dim2 = 4;

        final MultidimensionalCounter.Iterator iter
            = new MultidimensionalCounter(dim1, dim2).iterator();

        final int max = dim1 * dim2;
        for (int i = 0; i < max; i++) {
            Assertions.assertTrue(iter.hasNext());

            // Should not throw.
            iter.next();
        }

        Assertions.assertFalse(iter.hasNext());
    }

    @Test
    public void testIteratorNoMoreElements() {
        assertThrows(NoSuchElementException.class, () -> {
            final MultidimensionalCounter.Iterator iter
                = new MultidimensionalCounter(4, 2).iterator();

            while (iter.hasNext()) {
                iter.next();
            }

            // No more elements: should throw.
            iter.next();
        });
    }

    @Test
    public void testMulti2UniConversion() {
        final MultidimensionalCounter c = new MultidimensionalCounter(2, 4, 5);
        Assertions.assertEquals(33, c.getCount(1, 2, 3));
    }

    @Test
    public void testAccessors() {
        final int[] originalSize = new int[] {2, 6, 5};
        final MultidimensionalCounter c = new MultidimensionalCounter(originalSize);
        final int nDim = c.getDimension();
        Assertions.assertEquals(nDim, originalSize.length);

        final int[] size = c.getSizes();
        for (int i = 0; i < nDim; i++) {
            Assertions.assertEquals(originalSize[i], size[i]);
        }
    }

    @Test
    public void testIterationConsistency() {
        final MultidimensionalCounter c = new MultidimensionalCounter(2, 3, 4);
        final int[][] expected = new int[][] {
            { 0, 0, 0 },
            { 0, 0, 1 },
            { 0, 0, 2 },
            { 0, 0, 3 },
            { 0, 1, 0 },
            { 0, 1, 1 },
            { 0, 1, 2 },
            { 0, 1, 3 },
            { 0, 2, 0 },
            { 0, 2, 1 },
            { 0, 2, 2 },
            { 0, 2, 3 },
            { 1, 0, 0 },
            { 1, 0, 1 },
            { 1, 0, 2 },
            { 1, 0, 3 },
            { 1, 1, 0 },
            { 1, 1, 1 },
            { 1, 1, 2 },
            { 1, 1, 3 },
            { 1, 2, 0 },
            { 1, 2, 1 },
            { 1, 2, 2 },
            { 1, 2, 3 }
        };

        final int totalSize = c.getSize();
        Assertions.assertEquals(expected.length, totalSize);

        final int nDim = c.getDimension();
        final MultidimensionalCounter.Iterator iter = c.iterator();
        for (int i = 0; i < totalSize; i++) {
            if (!iter.hasNext()) {
                Assertions.fail("Too short");
            }
            final int uniDimIndex = iter.next().intValue();
            Assertions.assertEquals(i, uniDimIndex, "Wrong iteration at " + i);

            for (int dimIndex = 0; dimIndex < nDim; dimIndex++) {
                Assertions.assertEquals(expected[i][dimIndex], iter.getCount(dimIndex), "Wrong multidimensional index for [" + i + "][" + dimIndex + "]");
            }

            Assertions.assertEquals(c.getCount(expected[i]), uniDimIndex, "Wrong unidimensional index for [" + i + "]");

            final int[] indices = c.getCounts(uniDimIndex);
            for (int dimIndex = 0; dimIndex < nDim; dimIndex++) {
                Assertions.assertEquals(expected[i][dimIndex], indices[dimIndex], "Wrong multidimensional index for [" + i + "][" + dimIndex + "]");
            }
        }

        if (iter.hasNext()) {
            Assertions.fail("Too long");
        }
    }
}
