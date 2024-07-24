/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.util;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.NoSuchElementException;

/**
 * Test cases for the {@link RosenNumberPartitionIterator} class.
 */
public class RosenNumberPartitionIteratorTest {

    @Test
    public void testRosenPartitionNegativeK() {
        try {
         new RosenNumberPartitionIterator(4, -1);
         Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, miae.getSpecifier());
            Assertions.assertEquals(-1, ((Integer) miae.getParts()[0]).intValue());
            Assertions.assertEquals( 1, ((Integer) miae.getParts()[1]).intValue());
            Assertions.assertEquals( 4, ((Integer) miae.getParts()[2]).intValue());
        }
    }

    @Test
    public void testRosenPartitionKGreaterThanN() {
        try {
         new RosenNumberPartitionIterator(4, 5);
         Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, miae.getSpecifier());
            Assertions.assertEquals( 5, ((Integer) miae.getParts()[0]).intValue());
            Assertions.assertEquals( 1, ((Integer) miae.getParts()[1]).intValue());
            Assertions.assertEquals( 4, ((Integer) miae.getParts()[2]).intValue());
        }
    }

    @Test
    public void testRosenPartition42() {
        RosenNumberPartitionIterator i = new RosenNumberPartitionIterator(4, 2);
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 1,  3 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 2,  2 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 3, 1 }, i.next());
        Assertions.assertFalse(i.hasNext());
        try {
            i.next();
            Assertions.fail("an exception should have been thrown");
        } catch (NoSuchElementException e) {
            // expected
        }
    }

    @Test
    public void testRosenPartition103() {
        RosenNumberPartitionIterator i = new RosenNumberPartitionIterator(10, 3);
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 1, 1, 8 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 1, 2, 7 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 1, 3, 6 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 1, 4, 5 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 1, 5, 4 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 1, 6, 3 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 1, 7, 2 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 1, 8, 1 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 2, 1, 7 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 2, 2, 6 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 2, 3, 5 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 2, 4, 4 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 2, 5, 3 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 2, 6, 2 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 2, 7, 1 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 3, 1, 6 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 3, 2, 5 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 3, 3, 4 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 3, 4, 3 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 3, 5, 2 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 3, 6, 1 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 4, 1, 5 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 4, 2, 4 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 4, 3, 3 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 4, 4, 2 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 4, 5, 1 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 5, 1, 4 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 5, 2, 3 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 5, 3, 2 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 5, 4, 1 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 6, 1, 3 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 6, 2, 2 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 6, 3, 1 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 7, 1, 2 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 7, 2, 1 }, i.next());
        Assertions.assertTrue(i.hasNext());
        Assertions.assertArrayEquals(new int[] { 8, 1, 1 }, i.next());
        Assertions.assertFalse(i.hasNext());
    }

}
