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

import java.util.NoSuchElementException;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test cases for the {@link RosenNumberPartitionIterator} class.
 */
public class RosenNumberPartitionIteratorTest {

    @Test
    public void testRosenPartitionNegativeK() {
        try {
         new RosenNumberPartitionIterator(4, -1);
         Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, miae.getSpecifier());
            Assert.assertEquals(-1, ((Integer) miae.getParts()[0]).intValue());
            Assert.assertEquals( 1, ((Integer) miae.getParts()[1]).intValue());
            Assert.assertEquals( 4, ((Integer) miae.getParts()[2]).intValue());
        }
    }

    @Test
    public void testRosenPartitionKGreaterThanN() {
        try {
         new RosenNumberPartitionIterator(4, 5);
         Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, miae.getSpecifier());
            Assert.assertEquals( 5, ((Integer) miae.getParts()[0]).intValue());
            Assert.assertEquals( 1, ((Integer) miae.getParts()[1]).intValue());
            Assert.assertEquals( 4, ((Integer) miae.getParts()[2]).intValue());
        }
    }

    @Test
    public void testRosenPartition42() {
        RosenNumberPartitionIterator i = new RosenNumberPartitionIterator(4, 2);
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 1,  3 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 2,  2 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 3, 1 }, i.next());
        Assert.assertFalse(i.hasNext());
        try {
            i.next();
            Assert.fail("an exception should have been thrown");
        } catch (NoSuchElementException e) {
            // expected
        }
    }

    @Test
    public void testRosenPartition103() {
        RosenNumberPartitionIterator i = new RosenNumberPartitionIterator(10, 3);
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 1, 1, 8 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 1, 2, 7 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 1, 3, 6 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 1, 4, 5 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 1, 5, 4 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 1, 6, 3 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 1, 7, 2 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 1, 8, 1 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 2, 1, 7 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 2, 2, 6 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 2, 3, 5 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 2, 4, 4 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 2, 5, 3 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 2, 6, 2 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 2, 7, 1 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 3, 1, 6 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 3, 2, 5 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 3, 3, 4 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 3, 4, 3 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 3, 5, 2 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 3, 6, 1 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 4, 1, 5 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 4, 2, 4 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 4, 3, 3 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 4, 4, 2 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 4, 5, 1 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 5, 1, 4 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 5, 2, 3 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 5, 3, 2 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 5, 4, 1 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 6, 1, 3 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 6, 2, 2 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 6, 3, 1 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 7, 1, 2 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 7, 2, 1 }, i.next());
        Assert.assertTrue(i.hasNext());
        Assert.assertArrayEquals(new int[] { 8, 1, 1 }, i.next());
        Assert.assertFalse(i.hasNext());
    }

}
