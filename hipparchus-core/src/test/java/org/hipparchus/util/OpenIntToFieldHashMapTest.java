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

import java.util.ConcurrentModificationException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NoSuchElementException;
import java.util.Random;
import java.util.Set;

import org.hipparchus.Field;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.fraction.BigFraction;
import org.hipparchus.fraction.BigFractionField;
import org.hipparchus.fraction.Fraction;
import org.hipparchus.fraction.FractionField;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;


public class OpenIntToFieldHashMapTest {

    private final Map<Integer, Fraction> javaMap = new HashMap<>();
    private final FractionField field = FractionField.getInstance();

    @Before
    public void setUp() {
        javaMap.put(50, new Fraction(100.0));
        javaMap.put(75, new Fraction(75.0));
        javaMap.put(25, new Fraction(500.0));
        javaMap.put(Integer.MAX_VALUE, new Fraction(Integer.MAX_VALUE));
        javaMap.put(0, new Fraction(-1.0));
        javaMap.put(1, new Fraction(0.0));
        javaMap.put(33, new Fraction(-0.1));
        javaMap.put(23234234, new Fraction(-242343.0));
        javaMap.put(23321, new Fraction (Integer.MIN_VALUE));
        javaMap.put(-4444, new Fraction(332.0));
        javaMap.put(-1, new Fraction(-2323.0));
        javaMap.put(Integer.MIN_VALUE, new Fraction(44.0));

        /* Add a few more to cause the table to rehash */
        javaMap.putAll(generate());

    }

    private Map<Integer, Fraction> generate() {
        Map<Integer, Fraction> map = new HashMap<>();
        Random r = new Random();
        double dd;
        for (int i = 0; i < 2000; ++i) {
            dd = r.nextDouble();
            try {
                map.put(r.nextInt(), new Fraction(dd));
            } catch (MathIllegalStateException e) {
                throw new IllegalStateException("Invalid :" + dd, e);
            }
        }
        return map;
    }

    private OpenIntToFieldHashMap<Fraction> createFromJavaMap(Field<Fraction> field) {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(field);
        for (Map.Entry<Integer, Fraction> mapEntry : javaMap.entrySet()) {
            map.put(mapEntry.getKey(), mapEntry.getValue());
        }
        return map;
    }

    @Test
    public void testPutAndGetWith0ExpectedSize() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(field,0);
        assertPutAndGet(map);
    }

    @Test
    public void testPutAndGetWithExpectedSize() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(field,500);
        assertPutAndGet(map);
    }

    @Test
    public void testPutAndGet() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(field);
        assertPutAndGet(map);
    }

    private void assertPutAndGet(OpenIntToFieldHashMap<Fraction> map) {
        assertPutAndGet(map, 0, new HashSet<>());
    }

    private void assertPutAndGet(OpenIntToFieldHashMap<Fraction> map, int mapSize,
            Set<Integer> keysInMap) {
        Assert.assertEquals(mapSize, map.size());
        for (Map.Entry<Integer, Fraction> mapEntry : javaMap.entrySet()) {
            map.put(mapEntry.getKey(), mapEntry.getValue());
            if (!keysInMap.contains(mapEntry.getKey()))
                ++mapSize;
            Assert.assertEquals(mapSize, map.size());
            Assert.assertEquals(mapEntry.getValue(), map.get(mapEntry.getKey()));
        }
    }

    @Test
    public void testPutAbsentOnExisting() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        int size = javaMap.size();
        for (Map.Entry<Integer, Fraction> mapEntry : generateAbsent().entrySet()) {
            map.put(mapEntry.getKey(), mapEntry.getValue());
            Assert.assertEquals(++size, map.size());
            Assert.assertEquals(mapEntry.getValue(), map.get(mapEntry.getKey()));
        }
    }

    @Test
    public void testPutOnExisting() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        for (Map.Entry<Integer, Fraction> mapEntry : javaMap.entrySet()) {
            map.put(mapEntry.getKey(), mapEntry.getValue());
            Assert.assertEquals(javaMap.size(), map.size());
            Assert.assertEquals(mapEntry.getValue(), map.get(mapEntry.getKey()));
        }
    }

    @Test
    public void testGetAbsent() {
        Map<Integer, Fraction> generated = generateAbsent();
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);

        for (Map.Entry<Integer, Fraction> mapEntry : generated.entrySet())
            Assert.assertEquals(field.getZero(), map.get(mapEntry.getKey()));
    }

    @Test
    public void testGetFromEmpty() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(field);
        Assert.assertEquals(field.getZero(), map.get(5));
        Assert.assertEquals(field.getZero(), map.get(0));
        Assert.assertEquals(field.getZero(), map.get(50));
    }

    @Test
    public void testRemove() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        int mapSize = javaMap.size();
        Assert.assertEquals(mapSize, map.size());
        for (Map.Entry<Integer, Fraction> mapEntry : javaMap.entrySet()) {
            map.remove(mapEntry.getKey());
            Assert.assertEquals(--mapSize, map.size());
            Assert.assertEquals(field.getZero(), map.get(mapEntry.getKey()));
        }

        /* Ensure that put and get still work correctly after removals */
        assertPutAndGet(map);
    }

    /* This time only remove some entries */
    @Test
    public void testRemove2() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        int mapSize = javaMap.size();
        int count = 0;
        Set<Integer> keysInMap = new HashSet<>(javaMap.keySet());
        for (Map.Entry<Integer, Fraction> mapEntry : javaMap.entrySet()) {
            keysInMap.remove(mapEntry.getKey());
            map.remove(mapEntry.getKey());
            Assert.assertEquals(--mapSize, map.size());
            Assert.assertEquals(field.getZero(), map.get(mapEntry.getKey()));
            if (count++ > 5)
                break;
        }

        /* Ensure that put and get still work correctly after removals */
        assertPutAndGet(map, mapSize, keysInMap);
    }

    @Test
    public void testRemoveFromEmpty() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(field);
        Assert.assertEquals(field.getZero(), map.remove(50));
    }

    @Test
    public void testRemoveAbsent() {
        Map<Integer, Fraction> generated = generateAbsent();

        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        int mapSize = map.size();

        for (Map.Entry<Integer, Fraction> mapEntry : generated.entrySet()) {
            map.remove(mapEntry.getKey());
            Assert.assertEquals(mapSize, map.size());
            Assert.assertEquals(field.getZero(), map.get(mapEntry.getKey()));
        }
    }

    /**
     * Returns a map with at least 100 elements where each element is absent from javaMap.
     */
    private Map<Integer, Fraction> generateAbsent() {
        Map<Integer, Fraction> generated = new HashMap<>();
        do {
            generated.putAll(generate());
            for (Integer key : javaMap.keySet())
                generated.remove(key);
        } while (generated.size() < 100);
        return generated;
    }

    @Test
    public void testCopy() {
        OpenIntToFieldHashMap<Fraction> copy =
            new OpenIntToFieldHashMap<>(createFromJavaMap(field));
        Assert.assertEquals(javaMap.size(), copy.size());

        for (Map.Entry<Integer, Fraction> mapEntry : javaMap.entrySet())
            Assert.assertEquals(mapEntry.getValue(), copy.get(mapEntry.getKey()));
    }

    @Test
    public void testContainsKey() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        for (Entry<Integer, Fraction> mapEntry : javaMap.entrySet()) {
            Assert.assertTrue(map.containsKey(mapEntry.getKey()));
        }
        for (Map.Entry<Integer, Fraction> mapEntry : generateAbsent().entrySet()) {
            Assert.assertFalse(map.containsKey(mapEntry.getKey()));
        }
        for (Entry<Integer, Fraction> mapEntry : javaMap.entrySet()) {
            int key = mapEntry.getKey();
            Assert.assertTrue(map.containsKey(key));
            map.remove(key);
            Assert.assertFalse(map.containsKey(key));
        }
    }

    @Test
    public void testIterator() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        OpenIntToFieldHashMap<Fraction>.Iterator iterator = map.iterator();
        for (int i = 0; i < map.size(); ++i) {
            Assert.assertTrue(iterator.hasNext());
            iterator.advance();
            int key = iterator.key();
            Assert.assertTrue(map.containsKey(key));
            Assert.assertEquals(javaMap.get(key), map.get(key));
            Assert.assertEquals(javaMap.get(key), iterator.value());
            Assert.assertTrue(javaMap.containsKey(key));
        }
        Assert.assertFalse(iterator.hasNext());
        try {
            iterator.advance();
            Assert.fail("an exception should have been thrown");
        } catch (NoSuchElementException nsee) {
            // expected
        }
    }

    @Test
    public void testEquals() {
        OpenIntToFieldHashMap<Fraction> map1 = new OpenIntToFieldHashMap<>(FractionField.getInstance());
        map1.put(2,   new Fraction( 5, 2));
        map1.put(17,  new Fraction(-1, 2));
        map1.put(16,  Fraction.ZERO);
        Assert.assertEquals(map1, map1);
        OpenIntToFieldHashMap<Fraction> map2 = new OpenIntToFieldHashMap<>(FractionField.getInstance());
        map2.put(17, new Fraction(-1, 2));
        map2.put(2,  new Fraction( 5, 2));
        map2.put(16, new Fraction(0));
        Assert.assertEquals(map1, map2);
        map2.put(16,  new Fraction( 1, 4));
        Assert.assertNotEquals(map1, map2);
        map2.put(16,  Fraction.ZERO);
        Assert.assertEquals(map1, map2);
        OpenIntToFieldHashMap<BigFraction> map3 = new OpenIntToFieldHashMap<>(BigFractionField.getInstance());
        map3.put(2,   new BigFraction( 5, 2));
        map3.put(17,  new BigFraction(-1, 2));
        map3.put(16,  BigFraction.ZERO);
        Assert.assertNotEquals(map1, map3);
        Assert.assertNotEquals(map1, "");
        Assert.assertNotEquals(map1, null);
    }

    @Test
    public void testHashcode() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(FractionField.getInstance());
        map.put(2,   new Fraction( 5, 2));
        map.put(17,  new Fraction(-1, 2));
        map.put(16,  Fraction.ZERO);
        Assert.assertEquals(-528348218, map.hashCode());
    }

    @Test
    public void testConcurrentModification() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        OpenIntToFieldHashMap<Fraction>.Iterator iterator = map.iterator();
        map.put(3, new Fraction(3));
        try {
            iterator.advance();
            Assert.fail("an exception should have been thrown");
        } catch (ConcurrentModificationException cme) {
            // expected
        }
    }

    /**
     * Regression test for a bug in findInsertionIndex where the hashing in the second probing
     * loop was inconsistent with the first causing duplicate keys after the right sequence
     * of puts and removes.
     */
    @Test
    public void testPutKeysWithCollisions() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(field);
        int key1 = -1996012590;
        Fraction value1 = new Fraction(1);
        map.put(key1, value1);
        int key2 = 835099822;
        map.put(key2, value1);
        int key3 = 1008859686;
        map.put(key3, value1);
        Assert.assertEquals(value1, map.get(key3));
        Assert.assertEquals(3, map.size());

        map.remove(key2);
        Fraction value2 = new Fraction(2);
        map.put(key3, value2);
        Assert.assertEquals(value2, map.get(key3));
        Assert.assertEquals(2, map.size());
    }

    /**
     * Similar to testPutKeysWithCollisions() but exercises the codepaths in a slightly
     * different manner.
     */
    @Test
    public void testPutKeysWithCollision2() {
        OpenIntToFieldHashMap<Fraction>map = new OpenIntToFieldHashMap<>(field);
        int key1 = 837989881;
        Fraction value1 = new Fraction(1);
        map.put(key1, value1);
        int key2 = 476463321;
        map.put(key2, value1);
        Assert.assertEquals(2, map.size());
        Assert.assertEquals(value1, map.get(key2));

        map.remove(key1);
        Fraction value2 = new Fraction(2);
        map.put(key2, value2);
        Assert.assertEquals(1, map.size());
        Assert.assertEquals(value2, map.get(key2));
    }


}
