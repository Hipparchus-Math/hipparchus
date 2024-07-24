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

import org.hipparchus.Field;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.fraction.BigFraction;
import org.hipparchus.fraction.BigFractionField;
import org.hipparchus.fraction.Fraction;
import org.hipparchus.fraction.FractionField;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ConcurrentModificationException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NoSuchElementException;
import java.util.Random;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;


class OpenIntToFieldHashMapTest {

    private final Map<Integer, Fraction> javaMap = new HashMap<>();
    private final FractionField field = FractionField.getInstance();

    @BeforeEach
    void setUp() {
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
    void testPutAndGetWith0ExpectedSize() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(field,0);
        customAssertPutAndGet(map);
    }

    @Test
    void testPutAndGetWithExpectedSize() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(field,500);
        customAssertPutAndGet(map);
    }

    @Test
    void testPutAndGet() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(field);
        customAssertPutAndGet(map);
    }

    private void customAssertPutAndGet(OpenIntToFieldHashMap<Fraction> map) {
        customAssertPutAndGet(map, 0, new HashSet<>());
    }

    private void customAssertPutAndGet(OpenIntToFieldHashMap<Fraction> map, int mapSize,
                                       Set<Integer> keysInMap) {
        assertEquals(mapSize, map.size());
        for (Map.Entry<Integer, Fraction> mapEntry : javaMap.entrySet()) {
            map.put(mapEntry.getKey(), mapEntry.getValue());
            if (!keysInMap.contains(mapEntry.getKey()))
                ++mapSize;
            assertEquals(mapSize, map.size());
            assertEquals(mapEntry.getValue(), map.get(mapEntry.getKey()));
        }
    }

    @Test
    void testPutAbsentOnExisting() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        int size = javaMap.size();
        for (Map.Entry<Integer, Fraction> mapEntry : generateAbsent().entrySet()) {
            map.put(mapEntry.getKey(), mapEntry.getValue());
            assertEquals(++size, map.size());
            assertEquals(mapEntry.getValue(), map.get(mapEntry.getKey()));
        }
    }

    @Test
    void testPutOnExisting() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        for (Map.Entry<Integer, Fraction> mapEntry : javaMap.entrySet()) {
            map.put(mapEntry.getKey(), mapEntry.getValue());
            assertEquals(javaMap.size(), map.size());
            assertEquals(mapEntry.getValue(), map.get(mapEntry.getKey()));
        }
    }

    @Test
    void testGetAbsent() {
        Map<Integer, Fraction> generated = generateAbsent();
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);

        for (Map.Entry<Integer, Fraction> mapEntry : generated.entrySet())
            assertEquals(field.getZero(), map.get(mapEntry.getKey()));
    }

    @Test
    void testGetFromEmpty() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(field);
        assertEquals(field.getZero(), map.get(5));
        assertEquals(field.getZero(), map.get(0));
        assertEquals(field.getZero(), map.get(50));
    }

    @Test
    void testRemove() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        int mapSize = javaMap.size();
        assertEquals(mapSize, map.size());
        for (Map.Entry<Integer, Fraction> mapEntry : javaMap.entrySet()) {
            map.remove(mapEntry.getKey());
            assertEquals(--mapSize, map.size());
            assertEquals(field.getZero(), map.get(mapEntry.getKey()));
        }

        /* Ensure that put and get still work correctly after removals */
        customAssertPutAndGet(map);
    }

    /* This time only remove some entries */
    @Test
    void testRemove2() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        int mapSize = javaMap.size();
        int count = 0;
        Set<Integer> keysInMap = new HashSet<>(javaMap.keySet());
        for (Map.Entry<Integer, Fraction> mapEntry : javaMap.entrySet()) {
            keysInMap.remove(mapEntry.getKey());
            map.remove(mapEntry.getKey());
            assertEquals(--mapSize, map.size());
            assertEquals(field.getZero(), map.get(mapEntry.getKey()));
            if (count++ > 5)
                break;
        }

        /* Ensure that put and get still work correctly after removals */
        customAssertPutAndGet(map, mapSize, keysInMap);
    }

    @Test
    void testRemoveFromEmpty() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(field);
        assertEquals(field.getZero(), map.remove(50));
    }

    @Test
    void testRemoveAbsent() {
        Map<Integer, Fraction> generated = generateAbsent();

        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        int mapSize = map.size();

        for (Map.Entry<Integer, Fraction> mapEntry : generated.entrySet()) {
            map.remove(mapEntry.getKey());
            assertEquals(mapSize, map.size());
            assertEquals(field.getZero(), map.get(mapEntry.getKey()));
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
    void testCopy() {
        OpenIntToFieldHashMap<Fraction> copy =
            new OpenIntToFieldHashMap<>(createFromJavaMap(field));
        assertEquals(javaMap.size(), copy.size());

        for (Map.Entry<Integer, Fraction> mapEntry : javaMap.entrySet())
            assertEquals(mapEntry.getValue(), copy.get(mapEntry.getKey()));
    }

    @Test
    void testContainsKey() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        for (Entry<Integer, Fraction> mapEntry : javaMap.entrySet()) {
            assertTrue(map.containsKey(mapEntry.getKey()));
        }
        for (Map.Entry<Integer, Fraction> mapEntry : generateAbsent().entrySet()) {
            assertFalse(map.containsKey(mapEntry.getKey()));
        }
        for (Entry<Integer, Fraction> mapEntry : javaMap.entrySet()) {
            int key = mapEntry.getKey();
            assertTrue(map.containsKey(key));
            map.remove(key);
            assertFalse(map.containsKey(key));
        }
    }

    @Test
    void testIterator() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        OpenIntToFieldHashMap<Fraction>.Iterator iterator = map.iterator();
        for (int i = 0; i < map.size(); ++i) {
            assertTrue(iterator.hasNext());
            iterator.advance();
            int key = iterator.key();
            assertTrue(map.containsKey(key));
            assertEquals(javaMap.get(key), map.get(key));
            assertEquals(javaMap.get(key), iterator.value());
            assertTrue(javaMap.containsKey(key));
        }
        assertFalse(iterator.hasNext());
        try {
            iterator.advance();
            fail("an exception should have been thrown");
        } catch (NoSuchElementException nsee) {
            // expected
        }
    }

    @Test
    void testEquals() {
        OpenIntToFieldHashMap<Fraction> map1 = new OpenIntToFieldHashMap<>(FractionField.getInstance());
        map1.put(2,   new Fraction( 5, 2));
        map1.put(17,  new Fraction(-1, 2));
        map1.put(16,  Fraction.ZERO);
        assertEquals(map1, map1);
        OpenIntToFieldHashMap<Fraction> map2 = new OpenIntToFieldHashMap<>(FractionField.getInstance());
        map2.put(17, new Fraction(-1, 2));
        map2.put(2,  new Fraction( 5, 2));
        map2.put(16, new Fraction(0));
        assertEquals(map1, map2);
        map2.put(16,  new Fraction( 1, 4));
        assertNotEquals(map1, map2);
        map2.put(16,  Fraction.ZERO);
        assertEquals(map1, map2);
        OpenIntToFieldHashMap<BigFraction> map3 = new OpenIntToFieldHashMap<>(BigFractionField.getInstance());
        map3.put(2,   new BigFraction( 5, 2));
        map3.put(17,  new BigFraction(-1, 2));
        map3.put(16,  BigFraction.ZERO);
        assertNotEquals(map1, map3);
        assertNotEquals("", map1);
        assertNotEquals(null, map1);
    }

    @Test
    void testHashcode() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(FractionField.getInstance());
        map.put(2,   new Fraction( 5, 2));
        map.put(17,  new Fraction(-1, 2));
        map.put(16,  Fraction.ZERO);
        assertEquals(-528348218, map.hashCode());
    }

    @Test
    void testConcurrentModification() {
        OpenIntToFieldHashMap<Fraction> map = createFromJavaMap(field);
        OpenIntToFieldHashMap<Fraction>.Iterator iterator = map.iterator();
        map.put(3, new Fraction(3));
        try {
            iterator.advance();
            fail("an exception should have been thrown");
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
    void testPutKeysWithCollisions() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<>(field);
        int key1 = -1996012590;
        Fraction value1 = new Fraction(1);
        map.put(key1, value1);
        int key2 = 835099822;
        map.put(key2, value1);
        int key3 = 1008859686;
        map.put(key3, value1);
        assertEquals(value1, map.get(key3));
        assertEquals(3, map.size());

        map.remove(key2);
        Fraction value2 = new Fraction(2);
        map.put(key3, value2);
        assertEquals(value2, map.get(key3));
        assertEquals(2, map.size());
    }

    /**
     * Similar to testPutKeysWithCollisions() but exercises the codepaths in a slightly
     * different manner.
     */
    @Test
    void testPutKeysWithCollision2() {
        OpenIntToFieldHashMap<Fraction>map = new OpenIntToFieldHashMap<>(field);
        int key1 = 837989881;
        Fraction value1 = new Fraction(1);
        map.put(key1, value1);
        int key2 = 476463321;
        map.put(key2, value1);
        assertEquals(2, map.size());
        assertEquals(value1, map.get(key2));

        map.remove(key1);
        Fraction value2 = new Fraction(2);
        map.put(key2, value2);
        assertEquals(1, map.size());
        assertEquals(value2, map.get(key2));
    }


}
