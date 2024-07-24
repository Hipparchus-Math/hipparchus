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

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ConcurrentModificationException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Random;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;


/**
 * Test cases for the {@link OpenIntToDoubleHashMap}.
 */
class OpenIntToDoubleHashMapTest {

    private final Map<Integer, Double> javaMap = new HashMap<>();

    @BeforeEach
    void setUp() throws Exception {
        javaMap.put(50, 100.0);
        javaMap.put(75, 75.0);
        javaMap.put(25, 500.0);
        javaMap.put(Integer.MAX_VALUE, Double.MAX_VALUE);
        javaMap.put(0, -1.0);
        javaMap.put(1, 0.0);
        javaMap.put(33, -0.1);
        javaMap.put(23234234, -242343.0);
        javaMap.put(23321, Double.MIN_VALUE);
        javaMap.put(-4444, 332.0);
        javaMap.put(-1, -2323.0);
        javaMap.put(Integer.MIN_VALUE, 44.0);

        /* Add a few more to cause the table to rehash */
        javaMap.putAll(generate());

    }

    private Map<Integer, Double> generate() {
        Map<Integer, Double> map = new HashMap<>();
        Random r = new Random();
        for (int i = 0; i < 2000; ++i)
            map.put(r.nextInt(), r.nextDouble());
        return map;
    }

    private OpenIntToDoubleHashMap createFromJavaMap() {
        OpenIntToDoubleHashMap map = new OpenIntToDoubleHashMap();
        for (Map.Entry<Integer, Double> mapEntry : javaMap.entrySet()) {
            map.put(mapEntry.getKey(), mapEntry.getValue());
        }
        return map;
    }

    @Test
    void testPutAndGetWith0ExpectedSize() {
        OpenIntToDoubleHashMap map = new OpenIntToDoubleHashMap(0);
        customAssertPutAndGet(map);
    }

    @Test
    void testPutAndGetWithExpectedSize() {
        OpenIntToDoubleHashMap map = new OpenIntToDoubleHashMap(500);
        customAssertPutAndGet(map);
    }

    @Test
    void testPutAndGet() {
        OpenIntToDoubleHashMap map = new OpenIntToDoubleHashMap();
        customAssertPutAndGet(map);
    }

    private void customAssertPutAndGet(OpenIntToDoubleHashMap map) {
        customAssertPutAndGet(map, 0, new HashSet<>());
    }

    private void customAssertPutAndGet(OpenIntToDoubleHashMap map, int mapSize,
                                       Set<Integer> keysInMap) {
        assertEquals(mapSize, map.size());
        for (Map.Entry<Integer, Double> mapEntry : javaMap.entrySet()) {
            map.put(mapEntry.getKey(), mapEntry.getValue());
            if (!keysInMap.contains(mapEntry.getKey()))
                ++mapSize;
            assertEquals(mapSize, map.size());
            assertTrue(Precision.equals(mapEntry.getValue(), map.get(mapEntry.getKey()), 1));
        }
    }

    @Test
    void testPutAbsentOnExisting() {
        OpenIntToDoubleHashMap map = createFromJavaMap();
        int size = javaMap.size();
        for (Map.Entry<Integer, Double> mapEntry : generateAbsent().entrySet()) {
            map.put(mapEntry.getKey(), mapEntry.getValue());
            assertEquals(++size, map.size());
            assertTrue(Precision.equals(mapEntry.getValue(), map.get(mapEntry.getKey()), 1));
        }
    }

    @Test
    void testPutOnExisting() {
        OpenIntToDoubleHashMap map = createFromJavaMap();
        for (Map.Entry<Integer, Double> mapEntry : javaMap.entrySet()) {
            map.put(mapEntry.getKey(), mapEntry.getValue());
            assertEquals(javaMap.size(), map.size());
            assertTrue(Precision.equals(mapEntry.getValue(), map.get(mapEntry.getKey()), 1));
        }
    }

    @Test
    void testGetAbsent() {
        Map<Integer, Double> generated = generateAbsent();
        OpenIntToDoubleHashMap map = createFromJavaMap();

        for (Map.Entry<Integer, Double> mapEntry : generated.entrySet())
            assertTrue(Double.isNaN(map.get(mapEntry.getKey())));
    }

    @Test
    void testGetFromEmpty() {
        OpenIntToDoubleHashMap map = new OpenIntToDoubleHashMap();
        assertTrue(Double.isNaN(map.get(5)));
        assertTrue(Double.isNaN(map.get(0)));
        assertTrue(Double.isNaN(map.get(50)));
    }

    @Test
    void testRemove() {
        OpenIntToDoubleHashMap map = createFromJavaMap();
        int mapSize = javaMap.size();
        assertEquals(mapSize, map.size());
        for (Map.Entry<Integer, Double> mapEntry : javaMap.entrySet()) {
            map.remove(mapEntry.getKey());
            assertEquals(--mapSize, map.size());
            assertTrue(Double.isNaN(map.get(mapEntry.getKey())));
        }

        /* Ensure that put and get still work correctly after removals */
        customAssertPutAndGet(map);
    }

    /* This time only remove some entries */
    @Test
    void testRemove2() {
        OpenIntToDoubleHashMap map = createFromJavaMap();
        int mapSize = javaMap.size();
        int count = 0;
        Set<Integer> keysInMap = new HashSet<>(javaMap.keySet());
        for (Map.Entry<Integer, Double> mapEntry : javaMap.entrySet()) {
            keysInMap.remove(mapEntry.getKey());
            map.remove(mapEntry.getKey());
            assertEquals(--mapSize, map.size());
            assertTrue(Double.isNaN(map.get(mapEntry.getKey())));
            if (count++ > 5)
                break;
        }

        /* Ensure that put and get still work correctly after removals */
        customAssertPutAndGet(map, mapSize, keysInMap);
    }

    @Test
    void testRemoveFromEmpty() {
        OpenIntToDoubleHashMap map = new OpenIntToDoubleHashMap();
        assertTrue(Double.isNaN(map.remove(50)));
    }

    @Test
    void testRemoveAbsent() {
        Map<Integer, Double> generated = generateAbsent();

        OpenIntToDoubleHashMap map = createFromJavaMap();
        int mapSize = map.size();

        for (Map.Entry<Integer, Double> mapEntry : generated.entrySet()) {
            map.remove(mapEntry.getKey());
            assertEquals(mapSize, map.size());
            assertTrue(Double.isNaN(map.get(mapEntry.getKey())));
        }
    }

    /**
     * Returns a map with at least 100 elements where each element is absent from javaMap.
     */
    private Map<Integer, Double> generateAbsent() {
        Map<Integer, Double> generated = new HashMap<>();
        do {
            generated.putAll(generate());
            for (Integer key : javaMap.keySet())
                generated.remove(key);
        } while (generated.size() < 100);
        return generated;
    }

    @Test
    void testCopy() {
        OpenIntToDoubleHashMap copy =
            new OpenIntToDoubleHashMap(createFromJavaMap());
        assertEquals(javaMap.size(), copy.size());

        for (Map.Entry<Integer, Double> mapEntry : javaMap.entrySet())
            assertTrue(Precision.equals(mapEntry.getValue(), copy.get(mapEntry.getKey()), 1));
    }

    @Test
    void testContainsKey() {
        OpenIntToDoubleHashMap map = createFromJavaMap();
        for (Map.Entry<Integer, Double> mapEntry : javaMap.entrySet()) {
            assertTrue(map.containsKey(mapEntry.getKey()));
        }
        for (Map.Entry<Integer, Double> mapEntry : generateAbsent().entrySet()) {
            assertFalse(map.containsKey(mapEntry.getKey()));
        }
        for (Map.Entry<Integer, Double> mapEntry : javaMap.entrySet()) {
            int key = mapEntry.getKey();
            assertTrue(map.containsKey(key));
            map.remove(key);
            assertFalse(map.containsKey(key));
        }
    }

    @Test
    void testIterator() {
        OpenIntToDoubleHashMap map = createFromJavaMap();
        OpenIntToDoubleHashMap.Iterator iterator = map.iterator();
        for (int i = 0; i < map.size(); ++i) {
            assertTrue(iterator.hasNext());
            iterator.advance();
            int key = iterator.key();
            assertTrue(map.containsKey(key));
            assertEquals(javaMap.get(key), map.get(key), 0);
            assertEquals(javaMap.get(key), iterator.value(), 0);
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
        OpenIntToDoubleHashMap map1 = new OpenIntToDoubleHashMap();
        map1.put(2,   2.5);
        map1.put(17, -0.5);
        map1.put(16,  0.0);
        assertEquals(map1, map1);
        OpenIntToDoubleHashMap map2 = new OpenIntToDoubleHashMap();
        map2.put(17, -0.5);
        map2.put(2,   2.5);
        map2.put(16,  0.0);
        assertEquals(map1, map2);
        map2.put(16,  0.25);
        assertNotEquals(map1, map2);
        map2.put(16,  0.0);
        assertEquals(map1, map2);
        assertNotEquals("", map1);
        assertNotEquals(null, map1);
    }

    @Test
    void testHashcode() {
        OpenIntToDoubleHashMap map = new OpenIntToDoubleHashMap();
        map.put(2,   2.5);
        map.put(17, -0.5);
        map.put(16,  0.0);
        assertEquals(-686630537, map.hashCode());
    }

    @Test
    void testConcurrentModification() {
        OpenIntToDoubleHashMap map = createFromJavaMap();
        OpenIntToDoubleHashMap.Iterator iterator = map.iterator();
        map.put(3, 3);
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
        OpenIntToDoubleHashMap map = new OpenIntToDoubleHashMap();
        int key1 = -1996012590;
        double value1 = 1.0;
        map.put(key1, value1);
        int key2 = 835099822;
        map.put(key2, value1);
        int key3 = 1008859686;
        map.put(key3, value1);
        assertTrue(Precision.equals(value1, map.get(key3), 1));
        assertEquals(3, map.size());

        map.remove(key2);
        double value2 = 2.0;
        map.put(key3, value2);
        assertTrue(Precision.equals(value2, map.get(key3), 1));
        assertEquals(2, map.size());
    }

    /**
     * Similar to testPutKeysWithCollisions() but exercises the codepaths in a slightly
     * different manner.
     */
    @Test
    void testPutKeysWithCollision2() {
        OpenIntToDoubleHashMap map = new OpenIntToDoubleHashMap();
        int key1 = 837989881;
        double value1 = 1.0;
        map.put(key1, value1);
        int key2 = 476463321;
        map.put(key2, value1);
        assertEquals(2, map.size());
        assertTrue(Precision.equals(value1, map.get(key2), 1));

        map.remove(key1);
        double value2 = 2.0;
        map.put(key2, value2);
        assertEquals(1, map.size());
        assertTrue(Precision.equals(value2, map.get(key2), 1));
    }

}
