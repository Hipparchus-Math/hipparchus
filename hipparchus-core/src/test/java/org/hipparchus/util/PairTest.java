/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements. See the NOTICE file distributed with this
 * work for additional information regarding copyright ownership. The ASF
 * licenses this file to You under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * https://www.apache.org/licenses/LICENSE-2.0 Unless required by applicable law
 * or agreed to in writing, software distributed under the License is
 * distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied. See the License for the specific language
 * governing permissions and limitations under the License.
 */
package org.hipparchus.util;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertSame;

/**
 * Test for {@link Pair}.
 */
class PairTest {

    @Test
    void testAccessor() {
        final Pair<Integer, Double> p = new Pair<>(Integer.valueOf(1), Double.valueOf(2));
        assertEquals(Integer.valueOf(1), p.getKey());
        assertEquals(2, p.getValue().doubleValue(), Math.ulp(1d));
    }

    @Test
    void testAccessor2() {
        final Pair<Integer, Double> p = new Pair<>(Integer.valueOf(1), Double.valueOf(2));

        // Check that both APIs refer to the same data.

        assertSame(p.getFirst(), p.getKey());
        assertSame(p.getSecond(), p.getValue());
    }

    @Test
    void testEquals() {
        Pair<Integer, Double> p1 = new Pair<>(null, null);
        assertNotEquals(null, p1);

        Pair<Integer, Double> p2 = new Pair<>(null, null);
        assertEquals(p1, p2);

        p1 = new Pair<>(Integer.valueOf(1), Double.valueOf(2));
        assertNotEquals(p1, p2);

        p2 = new Pair<>(Integer.valueOf(1), Double.valueOf(2));
        assertEquals(p1, p2);

        Pair<Integer, Float> p3 = new Pair<>(Integer.valueOf(1), Float.valueOf(2));
        assertNotEquals(p1, p3);
    }

    @Test
    void testHashCode() {
        final MyInteger m1 = new MyInteger(1);
        final MyInteger m2 = new MyInteger(1);

        final Pair<MyInteger, MyInteger> p1 = new Pair<>(m1, m1);
        final Pair<MyInteger, MyInteger> p2 = new Pair<>(m2, m2);
        // Same contents, same hash code.
        assertEquals(p1.hashCode(), p2.hashCode());

        // Different contents, different hash codes.
        m2.set(2);
        assertNotEquals(p1.hashCode(), p2.hashCode());
    }

    @Test
    void testToString() {
        assertEquals("[null, null]", new Pair<>(null, null).toString());
        assertEquals("[foo, 3]", new Pair<>("foo", 3).toString());
    }

    @Test
    void testCreate() {
        final Pair<String, Integer> p1 = Pair.create("foo", 3);
        assertNotNull(p1);
        final Pair<String, Integer> p2 = new Pair<>("foo", 3);
        assertEquals(p2, p1);
    }

    /**
     * A mutable integer.
     */
    private static class MyInteger {
        private int i;

        public MyInteger(int i) {
            this.i = i;
        }

        public void set(int i) {
            this.i = i;
        }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof MyInteger)) {
                return false;
            } else {
                return i == ((MyInteger) o).i;
            }
        }

        @Override
        public int hashCode() {
            return i;
        }
    }
}
