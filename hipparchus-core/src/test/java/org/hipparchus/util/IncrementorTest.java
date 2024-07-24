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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test for {@link Incrementor}.
 */
class IncrementorTest {
    @Test
    void testConstructor1() {
        final Incrementor i = new Incrementor();
        assertEquals(Integer.MAX_VALUE, i.getMaximalCount());
        assertEquals(0, i.getCount());
    }

    @Test
    void testConstructor2() {
        final Incrementor i = new Incrementor(10);
        assertEquals(10, i.getMaximalCount());
        assertEquals(0, i.getCount());
    }

    @Test
    void testCanIncrement1() {
        final Incrementor i = new Incrementor(3);
        assertTrue(i.canIncrement());
        i.increment();
        assertTrue(i.canIncrement());
        i.increment();
        assertTrue(i.canIncrement());
        i.increment();
        assertFalse(i.canIncrement());
    }

    @Test
    void testCanIncrement2() {
        final Incrementor i = new Incrementor(3);
        while (i.canIncrement()) {
            i.increment();
        }

        // Must keep try/catch because the exception must be generated here,
        // and not in the previous loop.
        try {
            i.increment();
            fail("MathIllegalStateException expected");
        } catch (MathIllegalStateException e) {
            // Expected.
        }
    }

    @Test
    void testBulkCanIncrement() {
        final Incrementor i = new Incrementor(3);
        while (i.canIncrement(2)) {
            i.increment(2);
        }

        assertEquals(2, i.getCount());
    }

    @Test
    void testAccessor() {
        final Incrementor i = new Incrementor(10);

        assertEquals(10, i.getMaximalCount());
        assertEquals(0, i.getCount());
    }

    @Test
    void testBelowMaxCount() {
        final Incrementor i = new Incrementor(3);

        i.increment();
        i.increment();
        i.increment();

        assertEquals(3, i.getCount());
    }

    @Test
    void testAboveMaxCount() {
        assertThrows(MathIllegalStateException.class, () -> {
            final Incrementor i = new Incrementor(3);

            i.increment();
            i.increment();
            i.increment();
            i.increment();
        });
    }

    @Test
    void testAlternateException() {
        assertThrows(IllegalStateException.class, () -> {
            final Incrementor.MaxCountExceededCallback cb =
                (int max) -> {
                    throw new IllegalStateException();
                };

            final Incrementor i = new Incrementor(0, cb);
            i.increment();
        });
    }

    @Test
    void testReset() {
        final Incrementor i = new Incrementor(3);

        i.increment();
        i.increment();
        i.increment();
        assertEquals(3, i.getCount());
        i.reset();
        assertEquals(0, i.getCount());
    }

    @Test
    void testBulkIncrement() {
        final Incrementor i = new Incrementor(3);

        i.increment(2);
        assertEquals(2, i.getCount());
        i.increment(1);
        assertEquals(3, i.getCount());
    }

    @Test
    void testBulkIncrementExceeded() {
        assertThrows(MathIllegalStateException.class, () -> {
            final Incrementor i = new Incrementor(3);

            i.increment(2);
            assertEquals(2, i.getCount());
            i.increment(2);
        });
    }

    @Test
    void testWithMaximalValue()
    {
        final Incrementor i = new Incrementor(3);

        assertEquals(3, i.getMaximalCount());

        Incrementor i2 = i.withMaximalCount(10);

        assertNotEquals(i, i2);
        assertEquals(3, i.getMaximalCount());
        assertEquals(10, i2.getMaximalCount());
    }

    @Test
    void testMaxInt() {
       final Incrementor i = new Incrementor().withCount(Integer.MAX_VALUE - 2);
        i.increment();
        assertEquals(Integer.MAX_VALUE - 1, i.getCount());
        i.increment();
        assertEquals(Integer.MAX_VALUE, i.getCount());
        assertFalse(i.canIncrement());
        try {
            i.increment();
            fail("an exception should have been throwns");
        } catch (MathIllegalStateException mise) {
            assertEquals(LocalizedCoreFormats.MAX_COUNT_EXCEEDED, mise.getSpecifier());
        }
    }

}
