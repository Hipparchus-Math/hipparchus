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
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Test for {@link Incrementor}.
 */
public class IncrementorTest {
    @Test
    public void testConstructor1() {
        final Incrementor i = new Incrementor();
        Assertions.assertEquals(Integer.MAX_VALUE, i.getMaximalCount());
        Assertions.assertEquals(0, i.getCount());
    }

    @Test
    public void testConstructor2() {
        final Incrementor i = new Incrementor(10);
        Assertions.assertEquals(10, i.getMaximalCount());
        Assertions.assertEquals(0, i.getCount());
    }

    @Test
    public void testCanIncrement1() {
        final Incrementor i = new Incrementor(3);
        Assertions.assertTrue(i.canIncrement());
        i.increment();
        Assertions.assertTrue(i.canIncrement());
        i.increment();
        Assertions.assertTrue(i.canIncrement());
        i.increment();
        Assertions.assertFalse(i.canIncrement());
    }

    @Test
    public void testCanIncrement2() {
        final Incrementor i = new Incrementor(3);
        while (i.canIncrement()) {
            i.increment();
        }

        // Must keep try/catch because the exception must be generated here,
        // and not in the previous loop.
        try {
            i.increment();
            Assertions.fail("MathIllegalStateException expected");
        } catch (MathIllegalStateException e) {
            // Expected.
        }
    }

    @Test
    public void testBulkCanIncrement() {
        final Incrementor i = new Incrementor(3);
        while (i.canIncrement(2)) {
            i.increment(2);
        }

        Assertions.assertEquals(2, i.getCount());
    }

    @Test
    public void testAccessor() {
        final Incrementor i = new Incrementor(10);

        Assertions.assertEquals(10, i.getMaximalCount());
        Assertions.assertEquals(0, i.getCount());
    }

    @Test
    public void testBelowMaxCount() {
        final Incrementor i = new Incrementor(3);

        i.increment();
        i.increment();
        i.increment();

        Assertions.assertEquals(3, i.getCount());
    }

    @Test
    public void testAboveMaxCount() {
        assertThrows(MathIllegalStateException.class, () -> {
            final Incrementor i = new Incrementor(3);

            i.increment();
            i.increment();
            i.increment();
            i.increment();
        });
    }

    @Test
    public void testAlternateException() {
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
    public void testReset() {
        final Incrementor i = new Incrementor(3);

        i.increment();
        i.increment();
        i.increment();
        Assertions.assertEquals(3, i.getCount());
        i.reset();
        Assertions.assertEquals(0, i.getCount());
    }

    @Test
    public void testBulkIncrement() {
        final Incrementor i = new Incrementor(3);

        i.increment(2);
        Assertions.assertEquals(2, i.getCount());
        i.increment(1);
        Assertions.assertEquals(3, i.getCount());
    }

    @Test
    public void testBulkIncrementExceeded() {
        assertThrows(MathIllegalStateException.class, () -> {
            final Incrementor i = new Incrementor(3);

            i.increment(2);
            Assertions.assertEquals(2, i.getCount());
            i.increment(2);
        });
    }

    @Test
    public void testWithMaximalValue()
    {
        final Incrementor i = new Incrementor(3);

        Assertions.assertEquals(3, i.getMaximalCount());

        Incrementor i2 = i.withMaximalCount(10);

        Assertions.assertNotEquals(i, i2);
        Assertions.assertEquals(3, i.getMaximalCount());
        Assertions.assertEquals(10, i2.getMaximalCount());
    }

    @Test
    public void testMaxInt() {
       final Incrementor i = new Incrementor().withCount(Integer.MAX_VALUE - 2);
        i.increment();
        Assertions.assertEquals(Integer.MAX_VALUE - 1, i.getCount());
        i.increment();
        Assertions.assertEquals(Integer.MAX_VALUE, i.getCount());
        Assertions.assertFalse(i.canIncrement());
        try {
            i.increment();
            Assertions.fail("an exception should have been throwns");
        } catch (MathIllegalStateException mise) {
            Assertions.assertEquals(LocalizedCoreFormats.MAX_COUNT_EXCEEDED, mise.getSpecifier());
        }
    }

}
