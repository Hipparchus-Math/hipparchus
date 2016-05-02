/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements. See the NOTICE file distributed with this
 * work for additional information regarding copyright ownership. The ASF
 * licenses this file to You under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0 Unless required by applicable law
 * or agreed to in writing, software distributed under the License is
 * distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied. See the License for the specific language
 * governing permissions and limitations under the License.
 */
package org.hipparchus.util;

import org.hipparchus.exception.MathIllegalStateException;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test for {@link Incrementor}.
 */
public class IncrementorTest {
    @Test
    public void testConstructor1() {
        final Incrementor i = new Incrementor();
        Assert.assertEquals(Integer.MAX_VALUE, i.getMaximalCount());
        Assert.assertEquals(0, i.getCount());
    }

    @Test
    public void testConstructor2() {
        final Incrementor i = new Incrementor(10);
        Assert.assertEquals(10, i.getMaximalCount());
        Assert.assertEquals(0, i.getCount());
    }

    @Test
    public void testCanIncrement1() {
        final Incrementor i = new Incrementor(3);
        Assert.assertTrue(i.canIncrement());
        i.increment();
        Assert.assertTrue(i.canIncrement());
        i.increment();
        Assert.assertTrue(i.canIncrement());
        i.increment();
        Assert.assertFalse(i.canIncrement());
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
            Assert.fail("MathIllegalStateException expected");
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

        Assert.assertEquals(2, i.getCount());
    }

    @Test
    public void testAccessor() {
        final Incrementor i = new Incrementor(10);

        Assert.assertEquals(10, i.getMaximalCount());
        Assert.assertEquals(0, i.getCount());
    }

    @Test
    public void testBelowMaxCount() {
        final Incrementor i = new Incrementor(3);

        i.increment();
        i.increment();
        i.increment();

        Assert.assertEquals(3, i.getCount());
    }

    @Test(expected=MathIllegalStateException.class)
    public void testAboveMaxCount() {
        final Incrementor i = new Incrementor(3);

        i.increment();
        i.increment();
        i.increment();
        i.increment();
    }

    @Test(expected=IllegalStateException.class)
    public void testAlternateException() {
        final Incrementor.MaxCountExceededCallback cb =
            (int max) -> {
                throw new IllegalStateException();
            };

        final Incrementor i = new Incrementor(0, cb);
        i.increment();
    }

    @Test
    public void testReset() {
        final Incrementor i = new Incrementor(3);

        i.increment();
        i.increment();
        i.increment();
        Assert.assertEquals(3, i.getCount());
        i.reset();
        Assert.assertEquals(0, i.getCount());
    }

    @Test
    public void testBulkIncrement() {
        final Incrementor i = new Incrementor(3);

        i.increment(2);
        Assert.assertEquals(2, i.getCount());
        i.increment(1);
        Assert.assertEquals(3, i.getCount());
    }

    @Test(expected=MathIllegalStateException.class)
    public void testBulkIncrementExceeded() {
        final Incrementor i = new Incrementor(3);

        i.increment(2);
        Assert.assertEquals(2, i.getCount());
        i.increment(2);
    }

    @Test
    public void testWithMaximalValue()
    {
        final Incrementor i = new Incrementor(3);

        Assert.assertEquals(3, i.getMaximalCount());

        Incrementor i2 = i.withMaximalCount(10);

        Assert.assertNotEquals(i, i2);
        Assert.assertEquals(3, i.getMaximalCount());
        Assert.assertEquals(10, i2.getMaximalCount());
    }
}
