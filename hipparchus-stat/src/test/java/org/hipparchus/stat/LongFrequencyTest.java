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
package org.hipparchus.stat;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.hipparchus.UnitTestUtils;
import org.junit.Test;

/**
 * Test cases for the {@link Frequency} class.
 */
public final class LongFrequencyTest {
    private static final long ONE_LONG = 1L;
    private static final long TWO_LONG = 2L;
    private static final long THREE_LONG = 3L;
    private static final int  ONE = 1;
    private static final int  TWO = 2;
    private static final int  THREE = 3 ;

    private static final double TOLERANCE = 10E-15d;

    /** test freq counts */
    @Test
    public void testCounts() {
        LongFrequency f = new LongFrequency();

        assertEquals("total count", 0, f.getSumFreq());
        f.addValue(ONE_LONG);
        f.addValue(TWO_LONG);
        f.addValue(1L);
        f.addValue((long) ONE);
        assertEquals("one frequency count", 3, f.getCount(1L));
        assertEquals("two frequency count", 1, f.getCount(2L));
        assertEquals("three frequency count", 0, f.getCount(3L));
        assertEquals("total count", 4, f.getSumFreq());
        assertEquals("zero cumulative frequency", 0, f.getCumFreq(0L));
        assertEquals("one cumulative frequency", 3,  f.getCumFreq(1L));
        assertEquals("two cumulative frequency", 4,  f.getCumFreq(2L));
        assertEquals("five cumulative frequency", 4, f.getCumFreq(5L));

        f.clear();
        assertEquals("total count", 0, f.getSumFreq());
    }

    /** test pcts */
    @Test
    public void testPcts() {
        LongFrequency f = new LongFrequency();

        f.addValue(ONE_LONG);
        f.addValue(TWO_LONG);
        f.addValue((long) ONE);
        f.addValue((long) TWO);
        f.addValue(THREE_LONG);
        f.addValue(THREE_LONG);
        f.addValue(3L);
        f.addValue((long) THREE);
        assertEquals("one pct", 0.25, f.getPct(1L), TOLERANCE);
        assertEquals("two pct", 0.25, f.getPct(Long.valueOf(2)), TOLERANCE);
        assertEquals("three pct", 0.5, f.getPct(THREE_LONG), TOLERANCE);
        assertEquals("five pct", 0, f.getPct(5L), TOLERANCE);
        assertEquals("one cum pct", 0.25, f.getCumPct(1L), TOLERANCE);
        assertEquals("two cum pct", 0.50, f.getCumPct(Long.valueOf(2)), TOLERANCE);
        assertEquals("three cum pct", 1.0, f.getCumPct(THREE_LONG), TOLERANCE);
        assertEquals("five cum pct", 1.0, f.getCumPct(5L), TOLERANCE);
        assertEquals("zero cum pct", 0.0, f.getCumPct(0L), TOLERANCE);
    }

    /**
     * Tests toString()
     */
    @Test
    public void testToString() throws Exception {
        LongFrequency f = new LongFrequency();

        f.addValue(ONE_LONG);
        f.addValue(TWO_LONG);
        f.addValue((long) ONE);
        f.addValue((long) TWO);

        String s = f.toString();
        //System.out.println(s);
        assertNotNull(s);
        BufferedReader reader = new BufferedReader(new StringReader(s));
        String line = reader.readLine(); // header line
        assertNotNull(line);

        line = reader.readLine(); // one's or two's line
        assertNotNull(line);

        line = reader.readLine(); // one's or two's line
        assertNotNull(line);

        line = reader.readLine(); // no more elements
        assertNull(line);
    }

    @Test
    public void testLongValues() {
        LongFrequency f = new LongFrequency();

        Integer obj1 = null;
        obj1 = Integer.valueOf(1);
        Integer int1 = Integer.valueOf(1);
        f.addValue(obj1);
        f.addValue(int1);
        f.addValue(2);
        f.addValue(Long.valueOf(2).intValue());
        assertEquals("Integer 1 count", 2, f.getCount(1));
        assertEquals("Integer 1 count", 2, f.getCount(Integer.valueOf(1)));
        assertEquals("Integer 1 count", 2, f.getCount(Long.valueOf(1).intValue()));
        assertEquals("Integer 1 cumPct", 0.5, f.getCumPct(1), TOLERANCE);
        assertEquals("Integer 1 cumPct", 0.5, f.getCumPct(Long.valueOf(1).intValue()), TOLERANCE);
        assertEquals("Integer 1 cumPct", 0.5, f.getCumPct(Integer.valueOf(1)), TOLERANCE);

        f.incrementValue(ONE, -2);
        f.incrementValue(THREE, 5);

        assertEquals("Integer 1 count", 0, f.getCount(1));
        assertEquals("Integer 3 count", 5, f.getCount(3));

        Iterator<?> it = f.valuesIterator();
        while (it.hasNext()) {
            assertTrue(it.next() instanceof Long);
        }
    }

    @Test
    public void testSerial() {
        LongFrequency f = new LongFrequency();

        f.addValue(ONE_LONG);
        f.addValue(TWO_LONG);
        f.addValue((long) ONE);
        f.addValue((long) TWO);
        assertEquals(f, UnitTestUtils.serializeAndRecover(f));
    }

    @Test
    public void testGetUniqueCount() {
        LongFrequency f = new LongFrequency();

        assertEquals(0, f.getUniqueCount());
        f.addValue(ONE_LONG);
        assertEquals(1, f.getUniqueCount());
        f.addValue(ONE_LONG);
        assertEquals(1, f.getUniqueCount());
        f.addValue((long) TWO);
        assertEquals(2, f.getUniqueCount());
    }

    @Test
    public void testIncrement() {
        LongFrequency f = new LongFrequency();

        assertEquals(0, f.getUniqueCount());
        f.incrementValue(ONE_LONG, 1);
        assertEquals(1, f.getCount(ONE_LONG));

        f.incrementValue(ONE_LONG, 4);
        assertEquals(5, f.getCount(ONE_LONG));

        f.incrementValue(ONE_LONG, -5);
        assertEquals(0, f.getCount(ONE_LONG));

    }

    @Test
    public void testMerge() {
        LongFrequency f = new LongFrequency();

        assertEquals(0, f.getUniqueCount());
        f.addValue(ONE_LONG);
        f.addValue(TWO_LONG);
        f.addValue((long) ONE);
        f.addValue((long) TWO);

        assertEquals(2, f.getUniqueCount());
        assertEquals(2, f.getCount((long) ONE));
        assertEquals(2, f.getCount((long) TWO));

        LongFrequency g = new LongFrequency();
        g.addValue(ONE_LONG);
        g.addValue(THREE_LONG);
        g.addValue((long) THREE);

        assertEquals(2, g.getUniqueCount());
        assertEquals(1, g.getCount((long) ONE));
        assertEquals(2, g.getCount((long) THREE));

        f.merge(g);

        assertEquals(3, f.getUniqueCount());
        assertEquals(3, f.getCount((long) ONE));
        assertEquals(2, f.getCount((long) TWO));
        assertEquals(2, f.getCount((long) THREE));
    }

    @Test
    public void testMergeCollection() {
        LongFrequency f = new LongFrequency();

        assertEquals(0, f.getUniqueCount());
        f.addValue(ONE_LONG);

        assertEquals(1, f.getUniqueCount());
        assertEquals(1, f.getCount((long) ONE));
        assertEquals(0, f.getCount((long) TWO));

        LongFrequency g = new LongFrequency();
        g.addValue(TWO_LONG);

        LongFrequency h = new LongFrequency();
        h.addValue(THREE_LONG);

        List<LongFrequency> coll = new ArrayList<>();
        coll.add(g);
        coll.add(h);
        f.merge(coll);

        assertEquals(3, f.getUniqueCount());
        assertEquals(1, f.getCount((long) ONE));
        assertEquals(1, f.getCount((long) TWO));
        assertEquals(1, f.getCount((long) THREE));
    }

    @Test
    public void testMode() {
        LongFrequency f = new LongFrequency();

        List<Long> mode;
        mode = f.getMode();
        assertEquals(0, mode.size());

        f.addValue(3L);
        mode = f.getMode();
        assertEquals(1, mode.size());
        assertEquals(Long.valueOf(3L), mode.get(0));

        f.addValue(2L);
        mode = f.getMode();
        assertEquals(2, mode.size());
        assertEquals(Long.valueOf(2L), mode.get(0));
        assertEquals(Long.valueOf(3L),mode.get(1));

        f.addValue(2L);
        mode = f.getMode();
        assertEquals(1, mode.size());
        assertEquals(Long.valueOf(2L), mode.get(0));
        assertFalse(mode.contains(1L));
        assertTrue(mode.contains(2L));
    }

}
