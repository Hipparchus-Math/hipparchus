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

import org.hipparchus.UnitTestUtils;
import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test cases for the {@link Frequency} class.
 */
final class FrequencyTest {
    private static final long ONE_LONG = 1L;
    private static final long TWO_LONG = 2L;
    private static final long THREE_LONG = 3L;
    private static final int  ONE = 1;
    private static final int  TWO = 2;
    private static final int  THREE = 3 ;
    private static final char CHAR_A = 'a';

    private static final double TOLERANCE = 10E-15d;

    /** test freq counts */
    @Test
    void testCounts() {
        Frequency<Long> f = new Frequency<>();

        assertEquals(0, f.getSumFreq(), "total count");
        f.addValue(ONE_LONG);
        f.addValue(TWO_LONG);
        f.addValue(1L);
        f.addValue((long) ONE);
        assertEquals(3, f.getCount(1L), "one frequency count");
        assertEquals(1, f.getCount(2L), "two frequency count");
        assertEquals(0, f.getCount(3L), "three frequency count");
        assertEquals(4, f.getSumFreq(), "total count");
        assertEquals(0, f.getCumFreq(0L), "zero cumulative frequency");
        assertEquals(3,  f.getCumFreq(1L),  "one cumulative frequency");
        assertEquals(4,  f.getCumFreq(2L),  "two cumulative frequency");
        assertEquals(4, f.getCumFreq(5L), "five cumulative frequency");

        f.clear();
        assertEquals(0, f.getSumFreq(), "total count");
    }

    @Test
    void testCountsString() {
        Frequency<String> f = new Frequency<>();

        f.addValue("one");
        f.addValue("One");
        f.addValue("oNe");
        f.addValue("Z");

        assertEquals(1, f.getCount("one"), "one cumulative frequency");
        assertEquals(0.5, f.getCumPct("Z"), TOLERANCE, "Z cumulative pct");
        assertEquals(1.0, f.getCumPct("z"), TOLERANCE, "z cumulative pct");
        assertEquals(0.25, f.getCumPct("Ot"), TOLERANCE, "Ot cumulative pct");
        f.clear();
    }

    @Test
    void testCountsInteger() {
        Frequency<Integer> f = new Frequency<>();

        f.addValue(1);
        f.addValue(Integer.valueOf(1));
        f.addValue(Long.valueOf(1).intValue());
        f.addValue(2);
        f.addValue(Integer.valueOf(-1));

        assertEquals(3, f.getCount(1), "1 count");
        assertEquals(3, f.getCount(Integer.valueOf(1)), "1 count");
        assertEquals(0.2, f.getCumPct(0), TOLERANCE, "0 cum pct");
        assertEquals(0.6, f.getPct(Integer.valueOf(1)), TOLERANCE, "1 pct");
        assertEquals(0, f.getCumPct(-2), TOLERANCE, "-2 cum pct");
        assertEquals(1, f.getCumPct(10), TOLERANCE, "10 cum pct");
    }

    @Test
    void testCountsComparator() {
        Frequency<String> f = new Frequency<>(String.CASE_INSENSITIVE_ORDER);

        f.addValue("one");
        f.addValue("One");
        f.addValue("oNe");
        f.addValue("Z");
        assertEquals(3 ,  f.getCount("one"),  "one count");
        assertEquals(1 ,  f.getCumPct("Z"), TOLERANCE, "Z cumulative pct -- case insensitive");
        assertEquals(1 ,  f.getCumPct("z"), TOLERANCE, "z cumulative pct -- case insensitive");
    }

    @Test
    void testCountsCharacter() {
        Frequency<Character> f = new Frequency<>();

        assertEquals(0L, f.getCount('a'));
        assertEquals(0L, f.getCumFreq('b'));
        UnitTestUtils.customAssertEquals(Double.NaN, f.getPct('a'), 0.0);
        UnitTestUtils.customAssertEquals(Double.NaN, f.getCumPct('b'), 0.0);
        f.addValue('a');
        f.addValue('b');
        f.addValue('c');
        f.addValue('d');
        assertEquals(1L, f.getCount('a'));
        assertEquals(2L, f.getCumFreq('b'));
        assertEquals(0.25, f.getPct('a'), 0.0);
        assertEquals(0.5, f.getCumPct('b'), 0.0);
        assertEquals(1.0, f.getCumPct('e'), 0.0);
    }

    /** test pcts */
    @Test
    void testPcts() {
        Frequency<Long> f = new Frequency<>();

        f.addValue(ONE_LONG);
        f.addValue(TWO_LONG);
        f.addValue((long) ONE);
        f.addValue((long) TWO);
        f.addValue(THREE_LONG);
        f.addValue(THREE_LONG);
        f.addValue(3L);
        f.addValue((long) THREE);
        assertEquals(0.25, f.getPct(1L), TOLERANCE, "one pct");
        assertEquals(0.25, f.getPct(Long.valueOf(2)), TOLERANCE, "two pct");
        assertEquals(0.5, f.getPct(THREE_LONG), TOLERANCE, "three pct");
        assertEquals(0, f.getPct(5L), TOLERANCE, "five pct");
        assertEquals(0.25, f.getCumPct(1L), TOLERANCE, "one cum pct");
        assertEquals(0.50, f.getCumPct(Long.valueOf(2)), TOLERANCE, "two cum pct");
        assertEquals(1.0, f.getCumPct(THREE_LONG), TOLERANCE, "three cum pct");
        assertEquals(1.0, f.getCumPct(5L), TOLERANCE, "five cum pct");
        assertEquals(0.0, f.getCumPct(0L), TOLERANCE, "zero cum pct");
    }

    /** test empty table */
    @Test
    void testEmptyTable() {
        final Frequency<Integer> freq = new Frequency<>();

        assertEquals(0, freq.getSumFreq(), "freq sum, empty table");
        assertEquals(0, freq.getCount(0), "count, empty table");
        assertEquals(0, freq.getCount(Integer.valueOf(0)), "count, empty table");
        assertEquals(0, freq.getCumFreq(0), "cum freq, empty table");
        assertEquals(0, freq.getCumFreq(2), "cum freq, empty table");

        assertTrue(Double.isNaN(freq.getPct(0)), "pct, empty table");
        assertTrue(Double.isNaN(freq.getPct(Integer.valueOf(0))), "pct, empty table");
        assertTrue(Double.isNaN(freq.getCumPct(0)), "cum pct, empty table");
        assertTrue(Double.isNaN(freq.getCumPct(Integer.valueOf(0))), "cum pct, empty table");
    }

    @Test
    void testCumPct() {
        Frequency<Integer> f = new Frequency<>();
        f.addValue(3);
        f.addValue(4);
        f.addValue(5);
        f.addValue(6);

        assertEquals(0.25d, f.getCumPct(3), TOLERANCE, "cum freq, single entry");
        assertEquals(1.0d,  f.getCumPct(6), TOLERANCE, "cum freq, single entry");

        assertEquals(0.0d,  f.getCumPct(1), TOLERANCE, "cum freq, single entry");
        assertEquals(1.0d,  f.getCumPct(10), TOLERANCE, "cum freq, single entry");
    }

    /**
     * Tests toString()
     */
    @Test
    void testToString() throws Exception {
        Frequency<Long> f = new Frequency<>();

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
    void testIntegerValues() {
        Frequency<Integer> f = new Frequency<>();

        Integer obj1 = null;
        obj1 = Integer.valueOf(1);
        Integer int1 = Integer.valueOf(1);
        f.addValue(obj1);
        f.addValue(int1);
        f.addValue(2);
        f.addValue(Long.valueOf(2).intValue());
        assertEquals(2, f.getCount(1), "Integer 1 count");
        assertEquals(2, f.getCount(Integer.valueOf(1)), "Integer 1 count");
        assertEquals(2, f.getCount(Long.valueOf(1).intValue()), "Integer 1 count");
        assertEquals(0.5, f.getCumPct(1), TOLERANCE, "Integer 1 cumPct");
        assertEquals(0.5, f.getCumPct(Long.valueOf(1).intValue()), TOLERANCE, "Integer 1 cumPct");
        assertEquals(0.5, f.getCumPct(Integer.valueOf(1)), TOLERANCE, "Integer 1 cumPct");

        f.incrementValue(ONE, -2);
        f.incrementValue(THREE, 5);

        assertEquals(0, f.getCount(1), "Integer 1 count");
        assertEquals(5, f.getCount(3), "Integer 3 count");

        Iterator<?> it = f.valuesIterator();
        while (it.hasNext()) {
            assertInstanceOf(Integer.class, it.next());
        }
    }

    @Test
    void testSerial() {
        Frequency<Long> f = new Frequency<>();

        f.addValue(ONE_LONG);
        f.addValue(TWO_LONG);
        f.addValue((long) ONE);
        f.addValue((long) TWO);
        assertEquals(f, UnitTestUtils.serializeAndRecover(f));
    }

    @Test
    void testGetUniqueCount() {
        Frequency<Long> f = new Frequency<>();

        assertEquals(0, f.getUniqueCount());
        f.addValue(ONE_LONG);
        assertEquals(1, f.getUniqueCount());
        f.addValue(ONE_LONG);
        assertEquals(1, f.getUniqueCount());
        f.addValue((long) TWO);
        assertEquals(2, f.getUniqueCount());
    }

    @Test
    void testIncrement() {
        Frequency<Long> f = new Frequency<>();

        assertEquals(0, f.getUniqueCount());
        f.incrementValue(ONE_LONG, 1);
        assertEquals(1, f.getCount(ONE_LONG));

        f.incrementValue(ONE_LONG, 4);
        assertEquals(5, f.getCount(ONE_LONG));

        f.incrementValue(ONE_LONG, -5);
        assertEquals(0, f.getCount(ONE_LONG));

        Frequency<Character> f2 = new Frequency<>();
        f2.incrementValue(CHAR_A, 2);

        assertEquals(2, f2.getCount(CHAR_A));

        f2.incrementValue(CHAR_A, 3);
        assertEquals(5, f2.getCount(CHAR_A));
    }

    @Test
    void testMerge() {
        Frequency<Long> f = new Frequency<>();

        assertEquals(0, f.getUniqueCount());
        f.addValue(ONE_LONG);
        f.addValue(TWO_LONG);
        f.addValue((long) ONE);
        f.addValue((long) TWO);

        assertEquals(2, f.getUniqueCount());
        assertEquals(2, f.getCount((long) ONE));
        assertEquals(2, f.getCount((long) TWO));

        Frequency<Long> g = new Frequency<>();
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
    void testMergeCollection() {
        Frequency<Long> f = new Frequency<>();

        assertEquals(0, f.getUniqueCount());
        f.addValue(ONE_LONG);

        assertEquals(1, f.getUniqueCount());
        assertEquals(1, f.getCount((long) ONE));
        assertEquals(0, f.getCount((long) TWO));

        Frequency<Long> g = new Frequency<>();
        g.addValue(TWO_LONG);

        Frequency<Long> h = new Frequency<>();
        h.addValue(THREE_LONG);

        List<Frequency<Long>> coll = new ArrayList<>();
        coll.add(g);
        coll.add(h);
        f.merge(coll);

        assertEquals(3, f.getUniqueCount());
        assertEquals(1, f.getCount((long) ONE));
        assertEquals(1, f.getCount((long) TWO));
        assertEquals(1, f.getCount((long) THREE));
    }

    @Test
    void testMode() {
        Frequency<String> f = new Frequency<>();

        List<String> mode;
        mode = f.getMode();
        assertEquals(0, mode.size());

        f.addValue("3");
        mode = f.getMode();
        assertEquals(1, mode.size());
        assertEquals("3", mode.get(0));

        f.addValue("2");
        mode = f.getMode();
        assertEquals(2, mode.size());
        assertEquals("2", mode.get(0));
        assertEquals("3",mode.get(1));

        f.addValue("2");
        mode = f.getMode();
        assertEquals(1, mode.size());
        assertEquals("2", mode.get(0));
        assertFalse(mode.contains("1"));
        assertTrue(mode.contains("2"));
    }

    @Test
    void testModeDoubleNan() {
        Frequency<Double> f = new Frequency<>();

        List<Double> mode;
        f.addValue(Double.valueOf(Double.NaN));
        f.addValue(Double.valueOf(Double.NaN));
        f.addValue(Double.valueOf(Double.NaN));
        f.addValue(Double.valueOf(Double.NEGATIVE_INFINITY));
        f.addValue(Double.valueOf(Double.POSITIVE_INFINITY));
        f.addValue(Double.valueOf(Double.NEGATIVE_INFINITY));
        f.addValue(Double.valueOf(Double.POSITIVE_INFINITY));
        f.addValue(Double.valueOf(Double.NEGATIVE_INFINITY));
        f.addValue(Double.valueOf(Double.POSITIVE_INFINITY));
        mode = f.getMode();
        assertEquals(3, mode.size());
        assertEquals(Double.valueOf(Double.NEGATIVE_INFINITY), mode.get(0));
        assertEquals(Double.valueOf(Double.POSITIVE_INFINITY), mode.get(1));
        assertEquals(Double.valueOf(Double.NaN), mode.get(2));
    }

    @Test
    void testModeFloatNan() {
        Frequency<Float> f = new Frequency<>();

        List<Float> mode;
        f.addValue(Float.valueOf(Float.NaN));
        f.addValue(Float.valueOf(Float.NaN));
        f.addValue(Float.valueOf(Float.NaN));
        f.addValue(Float.valueOf(Float.NEGATIVE_INFINITY));
        f.addValue(Float.valueOf(Float.POSITIVE_INFINITY));
        f.addValue(Float.valueOf(Float.NEGATIVE_INFINITY));
        f.addValue(Float.valueOf(Float.POSITIVE_INFINITY));
        f.addValue(Float.valueOf(Float.NEGATIVE_INFINITY));
        f.addValue(Float.valueOf(Float.POSITIVE_INFINITY));
        mode = f.getMode();
        assertEquals(3, mode.size());
        assertEquals(Float.valueOf(Float.NEGATIVE_INFINITY), mode.get(0));
        assertEquals(Float.valueOf(Float.POSITIVE_INFINITY), mode.get(1));
        assertEquals(Float.valueOf(Float.NaN), mode.get(2));
    }

}
