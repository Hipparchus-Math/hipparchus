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
package org.hipparchus.stat.descriptive;

import org.hipparchus.UnitTestUtils;
import org.junit.jupiter.api.Test;

import java.util.Locale;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

/**
 * Test cases for the {@link StatisticalSummaryValues} class.
 */
final class StatisticalSummaryValuesTest {

    @Test
    void testSerialization() {
        StatisticalSummaryValues u = new StatisticalSummaryValues(1, 2, 3, 4, 5, 6);
        UnitTestUtils.checkSerializedEquality(u);
        StatisticalSummaryValues t = (StatisticalSummaryValues) UnitTestUtils.serializeAndRecover(u);
        verifyEquality(u, t);
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testEqualsAndHashCode() {
        StatisticalSummaryValues u  = new StatisticalSummaryValues(1, 2, 3, 4, 5, 6);
        StatisticalSummaryValues t = null;
        assertEquals(u, u, "reflexive");
        assertNotEquals(u, t, "non-null compared to null");
        assertNotEquals(u, Double.valueOf(0), "wrong type");
        t = new StatisticalSummaryValues(1, 2, 3, 4, 5, 6);
        assertEquals(t, u, "instances with same data should be equal");
        assertEquals(u.hashCode(), t.hashCode(), "hash code");

        u = new StatisticalSummaryValues(Double.NaN, 2, 3, 4, 5, 6);
        t = new StatisticalSummaryValues(1, Double.NaN, 3, 4, 5, 6);
        assertFalse((u.equals(t) ||t.equals(u)),
                "instances based on different data should be different");
    }

    private void verifyEquality(StatisticalSummaryValues s, StatisticalSummaryValues u) {
        assertEquals(s.getN(),u.getN(),"N");
        UnitTestUtils.customAssertEquals("sum", s.getSum(), u.getSum(), 0);
        UnitTestUtils.customAssertEquals("var", s.getVariance(), u.getVariance(), 0);
        UnitTestUtils.customAssertEquals("std", s.getStandardDeviation(), u.getStandardDeviation(), 0);
        UnitTestUtils.customAssertEquals("mean", s.getMean(), u.getMean(), 0);
        UnitTestUtils.customAssertEquals("min", s.getMin(), u.getMin(), 0);
        UnitTestUtils.customAssertEquals("max", s.getMax(), u.getMax(), 0);
    }

    @Test
    void testToString() {
        StatisticalSummaryValues u  = new StatisticalSummaryValues(4.5, 16, 10, 5, 4, 45);
        Locale d = Locale.getDefault();
        Locale.setDefault(Locale.US);
        assertEquals("StatisticalSummaryValues:\n" +
                     "n: 10\n" +
                     "min: 4.0\n" +
                     "max: 5.0\n" +
                     "mean: 4.5\n" +
                     "std dev: 4.0\n" +
                     "variance: 16.0\n" +
                     "sum: 45.0\n",  u.toString());
        Locale.setDefault(d);
    }
}
