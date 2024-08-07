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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.special.Gamma;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Test cases for the {@link CombinatoricsUtils.FactorialLog} class.
 */
class FactorialLogTest {

    @Test
    void testPrecondition1() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            CombinatoricsUtils.FactorialLog.create().withCache(-1);
        });
    }

    @Test
    void testNonPositiveArgument() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final CombinatoricsUtils.FactorialLog f = CombinatoricsUtils.FactorialLog.create();
            f.value(-1);
        });
    }

    @Test
    void testDelegation() {
        final CombinatoricsUtils.FactorialLog f = CombinatoricsUtils.FactorialLog.create();

        // Starting at 21 because for smaller arguments, there is no delegation to the
        // "Gamma" class.
        for (int i = 21; i < 10000; i++) {
            final double expected = Gamma.logGamma(i + 1);
            assertEquals(expected, f.value(i), 0d, i + "! ");
        }
    }

    @Test
    void testCompareDirectWithoutCache() {
        // This test shows that delegating to the "Gamma" class will also lead to a
        // less accurate result.

        final int max = 100;
        final CombinatoricsUtils.FactorialLog f = CombinatoricsUtils.FactorialLog.create();

        for (int i = 0; i < max; i++) {
            final double expected = factorialLog(i);
            assertEquals(expected, f.value(i), 2 * Math.ulp(expected), i + "! ");
        }
    }

    @Test
    void testCompareDirectWithCache() {
        final int max = 1000;
        final CombinatoricsUtils.FactorialLog f = CombinatoricsUtils.FactorialLog.create().withCache(max);

        for (int i = 0; i < max; i++) {
            final double expected = factorialLog(i);
            assertEquals(expected, f.value(i), 0d, i + "! ");
        }
    }

    @Test
    void testCacheIncrease() {
        final int max = 100;
        final CombinatoricsUtils.FactorialLog f1 = CombinatoricsUtils.FactorialLog.create().withCache(max);
        final CombinatoricsUtils.FactorialLog f2 = f1.withCache(2 * max);

        final int val = max + max / 2;
        final double expected = factorialLog(val);
        assertEquals(expected, f2.value(val), 0d);
    }

    @Test
    void testCacheDecrease() {
        final int max = 100;
        final CombinatoricsUtils.FactorialLog f1 = CombinatoricsUtils.FactorialLog.create().withCache(max);
        final CombinatoricsUtils.FactorialLog f2 = f1.withCache(max / 2);

        final int val = max / 4;
        final double expected = factorialLog(val);
        assertEquals(expected, f2.value(val), 0d);
    }

    // Direct implementation.
    private double factorialLog(final int n) {
        double logSum = 0;
        for (int i = 2; i <= n; i++) {
            logSum += FastMath.log(i);
        }
        return logSum;
    }
}
