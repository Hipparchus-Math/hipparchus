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
package org.hipparchus.stat.descriptive.summary;

import org.hipparchus.stat.descriptive.StorelessUnivariateStatistic;
import org.hipparchus.stat.descriptive.StorelessUnivariateStatisticAbstractTest;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test cases for the {@link SumOfSquares} class.
 */
public class SumSqTest extends StorelessUnivariateStatisticAbstractTest {

    @Override
    public SumOfSquares getUnivariateStatistic() {
        return new SumOfSquares();
    }

    @Override
    public double expectedValue() {
        return this.sumSq;
    }

    @Test
    void testSpecialValues() {
        SumOfSquares sumSq = getUnivariateStatistic();
        assertEquals(0, sumSq.getResult(), 0);
        sumSq.increment(2d);
        assertEquals(4d, sumSq.getResult(), 0);
        sumSq.increment(Double.POSITIVE_INFINITY);
        assertEquals(Double.POSITIVE_INFINITY, sumSq.getResult(), 0);
        sumSq.increment(Double.NEGATIVE_INFINITY);
        assertEquals(Double.POSITIVE_INFINITY, sumSq.getResult(), 0);
        sumSq.increment(Double.NaN);
        assertTrue(Double.isNaN(sumSq.getResult()));
        sumSq.increment(1);
        assertTrue(Double.isNaN(sumSq.getResult()));
    }

    @Override
    protected void checkClearValue(StorelessUnivariateStatistic statistic) {
        assertEquals(0, statistic.getResult(), 0);
    }

}
