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
 * Test cases for the {@link Sum} class.
 */
public class SumTest extends StorelessUnivariateStatisticAbstractTest {

    @Override
    public Sum getUnivariateStatistic() {
        return new Sum();
    }

    @Override
    public double expectedValue() {
        return this.sum;
    }

    /** Expected value for the testArray defined in UnivariateStatisticAbstractTest */
    public double expectedWeightedValue() {
        return this.weightedSum;
    }

    @Test
    void testSpecialValues() {
        Sum sum = getUnivariateStatistic();
        assertEquals(0, sum.getResult(), 0);
        sum.increment(1);
        assertEquals(1, sum.getResult(), 0);
        sum.increment(Double.POSITIVE_INFINITY);
        assertEquals(Double.POSITIVE_INFINITY, sum.getResult(), 0);
        sum.increment(Double.NEGATIVE_INFINITY);
        assertTrue(Double.isNaN(sum.getResult()));
        sum.increment(1);
        assertTrue(Double.isNaN(sum.getResult()));
    }

    @Test
    void testWeightedSum() {
        Sum sum = new Sum();
        assertEquals(expectedWeightedValue(),
                     sum.evaluate(testArray, testWeightsArray, 0, testArray.length), getTolerance());
        assertEquals(expectedValue(),
                     sum.evaluate(testArray, unitWeightsArray, 0, testArray.length), getTolerance());
    }

    @Override
    protected void checkClearValue(StorelessUnivariateStatistic statistic) {
        assertEquals(0, statistic.getResult(), 0);
    }

}
