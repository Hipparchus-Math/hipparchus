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
package org.hipparchus.stat.descriptive.rank;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.hipparchus.stat.descriptive.StorelessUnivariateStatisticAbstractTest;
import org.junit.Test;

/**
 * Test cases for the {@link Max} class.
 */
public class MaxTest extends StorelessUnivariateStatisticAbstractTest {

    @Override
    public Max getUnivariateStatistic() {
        return new Max();
    }

    @Override
    public double expectedValue() {
        return this.max;
    }

    @Test
    public void testSpecialValues() {
        double[] testArray = {0d, Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY};
        Max max = getUnivariateStatistic();
        assertTrue(Double.isNaN(max.getResult()));
        max.increment(testArray[0]);
        assertEquals(0d, max.getResult(), 0);
        max.increment(testArray[1]);
        assertEquals(0d, max.getResult(), 0);
        max.increment(testArray[2]);
        assertEquals(0d, max.getResult(), 0);
        max.increment(testArray[3]);
        assertEquals(Double.POSITIVE_INFINITY, max.getResult(), 0);
        assertEquals(Double.POSITIVE_INFINITY, max.evaluate(testArray), 0);
    }

    @Test
    public void testNaNs() {
        Max max = getUnivariateStatistic();
        double nan = Double.NaN;
        assertEquals(3d, max.evaluate(new double[]{nan, 2d, 3d}), 0);
        assertEquals(3d, max.evaluate(new double[]{1d, nan, 3d}), 0);
        assertEquals(2d, max.evaluate(new double[]{1d, 2d, nan}), 0);
        assertTrue(Double.isNaN(max.evaluate(new double[]{nan, nan, nan})));
    }

}
