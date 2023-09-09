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
package org.hipparchus.stat.descriptive.moment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.hipparchus.stat.descriptive.StorelessUnivariateStatisticAbstractTest;
import org.junit.Test;

/**
 * Test cases for the {@link Mean} class.
 */
public class MeanTest extends StorelessUnivariateStatisticAbstractTest{

    @Override
    public Mean getUnivariateStatistic() {
        return new Mean();
    }

    @Override
    public double expectedValue() {
        return this.mean;
    }

    /** Expected value for the testArray defined in UnivariateStatisticAbstractTest */
    public double expectedWeightedValue() {
        return this.weightedMean;
    }

    @Test
    public void testSmallSamples() {
        Mean mean = getUnivariateStatistic();
        assertTrue(Double.isNaN(mean.getResult()));
        mean.increment(1d);
        assertEquals(1d, mean.getResult(), 0);
    }

    @Test
    public void testWeightedMean() {
        Mean mean = getUnivariateStatistic();
        assertEquals(expectedWeightedValue(),
                     mean.evaluate(testArray, testWeightsArray, 0, testArray.length),
                     getTolerance());
        assertEquals(expectedValue(),
                     mean.evaluate(testArray, identicalWeightsArray, 0, testArray.length),
                     getTolerance());
    }

}
