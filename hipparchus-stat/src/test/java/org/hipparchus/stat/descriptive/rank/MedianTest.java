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

import org.hipparchus.stat.descriptive.UnivariateStatistic;
import org.hipparchus.stat.descriptive.UnivariateStatisticAbstractTest;
import org.hipparchus.stat.descriptive.rank.Percentile.EstimationType;
import org.hipparchus.stat.ranking.NaNStrategy;
import org.junit.jupiter.api.Test;

import static org.hipparchus.stat.descriptive.rank.Percentile.EstimationType.LEGACY;
import static org.hipparchus.stat.descriptive.rank.Percentile.EstimationType.R_1;
import static org.hipparchus.stat.descriptive.rank.Percentile.EstimationType.R_2;
import static org.hipparchus.stat.descriptive.rank.Percentile.EstimationType.R_3;
import static org.hipparchus.stat.descriptive.rank.Percentile.EstimationType.R_4;
import static org.hipparchus.stat.descriptive.rank.Percentile.EstimationType.R_5;
import static org.hipparchus.stat.descriptive.rank.Percentile.EstimationType.R_6;
import static org.hipparchus.stat.descriptive.rank.Percentile.EstimationType.R_7;
import static org.hipparchus.stat.descriptive.rank.Percentile.EstimationType.R_8;
import static org.hipparchus.stat.descriptive.rank.Percentile.EstimationType.R_9;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test cases for the {@link Median} class.
 */
public class MedianTest extends UnivariateStatisticAbstractTest{

    /**
     * {@link Percentile.EstimationType type} to be used while calling
     * {@link #getUnivariateStatistic()}
     */
    private final EstimationType estimationType = LEGACY;

    @Override
    public Median getUnivariateStatistic() {
        return new Median().withEstimationType(estimationType);
    }

    private Median getTestMedian(EstimationType type) {
        NaNStrategy strategy = (type == LEGACY) ? NaNStrategy.FIXED : NaNStrategy.REMOVED;
        return new Median().withEstimationType(type).withNaNStrategy(strategy);
    }

    @Override
    public double expectedValue() {
        return this.median;
    }

    @Test
    void testAllTechniquesSingleton() {
        double[] singletonArray = new double[] { 1d };
        for (EstimationType e : EstimationType.values()) {
            UnivariateStatistic percentile = getTestMedian(e);
            assertEquals(1d, percentile.evaluate(singletonArray), 0);
            assertEquals(1d, percentile.evaluate(singletonArray, 0, 1), 0);
            assertTrue(Double.isNaN(percentile.evaluate(singletonArray, 0, 0)));
        }
    }

    @Test
    void testAllTechniquesMedian() {
        double[] d = new double[] { 1, 3, 2, 4 };
        testAssertMappedValues(d, new Object[][] { { LEGACY, 2.5d },
            { R_1, 2d }, { R_2, 2.5d }, { R_3, 2d }, { R_4, 2d }, { R_5, 2.5 },
            { R_6, 2.5 },{ R_7, 2.5 },{ R_8, 2.5 }, { R_9 , 2.5 } },  1.0e-05);
    }

    /**
     * Simple test assertion utility method
     *
     * @param d input data
     * @param map of expected result against a {@link EstimationType}
     * @param tolerance the tolerance of difference allowed
     */
    protected void testAssertMappedValues(double[] d, Object[][] map, Double tolerance) {
        for (Object[] o : map) {
            EstimationType e = (EstimationType) o[0];
            double expected = (Double) o[1];
            double result = getTestMedian(e).evaluate(d);
            assertEquals(expected, result, tolerance, "expected[" + e + "] = " + expected +
                         " but was = " + result);
        }
    }
}
