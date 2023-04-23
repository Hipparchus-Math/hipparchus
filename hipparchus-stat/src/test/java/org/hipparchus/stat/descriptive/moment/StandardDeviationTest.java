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

import org.hipparchus.stat.StatUtils;
import org.hipparchus.stat.descriptive.StorelessUnivariateStatisticAbstractTest;
import org.hipparchus.util.FastMath;
import org.junit.Test;

/**
 * Test cases for the {@link StandardDeviation} class.
 */
public class StandardDeviationTest extends StorelessUnivariateStatisticAbstractTest{

    @Override
    public StandardDeviation getUnivariateStatistic() {
        return new StandardDeviation();
    }

    @Override
    public double expectedValue() {
        return this.std;
    }

    /**
     * Make sure Double.NaN is returned iff n = 0
     */
    @Test
    public void testNaN() {
        StandardDeviation std = getUnivariateStatistic();
        assertTrue(Double.isNaN(std.getResult()));
        std.increment(1d);
        assertEquals(0d, std.getResult(), 0);
    }

    /**
     * Test population version of variance
     */
    @Test
    public void testPopulation() {
        double[] values = { -1.0d, 3.1d, 4.0d, -2.1d, 22d, 11.7d, 3d, 14d };
        double sigma = populationStandardDeviation(values);
        SecondMoment m = new SecondMoment();
        m.incrementAll(values);  // side effect is to add values

        StandardDeviation s1 = getUnivariateStatistic();
        s1 = s1.withBiasCorrection(false);
        assertEquals(sigma, s1.evaluate(values), 1E-14);
        s1.incrementAll(values);
        assertEquals(sigma, s1.getResult(), 1E-14);
        s1 = new StandardDeviation(false, m);
        assertEquals(sigma, s1.getResult(), 1E-14);
        s1 = new StandardDeviation(false);
        assertEquals(sigma, s1.evaluate(values), 1E-14);
        s1.incrementAll(values);
        assertEquals(sigma, s1.getResult(), 1E-14);
    }

    /**
     * Definitional formula for population standard deviation
     */
    protected double populationStandardDeviation(double[] v) {
        double mean = StatUtils.mean(v);
        double sum = 0;
        for (double val : v) {
            sum += (val - mean) * (val - mean);
        }
        return FastMath.sqrt(sum / v.length);
    }

}
