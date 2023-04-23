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
import org.hipparchus.util.MathArrays;
import org.junit.Test;

/**
 * Test cases for the {@link Variance} class.
 */
public class VarianceTest extends StorelessUnivariateStatisticAbstractTest{

    @Override
    public Variance getUnivariateStatistic() {
        return new Variance();
    }

    @Override
    public double expectedValue() {
        return this.var;
    }

    /** Expected value for  the testArray defined in UnivariateStatisticAbstractTest */
    public double expectedWeightedValue() {
        return this.weightedVar;
    }

    /**
     * Make sure Double.NaN is returned iff n = 0
     */
    @Test
    public void testNaN() {
        Variance var = getUnivariateStatistic();
        assertTrue(Double.isNaN(var.getResult()));
        var.increment(1d);
        assertEquals(0d, var.getResult(), 0);
    }

    /**
     * Test population version of variance
     */
    @Test
    public void testPopulation() {
        double[] values = { -1.0d, 3.1d, 4.0d, -2.1d, 22d, 11.7d, 3d, 14d };
        SecondMoment m = new SecondMoment();
        m.incrementAll(values);  // side effect is to add values
        Variance v1 = new Variance();
        v1 = v1.withBiasCorrection(false);
        assertEquals(populationVariance(values), v1.evaluate(values), 1E-14);
        v1.incrementAll(values);
        assertEquals(populationVariance(values), v1.getResult(), 1E-14);
        v1 = new Variance(false, m);
        assertEquals(populationVariance(values), v1.getResult(), 1E-14);
        v1 = new Variance(false);
        assertEquals(populationVariance(values), v1.evaluate(values), 1E-14);
        v1.incrementAll(values);
        assertEquals(populationVariance(values), v1.getResult(), 1E-14);
    }

    /**
     * Definitional formula for population variance
     */
    protected double populationVariance(double[] v) {
        double mean = StatUtils.mean(v);
        double sum = 0;
        for (double val : v) {
           sum += (val - mean) * (val - mean);
        }
        return sum / v.length;
    }

    @Test
    public void testWeightedVariance() {
        Variance variance = getUnivariateStatistic();
        assertEquals(expectedWeightedValue(),
                     variance.evaluate(testArray, testWeightsArray, 0, testArray.length),
                     getTolerance());

        // All weights = 1 -> weighted variance = unweighted variance
        assertEquals(expectedValue(),
                     variance.evaluate(testArray, unitWeightsArray, 0, testArray.length),
                     getTolerance());

        // All weights the same -> when weights are normalized to sum to the length of the values array,
        // weighted variance = unweighted value
        assertEquals(expectedValue(),
                     variance.evaluate(testArray,
                                       MathArrays.normalizeArray(identicalWeightsArray, testArray.length),
                                       0, testArray.length),
                     getTolerance());

    }

}
