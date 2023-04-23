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
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.stat.StatUtils;
import org.hipparchus.stat.descriptive.UnivariateStatisticAbstractTest;
import org.junit.Before;
import org.junit.Test;

/**
 * Test cases for the {@link SemiVariance} class.
 */
public class SemiVarianceTest extends UnivariateStatisticAbstractTest {

    private double semiVariance;

    @Before
    public void setup() {
        // calculate the semivariance the same way as defined by
        // the SemiVariance class. This is not the same as calculating
        // the variance of the values that are below / above the cutoff.
        double val =
            Arrays.stream(testArray)
                  .filter(x -> x < this.mean)
                  .reduce(0, (x, y) -> x + (y - this.mean) * (y - this.mean));

        this.semiVariance = val / (testArray.length - 1);
    }

    @Override
    public SemiVariance getUnivariateStatistic() {
        return new SemiVariance();
    }

    @Override
    public double expectedValue() {
        return semiVariance;
    }

    @Test
    public void testInsufficientData() {
        SemiVariance sv = getUnivariateStatistic();
        try {
            sv.evaluate(null);
            fail("null is not a valid data array.");
        } catch (NullArgumentException nae) {
        }

        try {
            sv = sv.withVarianceDirection(SemiVariance.UPSIDE_VARIANCE);
            sv.evaluate(null);
            fail("null is not a valid data array.");
        } catch (NullArgumentException nae) {
        }
        assertTrue(Double.isNaN(sv.evaluate(new double[0])));
    }

    @Test
    public void testSingleDown() {
        SemiVariance sv = new SemiVariance();
        double[] values = { 50.0d };
        double singletest = sv.evaluate(values);
        assertEquals(0.0d, singletest, 0);
    }

    @Test
    public void testSingleUp() {
        SemiVariance sv = new SemiVariance(SemiVariance.UPSIDE_VARIANCE);
        double[] values = { 50.0d };
        double singletest = sv.evaluate(values);
        assertEquals(0.0d, singletest, 0);
    }

    @Test
    public void testSample() {
        final double[] values = { -2.0d, 2.0d, 4.0d, -2.0d, 22.0d, 11.0d, 3.0d, 14.0d, 5.0d };
        final int length = values.length;
        final double mean = StatUtils.mean(values); // 6.333...
        SemiVariance sv = getUnivariateStatistic(); // Default bias correction is true
        final double downsideSemiVariance = sv.evaluate(values); // Downside is the default
        assertEquals(UnitTestUtils.sumSquareDev(new double[] {-2d, 2d, 4d, -2d, 3d, 5d}, mean) / (length - 1),
                     downsideSemiVariance, 1E-14);

        sv = sv.withVarianceDirection(SemiVariance.UPSIDE_VARIANCE);
        final double upsideSemiVariance = sv.evaluate(values);
        assertEquals(UnitTestUtils.sumSquareDev(new double[] {22d, 11d, 14d}, mean) / (length - 1),
                     upsideSemiVariance, 1E-14);

        // Verify that upper + lower semivariance against the mean sum to variance
        assertEquals(StatUtils.variance(values), downsideSemiVariance + upsideSemiVariance, 10e-12);
    }

    @Test
    public void testPopulation() {
        double[] values = { -2.0d, 2.0d, 4.0d, -2.0d, 22.0d, 11.0d, 3.0d, 14.0d, 5.0d };
        SemiVariance sv = new SemiVariance(false);

        double singletest = sv.evaluate(values);
        assertEquals(19.556d, singletest, 0.01d);

        sv = sv.withVarianceDirection(SemiVariance.UPSIDE_VARIANCE);
        singletest = sv.evaluate(values);
        assertEquals(36.222d, singletest, 0.01d);
    }

    @Test
    public void testNonMeanCutoffs() {
        double[] values = { -2.0d, 2.0d, 4.0d, -2.0d, 22.0d, 11.0d, 3.0d, 14.0d, 5.0d };
        SemiVariance sv = new SemiVariance(false); // Turn off bias correction - use df = length

        double singletest = sv.evaluate(values, 1.0d, SemiVariance.DOWNSIDE_VARIANCE, false, 0, values.length);
        assertEquals(UnitTestUtils.sumSquareDev(new double[] { -2d, -2d }, 1.0d) / values.length,
                     singletest, 0.01d);

        singletest = sv.evaluate(values, 3.0d, SemiVariance.UPSIDE_VARIANCE, false, 0, values.length);
        assertEquals(UnitTestUtils.sumSquareDev(new double[] { 4d, 22d, 11d, 14d, 5d }, 3.0d) / values.length,
                     singletest,
                     0.01d);
    }

    /**
     * Check that the lower + upper semivariance against the mean sum to the variance.
     */
    @Test
    public void testVarianceDecompMeanCutoff() {
        double[] values = { -2.0d, 2.0d, 4.0d, -2.0d, 22.0d, 11.0d, 3.0d, 14.0d, 5.0d };
        double variance = StatUtils.variance(values);
        SemiVariance sv = new SemiVariance(true); // Bias corrected
        sv = sv.withVarianceDirection(SemiVariance.DOWNSIDE_VARIANCE);
        final double lower = sv.evaluate(values);
        sv = sv.withVarianceDirection(SemiVariance.UPSIDE_VARIANCE);
        final double upper = sv.evaluate(values);
        assertEquals(variance, lower + upper, 10e-12);
    }

    /**
     * Check that upper and lower semivariances against a cutoff sum to the sum
     * of squared deviations of the full set of values against the cutoff
     * divided by df = length - 1 (assuming bias-corrected).
     */
    @Test
    public void testVarianceDecompNonMeanCutoff() {
        double[] values = { -2.0d, 2.0d, 4.0d, -2.0d, 22.0d, 11.0d, 3.0d, 14.0d, 5.0d };
        double target = 0;
        double totalSumOfSquares = UnitTestUtils.sumSquareDev(values, target);
        SemiVariance sv = new SemiVariance(true); // Bias corrected
        sv = sv.withVarianceDirection(SemiVariance.DOWNSIDE_VARIANCE);
        double lower = sv.evaluate(values, target);
        sv = sv.withVarianceDirection(SemiVariance.UPSIDE_VARIANCE);
        double upper = sv.evaluate(values, target);
        assertEquals(totalSumOfSquares / (values.length - 1), lower + upper, 10e-12);
    }

    @Test
    public void testNoVariance() {
        final double[] values = {100d, 100d, 100d, 100d};
        SemiVariance sv = getUnivariateStatistic();
        assertEquals(0, sv.evaluate(values), 10E-12);
        assertEquals(0, sv.evaluate(values, 100d), 10E-12);
        assertEquals(0, sv.evaluate(values, 100d, SemiVariance.UPSIDE_VARIANCE, false, 0, values.length), 10E-12);
    }

}
