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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;
import java.util.Locale;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.stat.descriptive.moment.GeometricMean;
import org.hipparchus.stat.descriptive.moment.Mean;
import org.hipparchus.stat.descriptive.moment.Variance;
import org.hipparchus.stat.descriptive.rank.Max;
import org.hipparchus.stat.descriptive.rank.Min;
import org.hipparchus.stat.descriptive.summary.Sum;
import org.hipparchus.stat.descriptive.summary.SumOfSquares;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.Test;

/**
 * Test cases for the {@link DescriptiveStatistics} class.
 */
public class DescriptiveStatisticsTest {

    private final double[] testArray = new double[] { 1, 2, 2, 3 };

    private final double one = 1;
    private final float twoF = 2;
    private final long twoL = 2;
    private final int three = 3;
    private final double mean = 2;
    private final double sumSq = 18;
    private final double sum = 8;
    private final double var = 0.666666666666666666667;
    private final double popVar = 0.5;
    private final double std = FastMath.sqrt(var);
    private final double n = 4;
    private final double min = 1;
    private final double max = 3;
    private final double tolerance = 10E-15;

    protected DescriptiveStatistics createDescriptiveStatistics() {
        return new DescriptiveStatistics();
    }

    /** test stats */
    @Test
    public void testStats() {
        DescriptiveStatistics u = createDescriptiveStatistics();
        assertEquals("total count", 0, u.getN(), tolerance);
        u.addValue(one);
        u.addValue(twoF);
        u.addValue(twoL);
        u.addValue(three);
        assertEquals("N", n, u.getN(), tolerance);
        assertEquals("sum", sum, u.getSum(), tolerance);
        assertEquals("sumsq", sumSq, u.getSumOfSquares(), tolerance);
        assertEquals("var", var, u.getVariance(), tolerance);
        assertEquals("population var", popVar, u.getPopulationVariance(), tolerance);
        assertEquals("std", std, u.getStandardDeviation(), tolerance);
        assertEquals("mean", mean, u.getMean(), tolerance);
        assertEquals("min", min, u.getMin(), tolerance);
        assertEquals("max", max, u.getMax(), tolerance);
        u.clear();
        assertEquals("total count", 0, u.getN(), tolerance);
    }

    @Test
    public void testConsume() {
        DescriptiveStatistics u = createDescriptiveStatistics();
        assertEquals("total count", 0, u.getN(), tolerance);

        Arrays.stream(testArray)
              .forEach(u);

        assertEquals("N", n, u.getN(), tolerance);
        assertEquals("sum", sum, u.getSum(), tolerance);
        assertEquals("sumsq", sumSq, u.getSumOfSquares(), tolerance);
        assertEquals("var", var, u.getVariance(), tolerance);
        assertEquals("population var", popVar, u.getPopulationVariance(), tolerance);
        assertEquals("std", std, u.getStandardDeviation(), tolerance);
        assertEquals("mean", mean, u.getMean(), tolerance);
        assertEquals("min", min, u.getMin(), tolerance);
        assertEquals("max", max, u.getMax(), tolerance);
        u.clear();
        assertEquals("total count", 0, u.getN(), tolerance);
    }

    @Test
    public void testCopy() {
        DescriptiveStatistics stats = createDescriptiveStatistics();
        stats.addValue(1);
        stats.addValue(3);
        assertEquals(2, stats.getMean(), 1E-10);
        DescriptiveStatistics copy = stats.copy();
        assertEquals(2, copy.getMean(), 1E-10);
    }

    @Test
    public void testWindowSize() {
        DescriptiveStatistics stats = createDescriptiveStatistics();
        stats.setWindowSize(300);
        for (int i = 0; i < 100; ++i) {
            stats.addValue(i + 1);
        }
        int refSum = (100 * 101) / 2;
        assertEquals(refSum / 100.0, stats.getMean(), 1E-10);
        assertEquals(300, stats.getWindowSize());
        try {
            stats.setWindowSize(-3);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException iae) {
            // expected
        }
        assertEquals(300, stats.getWindowSize());
        stats.setWindowSize(50);
        assertEquals(50, stats.getWindowSize());
        int refSum2 = refSum - (50 * 51) / 2;
        assertEquals(refSum2 / 50.0, stats.getMean(), 1E-10);
    }

    @Test
    public void testGetValues() {
        DescriptiveStatistics stats = createDescriptiveStatistics();
        for (int i = 100; i > 0; --i) {
            stats.addValue(i);
        }
        int refSum = (100 * 101) / 2;
        assertEquals(refSum / 100.0, stats.getMean(), 1E-10);
        double[] v = stats.getValues();
        for (int i = 0; i < v.length; ++i) {
            assertEquals(100.0 - i, v[i], 1.0e-10);
        }
        double[] s = stats.getSortedValues();
        for (int i = 0; i < s.length; ++i) {
            assertEquals(i + 1.0, s[i], 1.0e-10);
        }
        assertEquals(12.0, stats.getElement(88), 1.0e-10);
    }

    @Test
    public void testQuadraticMean() {
        final double[] values = { 1.2, 3.4, 5.6, 7.89 };
        final DescriptiveStatistics stats = new DescriptiveStatistics(values);

        final int len = values.length;
        double expected = 0;
        for (int i = 0; i < len; i++) {
            final double v = values[i];
            expected += v * v / len;
        }
        expected = Math.sqrt(expected);

        assertEquals(expected, stats.getQuadraticMean(), Math.ulp(expected));
    }

    @Test
    public void testToString() {
        DescriptiveStatistics stats = createDescriptiveStatistics();
        stats.addValue(1);
        stats.addValue(2);
        stats.addValue(3);
        Locale d = Locale.getDefault();
        Locale.setDefault(Locale.US);
        assertEquals("DescriptiveStatistics:\n" +
                     "n: 3\n" +
                     "min: 1.0\n" +
                     "max: 3.0\n" +
                     "mean: 2.0\n" +
                     "std dev: 1.0\n" +
                     "median: 2.0\n" +
                     "skewness: 0.0\n" +
                     "kurtosis: NaN\n",  stats.toString());
        Locale.setDefault(d);
    }

    @Test
    public void testPercentile() {
        DescriptiveStatistics stats = createDescriptiveStatistics();

        stats.addValue(1);
        stats.addValue(2);
        stats.addValue(3);
        assertEquals(2, stats.getPercentile(50.0), 1E-10);
    }

    @Test
    public void test20090720() {
        DescriptiveStatistics descriptiveStatistics = new DescriptiveStatistics(100);
        for (int i = 0; i < 161; i++) {
            descriptiveStatistics.addValue(1.2);
        }
        descriptiveStatistics.clear();
        descriptiveStatistics.addValue(1.2);
        assertEquals(1, descriptiveStatistics.getN());
    }

    @Test
    public void testRemoval() {
        final DescriptiveStatistics dstat = createDescriptiveStatistics();

        checkRemoval(dstat, 1, 6.0, 0.0, Double.NaN);
        checkRemoval(dstat, 3, 5.0, 3.0, 4.5);
        checkRemoval(dstat, 6, 3.5, 2.5, 3.0);
        checkRemoval(dstat, 9, 3.5, 2.5, 3.0);
        checkRemoval(dstat, DescriptiveStatistics.INFINITE_WINDOW, 3.5, 2.5, 3.0);
    }

    @Test
    public void testSummaryConsistency() {
        final int windowSize = 5;
        final DescriptiveStatistics dstats = new DescriptiveStatistics(windowSize);
        final StreamingStatistics sstats = new StreamingStatistics();
        final double tol = 1E-12;
        for (int i = 0; i < 20; i++) {
            dstats.addValue(i);
            sstats.clear();
            double[] values = dstats.getValues();
            for (int j = 0; j < values.length; j++) {
                sstats.addValue(values[j]);
            }
            UnitTestUtils.assertEquals(dstats.getMean(), sstats.getMean(), tol);
            UnitTestUtils.assertEquals(new Mean().evaluate(values), dstats.getMean(), tol);
            UnitTestUtils.assertEquals(dstats.getMax(), sstats.getMax(), tol);
            UnitTestUtils.assertEquals(new Max().evaluate(values), dstats.getMax(), tol);
            UnitTestUtils.assertEquals(dstats.getGeometricMean(), sstats.getGeometricMean(), tol);
            UnitTestUtils.assertEquals(new GeometricMean().evaluate(values), dstats.getGeometricMean(), tol);
            UnitTestUtils.assertEquals(dstats.getMin(), sstats.getMin(), tol);
            UnitTestUtils.assertEquals(new Min().evaluate(values), dstats.getMin(), tol);
            UnitTestUtils.assertEquals(dstats.getStandardDeviation(), sstats.getStandardDeviation(), tol);
            UnitTestUtils.assertEquals(dstats.getVariance(), sstats.getVariance(), tol);
            UnitTestUtils.assertEquals(new Variance().evaluate(values), dstats.getVariance(), tol);
            UnitTestUtils.assertEquals(dstats.getSum(), sstats.getSum(), tol);
            UnitTestUtils.assertEquals(new Sum().evaluate(values), dstats.getSum(), tol);
            UnitTestUtils.assertEquals(dstats.getSumOfSquares(), sstats.getSumOfSquares(), tol);
            UnitTestUtils.assertEquals(new SumOfSquares().evaluate(values), dstats.getSumOfSquares(), tol);
            UnitTestUtils.assertEquals(dstats.getPopulationVariance(), sstats.getPopulationVariance(), tol);
            UnitTestUtils.assertEquals(new Variance(false).evaluate(values), dstats.getPopulationVariance(), tol);
        }
    }

    @Test
    public void testMath1129(){
        final double[] data = new double[] {
            -0.012086732064244697,
            -0.24975668704012527,
            0.5706168483164684,
            -0.322111769955327,
            0.24166759508327315,
            Double.NaN,
            0.16698443218942854,
            -0.10427763937565114,
            -0.15595963093172435,
            -0.028075857595882995,
            -0.24137994506058857,
            0.47543170476574426,
            -0.07495595384947631,
            0.37445697625436497,
            -0.09944199541668033
        };

        final DescriptiveStatistics ds = new DescriptiveStatistics(data);

        final double t = ds.getPercentile(75);
        final double o = ds.getPercentile(25);

        final double iqr = t - o;
        // System.out.println(String.format("25th percentile %s 75th percentile %s", o, t));
        assertTrue(iqr >= 0);
    }

    public void checkRemoval(DescriptiveStatistics dstat, int wsize,
                             double mean1, double mean2, double mean3) {

        dstat.setWindowSize(wsize);
        dstat.clear();

        for (int i = 1 ; i <= 6 ; ++i) {
            dstat.addValue(i);
        }

        assertTrue(Precision.equalsIncludingNaN(mean1, dstat.getMean()));
        dstat.replaceMostRecentValue(0);
        assertTrue(Precision.equalsIncludingNaN(mean2, dstat.getMean()));
        dstat.removeMostRecentValue();
        assertTrue(Precision.equalsIncludingNaN(mean3, dstat.getMean()));
    }

}
