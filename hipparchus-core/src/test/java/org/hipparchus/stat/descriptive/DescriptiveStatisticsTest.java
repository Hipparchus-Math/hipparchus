/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.stat.descriptive;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Locale;

import org.hipparchus.TestUtils;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.stat.descriptive.moment.GeometricMean;
import org.hipparchus.stat.descriptive.moment.Mean;
import org.hipparchus.stat.descriptive.moment.Variance;
import org.hipparchus.stat.descriptive.rank.Max;
import org.hipparchus.stat.descriptive.rank.Min;
import org.hipparchus.stat.descriptive.rank.Percentile;
import org.hipparchus.stat.descriptive.summary.Sum;
import org.hipparchus.stat.descriptive.summary.SumOfSquares;
import org.hipparchus.util.Precision;
import org.junit.Test;

/**
 * Test cases for the {@link DescriptiveStatistics} class.
 */
public class DescriptiveStatisticsTest {

    protected DescriptiveStatistics createDescriptiveStatistics() {
        return new DescriptiveStatistics();
    }

    protected DescriptiveStatistics.Builder createBuilder() {
        return DescriptiveStatistics.builder();
    }

    @Test
    public void testCustomMeanImpl() {
        DescriptiveStatistics stats = createDescriptiveStatistics();
        stats.addValue(1);
        stats.addValue(3);
        assertEquals(2, stats.getMean(), 1E-10);

        // Now lets try some new math
        DescriptiveStatistics stats2 =
                createBuilder().withMeanImpl(new DeepMean())
                               .withInitialValues(stats.getValues())
                               .build();

        assertEquals(42, stats2.getMean(), 1E-10);
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

        // Try wrapper impl
        DescriptiveStatistics stats1 =
            DescriptiveStatistics.builder()
                                 .withInitialValues(stats.getValues())
                                 .withPercentileImpl(new GoodPercentile())
                                 .build();

        assertEquals(2, stats1.getPercentile(50.0), 1E-10);

        // Try "new math" impl
        DescriptiveStatistics stats2 =
            DescriptiveStatistics.builder()
                                 .withInitialValues(stats.getValues())
                                 .withPercentileImpl(new SubPercentile())
                                 .build();

        assertEquals(10.0, stats2.getPercentile(10.0), 1E-10);
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
        final SummaryStatistics sstats = new SummaryStatistics();
        final double tol = 1E-12;
        for (int i = 0; i < 20; i++) {
            dstats.addValue(i);
            sstats.clear();
            double[] values = dstats.getValues();
            for (int j = 0; j < values.length; j++) {
                sstats.addValue(values[j]);
            }
            TestUtils.assertEquals(dstats.getMean(), sstats.getMean(), tol);
            TestUtils.assertEquals(new Mean().evaluate(values), dstats.getMean(), tol);
            TestUtils.assertEquals(dstats.getMax(), sstats.getMax(), tol);
            TestUtils.assertEquals(new Max().evaluate(values), dstats.getMax(), tol);
            TestUtils.assertEquals(dstats.getGeometricMean(), sstats.getGeometricMean(), tol);
            TestUtils.assertEquals(new GeometricMean().evaluate(values), dstats.getGeometricMean(), tol);
            TestUtils.assertEquals(dstats.getMin(), sstats.getMin(), tol);
            TestUtils.assertEquals(new Min().evaluate(values), dstats.getMin(), tol);
            TestUtils.assertEquals(dstats.getStandardDeviation(), sstats.getStandardDeviation(), tol);
            TestUtils.assertEquals(dstats.getVariance(), sstats.getVariance(), tol);
            TestUtils.assertEquals(new Variance().evaluate(values), dstats.getVariance(), tol);
            TestUtils.assertEquals(dstats.getSum(), sstats.getSum(), tol);
            TestUtils.assertEquals(new Sum().evaluate(values), dstats.getSum(), tol);
            TestUtils.assertEquals(dstats.getSumsq(), sstats.getSumsq(), tol);
            TestUtils.assertEquals(new SumOfSquares().evaluate(values), dstats.getSumsq(), tol);
            TestUtils.assertEquals(dstats.getPopulationVariance(), sstats.getPopulationVariance(), tol);
            TestUtils.assertEquals(new Variance(false).evaluate(values), dstats.getPopulationVariance(), tol);
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

    // Test UnivariateStatistics impls for custom statistics tests.

    /**
     * A new way to compute the mean.
     */
    static class DeepMean implements UnivariateStatistic {

        @Override
        public double evaluate(double[] values, int begin, int length) {
            return 42;
        }

        @Override
        public double evaluate(double[] values) {
            return 42;
        }
        @Override
        public UnivariateStatistic copy() {
            return new DeepMean();
        }
    }

    /**
     * Test percentile implementation - wraps a Percentile.
     */
    static class GoodPercentile implements QuantiledUnivariateStatistic {

        private final Percentile percentile = new Percentile();

        @Override
        public void setQuantile(double quantile) {
            percentile.setQuantile(quantile);
        }

        @Override
        public double getQuantile() {
            return percentile.getQuantile();
        }

        @Override
        public double evaluate(double[] values, int begin, int length, double p) {
            return percentile.evaluate(values, begin, length, p);
        }

        @Override
        public double evaluate(double[] values, int begin, int length) {
            return percentile.evaluate(values, begin, length);
        }

        @Override
        public QuantiledUnivariateStatistic copy() {
            GoodPercentile result = new GoodPercentile();
            result.setQuantile(percentile.getQuantile());
            return result;
        }
    }

    /**
     * Test percentile subclass - another "new math" impl
     * Always returns currently set quantile.
     */
    static class SubPercentile extends Percentile {
        private static final long serialVersionUID = 1L;

        @Override
        public double evaluate(double[] values, int begin, int length, double p) {
            return getQuantile();
        }

        @Override
        public double evaluate(double[] values, double p) {
            return getQuantile();
        }

        @Override
        public Percentile copy() {
            SubPercentile result = new SubPercentile();
            result.setQuantile(this.getQuantile());
            return result;
        }
    }
}
