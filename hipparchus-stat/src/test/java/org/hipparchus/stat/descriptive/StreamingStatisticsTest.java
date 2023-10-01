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
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.distribution.continuous.UniformRealDistribution;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.random.Well19937a;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test cases for the {@link StreamingStatistics} class.
 */
public class StreamingStatisticsTest {

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

    protected StreamingStatistics createStreamingStatistics() {
        return new StreamingStatistics();
    }

    /** test stats */
    @Test
    public void testStats() {
        StreamingStatistics u = createStreamingStatistics();
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
        StreamingStatistics u = createStreamingStatistics();
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
    public void testN0andN1Conditions() {
        StreamingStatistics u = createStreamingStatistics();
        assertTrue("Mean of n = 0 set should be NaN", Double.isNaN( u.getMean() ) );
        assertTrue("Standard Deviation of n = 0 set should be NaN",
                   Double.isNaN( u.getStandardDeviation() ) );
        assertTrue("Variance of n = 0 set should be NaN", Double.isNaN(u.getVariance() ) );

        /* n=1 */
        u.addValue(one);
        assertTrue("mean should be one (n = 1)", u.getMean() == one);
        assertTrue("geometric should be one (n = 1) instead it is " + u.getGeometricMean(),
                   u.getGeometricMean() == one);
        assertTrue("Std should be zero (n = 1)", u.getStandardDeviation() == 0.0);
        assertTrue("variance should be zero (n = 1)", u.getVariance() == 0.0);

        /* n=2 */
        u.addValue(twoF);
        assertTrue("Std should not be zero (n = 2)", u.getStandardDeviation() != 0.0);
        assertTrue("variance should not be zero (n = 2)", u.getVariance() != 0.0);
    }

    @Test
    public void testProductAndGeometricMean() {
        StreamingStatistics u = createStreamingStatistics();
        u.addValue( 1.0 );
        u.addValue( 2.0 );
        u.addValue( 3.0 );
        u.addValue( 4.0 );

        assertEquals( "Geometric mean not expected", 2.213364, u.getGeometricMean(), 0.00001 );
    }

    @Test
    public void testNaNContracts() {
        StreamingStatistics u = createStreamingStatistics();
        assertTrue("mean not NaN",Double.isNaN(u.getMean()));
        assertTrue("min not NaN",Double.isNaN(u.getMin()));
        assertTrue("std dev not NaN",Double.isNaN(u.getStandardDeviation()));
        assertTrue("var not NaN",Double.isNaN(u.getVariance()));
        assertTrue("geom mean not NaN",Double.isNaN(u.getGeometricMean()));

        u.addValue(1.0);

        assertEquals("mean", 1.0, u.getMean(), Double.MIN_VALUE);
        assertEquals("variance", 0.0, u.getVariance(), Double.MIN_VALUE);
        assertEquals("geometric mean", 1.0, u.getGeometricMean(), Double.MIN_VALUE);

        u.addValue(-1.0);

        assertTrue("geom mean not NaN",Double.isNaN(u.getGeometricMean()));

        u.addValue(0.0);

        assertTrue("geom mean not NaN",Double.isNaN(u.getGeometricMean()));

        //FiXME: test all other NaN contract specs
    }

    @Test
    public void testGetSummary() {
        StreamingStatistics u = createStreamingStatistics();
        StatisticalSummary summary = u.getSummary();
        verifySummary(u, summary);
        u.addValue(1d);
        summary = u.getSummary();
        verifySummary(u, summary);
        u.addValue(2d);
        summary = u.getSummary();
        verifySummary(u, summary);
        u.addValue(2d);
        summary = u.getSummary();
        verifySummary(u, summary);
    }

    @Test
    public void testSerialization() {
        StreamingStatistics u = createStreamingStatistics();
        // Empty test
        UnitTestUtils.checkSerializedEquality(u);
        StreamingStatistics s = (StreamingStatistics) UnitTestUtils.serializeAndRecover(u);
        StatisticalSummary summary = s.getSummary();
        verifySummary(u, summary);

        // Add some data
        u.addValue(2d);
        u.addValue(1d);
        u.addValue(3d);
        u.addValue(4d);
        u.addValue(5d);

        // Test again
        UnitTestUtils.checkSerializedEquality(u);
        s = (StreamingStatistics) UnitTestUtils.serializeAndRecover(u);
        summary = s.getSummary();
        verifySummary(u, summary);
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    public void testEqualsAndHashCode() {
        StreamingStatistics u = createStreamingStatistics();
        StreamingStatistics t = null;
        int emptyHash = u.hashCode();
        assertTrue("reflexive", u.equals(u));
        assertFalse("non-null compared to null", u.equals(t));
        assertFalse("wrong type", u.equals(Double.valueOf(0)));
        t = createStreamingStatistics();
        assertTrue("empty instances should be equal", t.equals(u));
        assertTrue("empty instances should be equal", u.equals(t));
        assertEquals("empty hash code", emptyHash, t.hashCode());

        // Add some data to u
        u.addValue(2d);
        u.addValue(1d);
        u.addValue(3d);
        u.addValue(4d);
        assertFalse("different n's should make instances not equal", t.equals(u));
        assertFalse("different n's should make instances not equal", u.equals(t));
        assertTrue("different n's should make hashcodes different", u.hashCode() != t.hashCode());

        //Add data in same order to t
        t.addValue(2d);
        t.addValue(1d);
        t.addValue(3d);
        t.addValue(4d);
        assertTrue("summaries based on same data should be equal", t.equals(u));
        assertTrue("summaries based on same data should be equal", u.equals(t));
        assertEquals("summaries based on same data should have same hashcodes",
                     u.hashCode(), t.hashCode());

        // Clear and make sure summaries are indistinguishable from empty summary
        u.clear();
        t.clear();
        assertTrue("empty instances should be equal", t.equals(u));
        assertTrue("empty instances should be equal", u.equals(t));
        assertEquals("empty hash code", emptyHash, t.hashCode());
        assertEquals("empty hash code", emptyHash, u.hashCode());
    }

    @Test
    public void testCopy() {
        StreamingStatistics u = createStreamingStatistics();
        u.addValue(2d);
        u.addValue(1d);
        u.addValue(3d);
        u.addValue(4d);
        StreamingStatistics v = u.copy();
        assertEquals(u, v);
        assertEquals(v, u);

        // Make sure both behave the same with additional values added
        u.addValue(7d);
        u.addValue(9d);
        u.addValue(11d);
        u.addValue(23d);
        v.addValue(7d);
        v.addValue(9d);
        v.addValue(11d);
        v.addValue(23d);
        assertEquals(u, v);
        assertEquals(v, u);
    }

    private void verifySummary(StreamingStatistics u, StatisticalSummary s) {
        assertEquals("N",s.getN(),u.getN());
        UnitTestUtils.assertEquals("sum",s.getSum(),u.getSum(),tolerance);
        UnitTestUtils.assertEquals("var",s.getVariance(),u.getVariance(),tolerance);
        UnitTestUtils.assertEquals("std",s.getStandardDeviation(),u.getStandardDeviation(),tolerance);
        UnitTestUtils.assertEquals("mean",s.getMean(),u.getMean(),tolerance);
        UnitTestUtils.assertEquals("min",s.getMin(),u.getMin(),tolerance);
        UnitTestUtils.assertEquals("max",s.getMax(),u.getMax(),tolerance);
    }

    @Test
    public void testQuadraticMean() {
        final double[] values = { 1.2, 3.4, 5.6, 7.89 };
        final StreamingStatistics stats = createStreamingStatistics();

        final int len = values.length;
        double expected = 0;
        for (int i = 0; i < len; i++) {
            final double v = values[i];
            expected += v * v / len;

            stats.addValue(v);
        }
        expected = Math.sqrt(expected);

        assertEquals(expected, stats.getQuadraticMean(), Math.ulp(expected));
    }

    @Test
    public void testToString() {
        StreamingStatistics u = createStreamingStatistics();
        for (int i = 0; i < 5; i++) {
            u.addValue(i);
        }
        final String[] labels = {
            "min", "max", "sum", "geometric mean", "variance", "population variance",
            "second moment", "sum of squares", "standard deviation", "sum of logs"
        };
        final double[] values = {
            u.getMin(), u.getMax(), u.getSum(), u.getGeometricMean(), u.getVariance(),
            u.getPopulationVariance(), u.getSecondMoment(), u.getSumOfSquares(),
            u.getStandardDeviation(), u.getSumOfLogs()
        };
        final String toString = u.toString();
        assertTrue(toString.indexOf("n: " + u.getN()) > 0); // getN() returns a long
        for (int i = 0; i < values.length; i++) {
            assertTrue(toString.indexOf(labels[i] + ": " + String.valueOf(values[i])) > 0);
        }
    }

    /**
     * Verify that aggregating over a partition gives the same results
     * as direct computation.
     *
     *  1) Randomly generate a dataset of 10-100 values
     *     from [-100, 100]
     *  2) Divide the dataset it into 2-5 partitions
     *  3) Create an AggregateSummaryStatistic and ContributingStatistics
     *     for each partition
     *  4) Compare results from the AggregateSummaryStatistic with values
     *     returned by a single SummaryStatistics instance that is provided
     *     the full dataset
     */
    @Test
    public void testAggregationConsistency() {

        // Generate a random sample and random partition
        double[] totalSample = generateSample();
        double[][] subSamples = generatePartition(totalSample);
        int nSamples = subSamples.length;

        // Create aggregator and total stats for comparison
        StreamingStatistics aggregate = new StreamingStatistics();
        StreamingStatistics totalStats = new StreamingStatistics();

        // Create array of component stats
        StreamingStatistics componentStats[] = new StreamingStatistics[nSamples];

        for (int i = 0; i < nSamples; i++) {

            // Make componentStats[i] a contributing statistic to aggregate
            componentStats[i] = new StreamingStatistics();

            // Add values from subsample
            for (int j = 0; j < subSamples[i].length; j++) {
                componentStats[i].addValue(subSamples[i][j]);
            }
        }

        aggregate.aggregate(componentStats);

        // Compute totalStats directly
        for (int i = 0; i < totalSample.length; i++) {
            totalStats.addValue(totalSample[i]);
        }

        /*
         * Compare statistics in totalStats with aggregate.
         * Note that guaranteed success of this comparison depends on the
         * fact that <aggregate> gets values in exactly the same order
         * as <totalStats>.
         */
        assertSummaryStatisticsEquals(totalStats, aggregate, 1e-10);

        // Check some percentiles
        final double tol = 1e-13;
        assertEquals(totalStats.getPercentile(10), aggregate.getPercentile(10), tol);
        assertEquals(totalStats.getPercentile(25), aggregate.getPercentile(25), tol);
        assertEquals(totalStats.getPercentile(50), aggregate.getPercentile(50), tol);
        assertEquals(totalStats.getPercentile(75), aggregate.getPercentile(75), tol);
        assertEquals(totalStats.getPercentile(90), aggregate.getPercentile(90), tol);
        assertEquals(totalStats.getPercentile(99), aggregate.getPercentile(99), tol);
    }

    @Test
    public void testAggregateDegenerate() {
        double[] totalSample = {1, 2, 3, 4, 5};
        double[][] subSamples = {{1}, {2}, {3}, {4}, {5}};

        // Compute combined stats directly
        StreamingStatistics totalStats = new StreamingStatistics();
        for (int i = 0; i < totalSample.length; i++) {
            totalStats.addValue(totalSample[i]);
        }

        // Now compute subsample stats individually and aggregate
        StreamingStatistics[] subSampleStats = new StreamingStatistics[5];
        for (int i = 0; i < 5; i++) {
            subSampleStats[i] = new StreamingStatistics();
        }
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < subSamples[i].length; j++) {
                subSampleStats[i].addValue(subSamples[i][j]);
            }
        }

        // Compare values
        StreamingStatistics aggregatedStats = new StreamingStatistics();
        aggregatedStats.aggregate(subSampleStats);

        assertSummaryStatisticsEquals(totalStats, aggregatedStats, 10e-10);
    }

    @Test
    public void testAggregateSpecialValues() {
        double[] totalSample = {Double.POSITIVE_INFINITY, 2, 3, Double.NaN, 5};
        double[][] subSamples = {{Double.POSITIVE_INFINITY, 2}, {3}, {Double.NaN}, {5}};

        // Compute combined stats directly
        StreamingStatistics totalStats = new StreamingStatistics();
        for (int i = 0; i < totalSample.length; i++) {
            totalStats.addValue(totalSample[i]);
        }

        // Now compute subsample stats individually and aggregate
        StreamingStatistics[] subSampleStats = new StreamingStatistics[4];
        for (int i = 0; i < 4; i++) {
            subSampleStats[i] = new StreamingStatistics();
        }
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < subSamples[i].length; j++) {
                subSampleStats[i].addValue(subSamples[i][j]);
            }
        }

        // Compare values
        StreamingStatistics aggregatedStats = new StreamingStatistics();
        aggregatedStats.aggregate(subSampleStats);

        assertSummaryStatisticsEquals(totalStats, aggregatedStats, 10e-10);
    }

    @Test
    public void testBuilderDefault() {
       StreamingStatistics stats = StreamingStatistics.builder().build();
       stats.addValue(10);
       stats.addValue(20);
       stats.addValue(30);
       // Percentiles should be NaN, all others should have values
       Assert.assertFalse(Double.isNaN(stats.getMax()));
       Assert.assertFalse(Double.isNaN(stats.getMin()));
       Assert.assertFalse(Double.isNaN(stats.getMean()));
       Assert.assertFalse(Double.isNaN(stats.getSum()));
       Assert.assertFalse(Double.isNaN(stats.getVariance()));
       Assert.assertFalse(Double.isNaN(stats.getPopulationVariance()));
       Assert.assertFalse(Double.isNaN(stats.getStandardDeviation()));
       Assert.assertFalse(Double.isNaN(stats.getGeometricMean()));
       Assert.assertFalse(Double.isNaN(stats.getQuadraticMean()));
       Assert.assertFalse(Double.isNaN(stats.getSumOfSquares()));
       Assert.assertFalse(Double.isNaN(stats.getSumOfLogs()));
       Assert.assertTrue(Double.isNaN(stats.getMedian()));
       Assert.assertTrue(Double.isNaN(stats.getPercentile(10)));
    }

    @Test
    public void testBuilderPercentilesOn() {
        StreamingStatistics stats = StreamingStatistics.
                builder().
                percentiles(1.0e-6, new Well19937a(0x9b1fadc49a76102al)).
                build();
        stats.addValue(10);
        stats.addValue(20);
        stats.addValue(30);
        Assert.assertFalse(Double.isNaN(stats.getMax()));
        Assert.assertFalse(Double.isNaN(stats.getMin()));
        Assert.assertFalse(Double.isNaN(stats.getMean()));
        Assert.assertFalse(Double.isNaN(stats.getSum()));
        Assert.assertFalse(Double.isNaN(stats.getVariance()));
        Assert.assertFalse(Double.isNaN(stats.getPopulationVariance()));
        Assert.assertFalse(Double.isNaN(stats.getStandardDeviation()));
        Assert.assertFalse(Double.isNaN(stats.getGeometricMean()));
        Assert.assertFalse(Double.isNaN(stats.getQuadraticMean()));
        Assert.assertFalse(Double.isNaN(stats.getSumOfSquares()));
        Assert.assertFalse(Double.isNaN(stats.getSumOfLogs()));
        Assert.assertFalse(Double.isNaN(stats.getMedian()));
        Assert.assertFalse(Double.isNaN(stats.getPercentile(10)));
        stats.clear();
        Assert.assertTrue(Double.isNaN(stats.getMax()));
        Assert.assertTrue(Double.isNaN(stats.getMin()));
        Assert.assertTrue(Double.isNaN(stats.getMean()));
        Assert.assertEquals(0.0, stats.getSum(), 1.0e-15);
        Assert.assertTrue(Double.isNaN(stats.getVariance()));
        Assert.assertTrue(Double.isNaN(stats.getPopulationVariance()));
        Assert.assertTrue(Double.isNaN(stats.getStandardDeviation()));
        Assert.assertTrue(Double.isNaN(stats.getGeometricMean()));
        Assert.assertTrue(Double.isNaN(stats.getQuadraticMean()));
        Assert.assertEquals(0.0, stats.getSumOfSquares(), 1.0e-15);
        Assert.assertEquals(0.0, stats.getSumOfLogs(), 1.0e-15);
        Assert.assertTrue(Double.isNaN(stats.getMedian()));
        Assert.assertTrue(Double.isNaN(stats.getPercentile(10)));
        Assert.assertTrue(Double.isNaN(stats.getPercentile(10)));
    }

    @Test
    public void testBuilderMomentsOff() {
        StreamingStatistics stats = StreamingStatistics.
                builder().
                percentiles(1.0e-6, new Well19937a(0x9b1fadc49a76102al)).
                moments(false).
                build();
        stats.addValue(10);
        stats.addValue(20);
        stats.addValue(30);
        Assert.assertFalse(Double.isNaN(stats.getMax()));
        Assert.assertFalse(Double.isNaN(stats.getMin()));
        Assert.assertTrue(Double.isNaN(stats.getMean()));
        Assert.assertTrue(Double.isNaN(stats.getSum()));
        Assert.assertTrue(Double.isNaN(stats.getVariance()));
        Assert.assertTrue(Double.isNaN(stats.getPopulationVariance()));
        Assert.assertTrue(Double.isNaN(stats.getStandardDeviation()));
        Assert.assertFalse(Double.isNaN(stats.getGeometricMean()));
        Assert.assertFalse(Double.isNaN(stats.getQuadraticMean()));
        Assert.assertFalse(Double.isNaN(stats.getSumOfSquares()));
        Assert.assertFalse(Double.isNaN(stats.getSumOfLogs()));
        Assert.assertFalse(Double.isNaN(stats.getMedian()));
        Assert.assertFalse(Double.isNaN(stats.getPercentile(10)));
        stats.clear();
        Assert.assertTrue(Double.isNaN(stats.getMax()));
        Assert.assertTrue(Double.isNaN(stats.getMin()));
        Assert.assertTrue(Double.isNaN(stats.getMean()));
        Assert.assertTrue(Double.isNaN(stats.getSum()));
        Assert.assertTrue(Double.isNaN(stats.getVariance()));
        Assert.assertTrue(Double.isNaN(stats.getPopulationVariance()));
        Assert.assertTrue(Double.isNaN(stats.getStandardDeviation()));
        Assert.assertTrue(Double.isNaN(stats.getGeometricMean()));
        Assert.assertTrue(Double.isNaN(stats.getQuadraticMean()));
        Assert.assertEquals(0.0, stats.getSumOfSquares(), 1.0e-15);
        Assert.assertEquals(0.0, stats.getSumOfLogs(), 1.0e-15);
        Assert.assertTrue(Double.isNaN(stats.getMedian()));
        Assert.assertTrue(Double.isNaN(stats.getPercentile(10)));
        Assert.assertTrue(Double.isNaN(stats.getPercentile(10)));
    }

    @Test
    public void testBuilderSumOfLogsOff() {
        StreamingStatistics stats = StreamingStatistics.
                builder().
                percentiles(1.0e-6, new Well19937a(0x9b1fadc49a76102al)).
                sumOfLogs(false).
                build();
        stats.addValue(10);
        stats.addValue(20);
        stats.addValue(30);
        Assert.assertFalse(Double.isNaN(stats.getMax()));
        Assert.assertFalse(Double.isNaN(stats.getMin()));
        Assert.assertFalse(Double.isNaN(stats.getMean()));
        Assert.assertFalse(Double.isNaN(stats.getSum()));
        Assert.assertFalse(Double.isNaN(stats.getVariance()));
        Assert.assertFalse(Double.isNaN(stats.getPopulationVariance()));
        Assert.assertFalse(Double.isNaN(stats.getStandardDeviation()));
        Assert.assertTrue(Double.isNaN(stats.getGeometricMean()));
        Assert.assertFalse(Double.isNaN(stats.getQuadraticMean()));
        Assert.assertFalse(Double.isNaN(stats.getSumOfSquares()));
        Assert.assertTrue(Double.isNaN(stats.getSumOfLogs()));
        Assert.assertFalse(Double.isNaN(stats.getMedian()));
        Assert.assertFalse(Double.isNaN(stats.getPercentile(10)));
        stats.clear();
        Assert.assertTrue(Double.isNaN(stats.getMax()));
        Assert.assertTrue(Double.isNaN(stats.getMin()));
        Assert.assertTrue(Double.isNaN(stats.getMean()));
        Assert.assertEquals(0.0, stats.getSum(), 1.0e-15);
        Assert.assertTrue(Double.isNaN(stats.getVariance()));
        Assert.assertTrue(Double.isNaN(stats.getPopulationVariance()));
        Assert.assertTrue(Double.isNaN(stats.getStandardDeviation()));
        Assert.assertTrue(Double.isNaN(stats.getGeometricMean()));
        Assert.assertTrue(Double.isNaN(stats.getQuadraticMean()));
        Assert.assertEquals(0.0, stats.getSumOfSquares(), 1.0e-15);
        Assert.assertTrue(Double.isNaN(stats.getSumOfLogs()));
        Assert.assertTrue(Double.isNaN(stats.getMedian()));
        Assert.assertTrue(Double.isNaN(stats.getPercentile(10)));
        Assert.assertTrue(Double.isNaN(stats.getPercentile(10)));
    }

    @Test
    public void testBuilderExtremaOff() {
        StreamingStatistics stats = StreamingStatistics.
                builder().
                percentiles(1.0e-6, new Well19937a(0x9b1fadc49a76102al)).
                extrema(false).
                build();
        stats.addValue(10);
        stats.addValue(20);
        stats.addValue(30);
        Assert.assertTrue(Double.isNaN(stats.getMax()));
        Assert.assertTrue(Double.isNaN(stats.getMin()));
        Assert.assertFalse(Double.isNaN(stats.getMean()));
        Assert.assertFalse(Double.isNaN(stats.getSum()));
        Assert.assertFalse(Double.isNaN(stats.getVariance()));
        Assert.assertFalse(Double.isNaN(stats.getPopulationVariance()));
        Assert.assertFalse(Double.isNaN(stats.getStandardDeviation()));
        Assert.assertFalse(Double.isNaN(stats.getGeometricMean()));
        Assert.assertFalse(Double.isNaN(stats.getQuadraticMean()));
        Assert.assertFalse(Double.isNaN(stats.getSumOfSquares()));
        Assert.assertFalse(Double.isNaN(stats.getSumOfLogs()));
        Assert.assertFalse(Double.isNaN(stats.getMedian()));
        Assert.assertFalse(Double.isNaN(stats.getPercentile(10)));
        stats.clear();
        Assert.assertTrue(Double.isNaN(stats.getMax()));
        Assert.assertTrue(Double.isNaN(stats.getMin()));
        Assert.assertTrue(Double.isNaN(stats.getMean()));
        Assert.assertEquals(0.0, stats.getSum(), 1.0e-15);
        Assert.assertTrue(Double.isNaN(stats.getVariance()));
        Assert.assertTrue(Double.isNaN(stats.getPopulationVariance()));
        Assert.assertTrue(Double.isNaN(stats.getStandardDeviation()));
        Assert.assertTrue(Double.isNaN(stats.getGeometricMean()));
        Assert.assertTrue(Double.isNaN(stats.getQuadraticMean()));
        Assert.assertEquals(0.0, stats.getSumOfSquares(), 1.0e-15);
        Assert.assertEquals(0.0, stats.getSumOfLogs(), 1.0e-15);
        Assert.assertTrue(Double.isNaN(stats.getMedian()));
        Assert.assertTrue(Double.isNaN(stats.getPercentile(10)));
        Assert.assertTrue(Double.isNaN(stats.getPercentile(10)));
    }

    @Test
    public void testBuilderSumOfSquares() {
        StreamingStatistics stats = StreamingStatistics.
                builder().
                percentiles(1.0e-6, new Well19937a(0x9b1fadc49a76102al)).
                sumOfSquares(false).
                build();
        stats.addValue(10);
        stats.addValue(20);
        stats.addValue(30);
        Assert.assertFalse(Double.isNaN(stats.getMax()));
        Assert.assertFalse(Double.isNaN(stats.getMin()));
        Assert.assertFalse(Double.isNaN(stats.getMean()));
        Assert.assertFalse(Double.isNaN(stats.getSum()));
        Assert.assertFalse(Double.isNaN(stats.getVariance()));
        Assert.assertFalse(Double.isNaN(stats.getPopulationVariance()));
        Assert.assertFalse(Double.isNaN(stats.getStandardDeviation()));
        Assert.assertFalse(Double.isNaN(stats.getGeometricMean()));
        Assert.assertTrue(Double.isNaN(stats.getQuadraticMean()));
        Assert.assertTrue(Double.isNaN(stats.getSumOfSquares()));
        Assert.assertFalse(Double.isNaN(stats.getSumOfLogs()));
        Assert.assertFalse(Double.isNaN(stats.getMedian()));
        Assert.assertFalse(Double.isNaN(stats.getPercentile(10)));
        stats.clear();
        Assert.assertTrue(Double.isNaN(stats.getMax()));
        Assert.assertTrue(Double.isNaN(stats.getMin()));
        Assert.assertTrue(Double.isNaN(stats.getMean()));
        Assert.assertEquals(0.0, stats.getSum(), 1.0e-15);
        Assert.assertTrue(Double.isNaN(stats.getVariance()));
        Assert.assertTrue(Double.isNaN(stats.getPopulationVariance()));
        Assert.assertTrue(Double.isNaN(stats.getStandardDeviation()));
        Assert.assertTrue(Double.isNaN(stats.getGeometricMean()));
        Assert.assertTrue(Double.isNaN(stats.getQuadraticMean()));
        Assert.assertTrue(Double.isNaN(stats.getSumOfSquares()));
        Assert.assertEquals(0.0, stats.getSumOfLogs(), 1.0e-15);
        Assert.assertTrue(Double.isNaN(stats.getMedian()));
        Assert.assertTrue(Double.isNaN(stats.getPercentile(10)));
        Assert.assertTrue(Double.isNaN(stats.getPercentile(10)));
    }

    /**
     * Verifies that a StatisticalSummary and a StatisticalSummaryValues are equal up
     * to delta, with NaNs, infinities returned in the same spots. For max, min, n, values
     * have to agree exactly, delta is used only for sum, mean, variance, std dev.
     */
    protected static void assertSummaryStatisticsEquals(StreamingStatistics expected,
                                                        StreamingStatistics observed,
                                                        double delta) {
        UnitTestUtils.assertEquals(expected.getMax(), observed.getMax(), 0);
        UnitTestUtils.assertEquals(expected.getMin(), observed.getMin(), 0);
        Assert.assertEquals(expected.getN(), observed.getN());
        UnitTestUtils.assertEquals(expected.getSum(), observed.getSum(), delta);
        UnitTestUtils.assertEquals(expected.getMean(), observed.getMean(), delta);
        UnitTestUtils.assertEquals(expected.getStandardDeviation(), observed.getStandardDeviation(), delta);
        UnitTestUtils.assertEquals(expected.getVariance(), observed.getVariance(), delta);
    }


    /**
     * Generates a random sample of double values.
     * Sample size is random, between 10 and 100 and values are
     * uniformly distributed over [-100, 100].
     *
     * @return array of random double values
     */
    private double[] generateSample() {
        final RealDistribution uniformDist = new UniformRealDistribution(-100, 100);
        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(100);
        final int sampleSize = randomDataGenerator.nextInt(10, 100);
        final double[] out = randomDataGenerator.nextDeviates(uniformDist, sampleSize);
        return out;
    }

    /**
     * Generates a partition of <sample> into up to 5 sequentially selected
     * subsamples with randomly selected partition points.
     *
     * @param sample array to partition
     * @return rectangular array with rows = subsamples
     */
    private double[][] generatePartition(double[] sample) {
        final int length = sample.length;
        final double[][] out = new double[5][];
        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(100);
        int cur = 0;          // beginning of current partition segment
        int offset = 0;       // end of current partition segment
        int sampleCount = 0;  // number of segments defined
        for (int i = 0; i < 5; i++) {
            if (cur == length || offset == length) {
                break;
            }
            final int next;
            if (i == 4 || cur == length - 1) {
                next = length - 1;
            } else {
                next = randomDataGenerator.nextInt(cur, length - 1);
            }
            final int subLength = next - cur + 1;
            out[i] = new double[subLength];
            System.arraycopy(sample, offset, out[i], 0, subLength);
            cur = next + 1;
            sampleCount++;
            offset += subLength;
        }
        if (sampleCount < 5) {
            double[][] out2 = new double[sampleCount][];
            for (int j = 0; j < sampleCount; j++) {
                final int curSize = out[j].length;
                out2[j] = new double[curSize];
                System.arraycopy(out[j], 0, out2[j], 0, curSize);
            }
            return out2;
        } else {
            return out;
        }
    }

}
