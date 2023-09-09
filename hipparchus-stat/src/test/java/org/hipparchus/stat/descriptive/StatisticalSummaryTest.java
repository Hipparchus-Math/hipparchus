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

import java.util.ArrayList;
import java.util.Collection;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.distribution.continuous.UniformRealDistribution;
import org.hipparchus.random.RandomDataGenerator;
import org.junit.Assert;
import org.junit.Test;


/**
 * Test cases for {@link StatisticalSummary}.
 */
public class StatisticalSummaryTest {

    /**
     * Test aggregate function by randomly generating a dataset of 10-100 values
     * from [-100, 100], dividing it into 2-5 partitions, computing stats for each
     * partition and comparing the result of aggregate(...) applied to the collection
     * of per-partition SummaryStatistics with a single SummaryStatistics computed
     * over the full sample.
     */
    @Test
    public void testAggregate() {

        // Generate a random sample and random partition
        double[] totalSample = generateSample();
        double[][] subSamples = generatePartition(totalSample);
        int nSamples = subSamples.length;

        // Compute combined stats directly
        StreamingStatistics totalStats = new StreamingStatistics();
        for (int i = 0; i < totalSample.length; i++) {
            totalStats.addValue(totalSample[i]);
        }

        // Now compute subsample stats individually and aggregate
        StreamingStatistics[] subSampleStats = new StreamingStatistics[nSamples];
        for (int i = 0; i < nSamples; i++) {
            subSampleStats[i] = new StreamingStatistics();
        }
        Collection<StreamingStatistics> aggregate = new ArrayList<StreamingStatistics>();
        for (int i = 0; i < nSamples; i++) {
            for (int j = 0; j < subSamples[i].length; j++) {
                subSampleStats[i].addValue(subSamples[i][j]);
            }
            aggregate.add(subSampleStats[i]);
        }

        // Compare values
        StatisticalSummary aggregatedStats = StatisticalSummary.aggregate(aggregate);
        assertStatisticalSummaryEquals(totalStats.getSummary(), aggregatedStats, 10E-12);
    }

    /**
     * Similar to {@link #testAggregate()} but operating on
     * {@link StatisticalSummary} instead.
     */
    @Test
    public void testAggregateStatisticalSummary() {

        // Generate a random sample and random partition
        double[] totalSample = generateSample();
        double[][] subSamples = generatePartition(totalSample);
        int nSamples = subSamples.length;

        // Compute combined stats directly
        StreamingStatistics totalStats = new StreamingStatistics();
        for (int i = 0; i < totalSample.length; i++) {
            totalStats.addValue(totalSample[i]);
        }

        // Now compute subsample stats individually and aggregate
        StreamingStatistics[] subSampleStats = new StreamingStatistics[nSamples];
        for (int i = 0; i < nSamples; i++) {
            subSampleStats[i] = new StreamingStatistics();
        }
        Collection<StatisticalSummary> aggregate = new ArrayList<StatisticalSummary>();
        for (int i = 0; i < nSamples; i++) {
            for (int j = 0; j < subSamples[i].length; j++) {
                subSampleStats[i].addValue(subSamples[i][j]);
            }
            aggregate.add(subSampleStats[i].getSummary());
        }

        // Compare values
        StatisticalSummary aggregatedStats = StatisticalSummary.aggregate(aggregate);
        assertStatisticalSummaryEquals(totalStats.getSummary(), aggregatedStats, 10E-12);
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
        Collection<StreamingStatistics> aggregate = new ArrayList<StreamingStatistics>();
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < subSamples[i].length; j++) {
                subSampleStats[i].addValue(subSamples[i][j]);
            }
            aggregate.add(subSampleStats[i]);
        }

        // Compare values
        StatisticalSummary aggregatedStats = StatisticalSummary.aggregate(aggregate);
        assertStatisticalSummaryEquals(totalStats.getSummary(), aggregatedStats, 10E-12);
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
        StreamingStatistics[] subSampleStats = new StreamingStatistics[5];
        for (int i = 0; i < 4; i++) {
            subSampleStats[i] = new StreamingStatistics();
        }
        Collection<StreamingStatistics> aggregate = new ArrayList<StreamingStatistics>();
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < subSamples[i].length; j++) {
                subSampleStats[i].addValue(subSamples[i][j]);
            }
            aggregate.add(subSampleStats[i]);
        }

        // Compare values
        StatisticalSummary aggregatedStats = StatisticalSummary.aggregate(aggregate);
        assertStatisticalSummaryEquals(totalStats.getSummary(), aggregatedStats, 10E-12);
    }

    /**
     * Verifies that a StatisticalSummary and a StatisticalSummaryValues are equal up
     * to delta, with NaNs, infinities returned in the same spots. For max, min, n, values
     * have to agree exactly, delta is used only for sum, mean, variance, std dev.
     */
    protected static void assertStatisticalSummaryEquals(StatisticalSummary expected,
                                                         StatisticalSummary observed,
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
        final int sampleSize = randomDataGenerator.nextInt(10,  100);
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
        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(100);
        final int length = sample.length;
        final double[][] out = new double[5][];
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
