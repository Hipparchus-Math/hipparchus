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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.random;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.hipparchus.RetryRunner;
import org.hipparchus.UnitTestUtils;
import org.hipparchus.distribution.continuous.BetaDistribution;
import org.hipparchus.distribution.continuous.EnumeratedRealDistribution;
import org.hipparchus.distribution.continuous.ExponentialDistribution;
import org.hipparchus.distribution.continuous.GammaDistribution;
import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.distribution.discrete.EnumeratedIntegerDistribution;
import org.hipparchus.distribution.discrete.PoissonDistribution;
import org.hipparchus.distribution.discrete.ZipfDistribution;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;

/**
 * Test cases for the RandomDataGenerator class.
 *
 */
@RunWith(RetryRunner.class)
public class RandomDataGeneratorTest {

    public RandomDataGeneratorTest() {
        randomData = RandomDataGenerator.of(new Well19937c());
        randomData.setSeed(100);
    }

    protected final long smallSampleSize = 1000;
    protected final double[] expected = { 250, 250, 250, 250 };
    protected final int largeSampleSize = 10000;
    private final String[] hex = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
            "a", "b", "c", "d", "e", "f" };
    protected RandomDataGenerator randomData = null;

    @Test
    public void testNextIntExtremeValues() {
        int x = randomData.nextInt(Integer.MIN_VALUE, Integer.MAX_VALUE);
        int y = randomData.nextInt(Integer.MIN_VALUE, Integer.MAX_VALUE);
        Assert.assertFalse(x == y);
    }

    @Test
    public void testNextLongExtremeValues() {
        long x = randomData.nextLong(Long.MIN_VALUE, Long.MAX_VALUE);
        long y = randomData.nextLong(Long.MIN_VALUE, Long.MAX_VALUE);
        Assert.assertFalse(x == y);
    }

    @Test
    public void testNextUniformExtremeValues() {
        double x = randomData.nextUniform(-Double.MAX_VALUE, Double.MAX_VALUE);
        double y = randomData.nextUniform(-Double.MAX_VALUE, Double.MAX_VALUE);
        Assert.assertFalse(x == y);
        Assert.assertFalse(Double.isNaN(x));
        Assert.assertFalse(Double.isNaN(y));
        Assert.assertFalse(Double.isInfinite(x));
        Assert.assertFalse(Double.isInfinite(y));
    }

    @Test
    public void testNextIntIAE() {
        try {
            randomData.nextInt(4, 3);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testNextIntNegativeToPositiveRange() {
        for (int i = 0; i < 5; i++) {
            checkNextIntUniform(-3, 5);
            checkNextIntUniform(-3, 6);
        }
    }

    @Test
    public void testNextIntNegativeRange() {
        for (int i = 0; i < 5; i++) {
            checkNextIntUniform(-7, -4);
            checkNextIntUniform(-15, -2);
            checkNextIntUniform(Integer.MIN_VALUE + 1, Integer.MIN_VALUE + 12);
        }
    }

    @Test
    public void testNextIntPositiveRange() {
        for (int i = 0; i < 5; i++) {
            checkNextIntUniform(0, 3);
            checkNextIntUniform(2, 12);
            checkNextIntUniform(1,2);
            checkNextIntUniform(Integer.MAX_VALUE - 12, Integer.MAX_VALUE - 1);
        }
    }

    private void checkNextIntUniform(int min, int max) {
        final int len = max - min + 1;
        final UnitTestUtils.Frequency<Integer> freq = new UnitTestUtils.Frequency<Integer>();
        for (int i = 0; i < smallSampleSize; i++) {
            final int value = randomData.nextInt(min, max);
            Assert.assertTrue("nextInt range", (value >= min) && (value <= max));
            freq.addValue(value);
        }
        final long[] observed = new long[len];
        for (int i = 0; i < len; i++) {
            observed[i] = freq.getCount(min + i);
        }
        final double[] expected = new double[len];
        for (int i = 0; i < len; i++) {
            expected[i] = 1d / len;
        }

        UnitTestUtils.assertChiSquareAccept(expected, observed, 0.001);
    }

    @Test
    public void testNextIntWideRange() {
        int lower = -0x6543210F;
        int upper =  0x456789AB;
        int max   = Integer.MIN_VALUE;
        int min   = Integer.MAX_VALUE;
        for (int i = 0; i < 1000000; ++i) {
            int r = randomData.nextInt(lower, upper);
            max = FastMath.max(max, r);
            min = FastMath.min(min, r);
            Assert.assertTrue(r >= lower);
            Assert.assertTrue(r <= upper);
        }
        double ratio = (((double) max)   - ((double) min)) /
                       (((double) upper) - ((double) lower));
        Assert.assertTrue(ratio > 0.99999);
    }

    @Test
    public void testNextLongIAE() {
        try {
            randomData.nextLong(4, 3);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testNextLongNegativeToPositiveRange() {
        for (int i = 0; i < 5; i++) {
            checkNextLongUniform(-3, 5);
            checkNextLongUniform(-3, 6);
        }
    }

    @Test
    public void testNextLongNegativeRange() {
        for (int i = 0; i < 5; i++) {
            checkNextLongUniform(-7, -4);
            checkNextLongUniform(-15, -2);
            checkNextLongUniform(Long.MIN_VALUE + 1, Long.MIN_VALUE + 12);
        }
    }

    @Test
    public void testNextLongPositiveRange() {
        for (int i = 0; i < 5; i++) {
            checkNextLongUniform(0, 3);
            checkNextLongUniform(2, 12);
            checkNextLongUniform(Long.MAX_VALUE - 12, Long.MAX_VALUE - 1);
        }
    }

    private void checkNextLongUniform(long min, long max) {
        final int len = ((int) (max - min)) + 1;
        final UnitTestUtils.Frequency<Integer> freq = new UnitTestUtils.Frequency<Integer>();
        for (int i = 0; i < smallSampleSize; i++) {
            final long value = randomData.nextLong(min, max);
            Assert.assertTrue("nextLong range: " + value + " " + min + " " + max,
                              (value >= min) && (value <= max));
            freq.addValue((int)value);
        }
        final long[] observed = new long[len];
        for (int i = 0; i < len; i++) {
            observed[i] = freq.getCount((int) min + i);
        }
        final double[] expected = new double[len];
        for (int i = 0; i < len; i++) {
            expected[i] = 1d / len;
        }

        UnitTestUtils.assertChiSquareAccept(expected, observed, 0.01);
    }

    @Test
    public void testNextLongWideRange() {
        long lower = -0x6543210FEDCBA987L;
        long upper =  0x456789ABCDEF0123L;
        long max = Long.MIN_VALUE;
        long min = Long.MAX_VALUE;
        for (int i = 0; i < 10000000; ++i) {
            long r = randomData.nextLong(lower, upper);
            max = FastMath.max(max, r);
            min = FastMath.min(min, r);
            Assert.assertTrue(r >= lower);
            Assert.assertTrue(r <= upper);
        }
        double ratio = (((double) max)   - ((double) min)) /
                       (((double) upper) - ((double) lower));
        Assert.assertTrue(ratio > 0.99999);
    }

    /**
     * Make sure that empirical distribution of random Poisson(4)'s has P(X &lt;= 5)
     * close to actual cumulative Poisson probability and that nextPoisson
     * fails when mean is non-positive.
     */
    @Test
    public void testNextPoisson() {
        try {
            randomData.nextPoisson(0);
            Assert.fail("zero mean -- expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            randomData.nextPoisson(-1);
            Assert.fail("negative mean supplied -- MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            randomData.nextPoisson(0);
            Assert.fail("0 mean supplied -- MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }

        final double mean = 4.0d;
        final int len = 5;
        PoissonDistribution poissonDistribution = new PoissonDistribution(mean);
        final UnitTestUtils.Frequency<Integer> freq = new UnitTestUtils.Frequency<Integer>();
        randomData.setSeed(1000);
        for (int i = 0; i < largeSampleSize; i++) {
            freq.addValue(randomData.nextPoisson(mean));
        }
        final long[] observed = new long[len];
        for (int i = 0; i < len; i++) {
            observed[i] = freq.getCount(i + 1);
        }
        final double[] expected = new double[len];
        for (int i = 0; i < len; i++) {
            expected[i] = poissonDistribution.probability(i + 1) * largeSampleSize;
        }

        UnitTestUtils.assertChiSquareAccept(expected, observed, 0.0001);
    }

    @Test
    public void testNextPoissonConsistency() {

        // Small integral means
        for (int i = 1; i < 100; i++) {
            checkNextPoissonConsistency(i);
        }
        // non-integer means
        for (int i = 1; i < 10; i++) {
            checkNextPoissonConsistency(randomData.nextUniform(1, 1000));
        }
        // large means
        for (int i = 1; i < 10; i++) {
            checkNextPoissonConsistency(randomData.nextUniform(1000, 10000));
        }
    }

    /**
     * Verifies that nextPoisson(mean) generates an empirical distribution of values
     * consistent with PoissonDistributionImpl by generating 1000 values, computing a
     * grouped frequency distribution of the observed values and comparing this distribution
     * to the corresponding expected distribution computed using PoissonDistributionImpl.
     * Uses ChiSquare test of goodness of fit to evaluate the null hypothesis that the
     * distributions are the same. If the null hypothesis can be rejected with confidence
     * 1 - alpha, the check fails.
     */
    public void checkNextPoissonConsistency(double mean) {
        // Generate sample values
        final int sampleSize = 1000;        // Number of deviates to generate
        final int minExpectedCount = 7;     // Minimum size of expected bin count
        long maxObservedValue = 0;
        final double alpha = 0.001;         // Probability of false failure
        UnitTestUtils.Frequency<Long> frequency = new UnitTestUtils.Frequency<Long>();
        for (int i = 0; i < sampleSize; i++) {
            long value = randomData.nextPoisson(mean);
            if (value > maxObservedValue) {
                maxObservedValue = value;
            }
            frequency.addValue(value);
        }

        /*
         *  Set up bins for chi-square test.
         *  Ensure expected counts are all at least minExpectedCount.
         *  Start with upper and lower tail bins.
         *  Lower bin = [0, lower); Upper bin = [upper, +inf).
         */
        PoissonDistribution poissonDistribution = new PoissonDistribution(mean);
        int lower = 1;
        while (poissonDistribution.cumulativeProbability(lower - 1) * sampleSize < minExpectedCount) {
            lower++;
        }
        int upper = (int) (5 * mean);  // Even for mean = 1, not much mass beyond 5
        while ((1 - poissonDistribution.cumulativeProbability(upper - 1)) * sampleSize < minExpectedCount) {
            upper--;
        }

        // Set bin width for interior bins.  For poisson, only need to look at end bins.
        int binWidth = 0;
        boolean widthSufficient = false;
        double lowerBinMass = 0;
        double upperBinMass = 0;
        while (!widthSufficient) {
            binWidth++;
            lowerBinMass = poissonDistribution.probability(lower - 1, lower + binWidth - 1);
            upperBinMass = poissonDistribution.probability(upper - binWidth - 1, upper - 1);
            widthSufficient = FastMath.min(lowerBinMass, upperBinMass) * sampleSize >= minExpectedCount;
        }

        /*
         *  Determine interior bin bounds.  Bins are
         *  [1, lower = binBounds[0]), [lower, binBounds[1]), [binBounds[1], binBounds[2]), ... ,
         *    [binBounds[binCount - 2], upper = binBounds[binCount - 1]), [upper, +inf)
         *
         */
        List<Integer> binBounds = new ArrayList<Integer>();
        binBounds.add(lower);
        int bound = lower + binWidth;
        while (bound < upper - binWidth) {
            binBounds.add(bound);
            bound += binWidth;
        }
        binBounds.add(upper); // The size of bin [binBounds[binCount - 2], upper) satisfies binWidth <= size < 2*binWidth.

        // Compute observed and expected bin counts
        final int binCount = binBounds.size() + 1;
        long[] observed = new long[binCount];
        double[] expected = new double[binCount];

        // Bottom bin
        observed[0] = 0;
        for (int i = 0; i < lower; i++) {
            observed[0] += frequency.getCount((long)i);
        }
        expected[0] = poissonDistribution.cumulativeProbability(lower - 1) * sampleSize;

        // Top bin
        observed[binCount - 1] = 0;
        for (int i = upper; i <= maxObservedValue; i++) {
            observed[binCount - 1] += frequency.getCount((long)i);
        }
        expected[binCount - 1] = (1 - poissonDistribution.cumulativeProbability(upper - 1)) * sampleSize;

        // Interior bins
        for (int i = 1; i < binCount - 1; i++) {
            observed[i] = 0;
            for (int j = binBounds.get(i - 1); j < binBounds.get(i); j++) {
                observed[i] += frequency.getCount((long)j);
            } // Expected count is (mass in [binBounds[i-1], binBounds[i])) * sampleSize
            expected[i] = (poissonDistribution.cumulativeProbability(binBounds.get(i) - 1) -
                poissonDistribution.cumulativeProbability(binBounds.get(i - 1) -1)) * sampleSize;
        }

        // Use chisquare test to verify that generated values are poisson(mean)-distributed
        // Fail if we can reject null hypothesis that distributions are the same
        if (UnitTestUtils.chiSquareTest(expected, observed) <  alpha) {
            StringBuilder msgBuffer = new StringBuilder();
            DecimalFormat df = new DecimalFormat("#.##");
            msgBuffer.append("Chisquare test failed for mean = ");
            msgBuffer.append(mean);
            msgBuffer.append(" p-value = ");
            msgBuffer.append(UnitTestUtils.chiSquareTest(expected, observed));
            msgBuffer.append(" chisquare statistic = ");
            msgBuffer.append(UnitTestUtils.chiSquare(expected, observed));
            msgBuffer.append(". \n");
            msgBuffer.append("bin\t\texpected\tobserved\n");
            for (int i = 0; i < expected.length; i++) {
                msgBuffer.append("[");
                msgBuffer.append(i == 0 ? 1: binBounds.get(i - 1));
                msgBuffer.append(",");
                msgBuffer.append(i == binBounds.size() ? "inf": binBounds.get(i));
                msgBuffer.append(")");
                msgBuffer.append("\t\t");
                msgBuffer.append(df.format(expected[i]));
                msgBuffer.append("\t\t");
                msgBuffer.append(observed[i]);
                msgBuffer.append("\n");
            }
            msgBuffer.append("This test can fail randomly due to sampling error with probability ");
            msgBuffer.append(alpha);
            msgBuffer.append(".");
            Assert.fail(msgBuffer.toString());
        }
    }

    /** test dispersion and failure modes for nextHex() */
    @Test
    public void testNextHex() {
        try {
            randomData.nextHexString(-1);
            Assert.fail("negative length supplied -- MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            randomData.nextHexString(0);
            Assert.fail("zero length supplied -- MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        String hexString = randomData.nextHexString(3);
        if (hexString.length() != 3) {
            Assert.fail("incorrect length for generated string");
        }
        hexString = randomData.nextHexString(1);
        if (hexString.length() != 1) {
            Assert.fail("incorrect length for generated string");
        }
        try {
            hexString = randomData.nextHexString(0);
            Assert.fail("zero length requested -- expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        UnitTestUtils.Frequency<String> f = new UnitTestUtils.Frequency<String> ();
        for (int i = 0; i < smallSampleSize; i++) {
            hexString = randomData.nextHexString(100);
            if (hexString.length() != 100) {
                Assert.fail("incorrect length for generated string");
            }
            for (int j = 0; j < hexString.length(); j++) {
                f.addValue(hexString.substring(j, j + 1));
            }
        }
        double[] expected = new double[16];
        long[] observed = new long[16];
        for (int i = 0; i < 16; i++) {
            expected[i] = (double) smallSampleSize * 100 / 16;
            observed[i] = f.getCount(hex[i]);
        }
        UnitTestUtils.assertChiSquareAccept(expected, observed, 0.001);
    }

    @Test
    public void testNextUniformIAE() {
        try {
            randomData.nextUniform(4, 3);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            randomData.nextUniform(0, Double.POSITIVE_INFINITY);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            randomData.nextUniform(Double.NEGATIVE_INFINITY, 0);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            randomData.nextUniform(0, Double.NaN);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            randomData.nextUniform(Double.NaN, 0);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testNextUniformUniformPositiveBounds() {
        for (int i = 0; i < 5; i++) {
            checkNextUniformUniform(0, 10);
        }
    }

    @Test
    public void testNextUniformUniformNegativeToPositiveBounds() {
        for (int i = 0; i < 5; i++) {
            checkNextUniformUniform(-3, 5);
        }
    }

    @Test
    public void testNextUniformUniformNegaiveBounds() {
        for (int i = 0; i < 5; i++) {
            checkNextUniformUniform(-7, -3);
        }
    }

    @Test
    public void testNextUniformUniformMaximalInterval() {
        for (int i = 0; i < 5; i++) {
            checkNextUniformUniform(-Double.MAX_VALUE, Double.MAX_VALUE);
        }
    }

    private void checkNextUniformUniform(double min, double max) {
        // Set up bin bounds - min, binBound[0], ..., binBound[binCount-2], max
        final int binCount = 5;
        final double binSize = max / binCount - min/binCount; // Prevent overflow in extreme value case
        final double[] binBounds = new double[binCount - 1];
        binBounds[0] = min + binSize;
        for (int i = 1; i < binCount - 1; i++) {
            binBounds[i] = binBounds[i - 1] + binSize;  // + instead of * to avoid overflow in extreme case
        }

        UnitTestUtils.Frequency<Integer> freq = new UnitTestUtils.Frequency<Integer>();
        for (int i = 0; i < smallSampleSize; i++) {
            final double value = randomData.nextUniform(min, max);
            Assert.assertTrue("nextUniform range", (value > min) && (value < max));
            // Find bin
            int j = 0;
            while (j < binCount - 1 && value > binBounds[j]) {
                j++;
            }
            freq.addValue(j);
        }

        final long[] observed = new long[binCount];
        for (int i = 0; i < binCount; i++) {
            observed[i] = freq.getCount(i);
        }
        final double[] expected = new double[binCount];
        for (int i = 0; i < binCount; i++) {
            expected[i] = 1d / binCount;
        }

        UnitTestUtils.assertChiSquareAccept(expected, observed, 0.01);
    }

    /** test exclusive endpoints of nextUniform **/
    @Test
    public void testNextUniformExclusiveEndpoints() {
        for (int i = 0; i < 1000; i++) {
            double u = randomData.nextUniform(0.99, 1);
            Assert.assertTrue(u > 0.99 && u < 1);
        }
    }

    /** test failure modes and distribution of nextGaussian() */
    @Test
    public void testNextGaussian() {
        try {
            randomData.nextNormal(0, 0);
            Assert.fail("zero sigma -- MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        double[] quartiles = UnitTestUtils.getDistributionQuartiles(new NormalDistribution(0,1));
        long[] counts = new long[4];
        randomData.setSeed(1000);
        for (int i = 0; i < 1000; i++) {
            double value = randomData.nextNormal(0, 1);
            UnitTestUtils.updateCounts(value, counts, quartiles);
        }
        UnitTestUtils.assertChiSquareAccept(expected, counts, 0.001);
    }

    /** test failure modes and distribution of nextExponential() */
    @Test
    public void testNextExponential() {
        try {
            randomData.nextExponential(-1);
            Assert.fail("negative mean -- expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        try {
            randomData.nextExponential(0);
            Assert.fail("zero mean -- expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
        double[] quartiles;
        long[] counts;

        // Mean 1
        quartiles = UnitTestUtils.getDistributionQuartiles(new ExponentialDistribution(1));
        counts = new long[4];
        randomData.setSeed(1000);
        for (int i = 0; i < 1000; i++) {
            double value = randomData.nextExponential(1);
            UnitTestUtils.updateCounts(value, counts, quartiles);
        }
        UnitTestUtils.assertChiSquareAccept(expected, counts, 0.001);

        // Mean 5
        quartiles = UnitTestUtils.getDistributionQuartiles(new ExponentialDistribution(5));
        counts = new long[4];
        randomData.setSeed(1000);
        for (int i = 0; i < 1000; i++) {
            double value = randomData.nextExponential(5);
            UnitTestUtils.updateCounts(value, counts, quartiles);
        }
        UnitTestUtils.assertChiSquareAccept(expected, counts, 0.001);
    }

    /** test reseeding, algorithm/provider games */
    @Test
    public void testConfig() {
        randomData.setSeed(1000);
        double v = randomData.nextUniform(0, 1);
        randomData.setSeed(System.currentTimeMillis());
        Assert.assertTrue("different seeds", FastMath.abs(v - randomData.nextUniform(0, 1)) > 10E-12);
        randomData.setSeed(1000);
        Assert.assertEquals("same seeds", v, randomData.nextUniform(0, 1), 10E-12);
    }

    /** tests for nextSample() sampling from Collection */
    @Test
    public void testNextSample() {
        Object[][] c = { { "0", "1" }, { "0", "2" }, { "0", "3" },
                { "0", "4" }, { "1", "2" }, { "1", "3" }, { "1", "4" },
                { "2", "3" }, { "2", "4" }, { "3", "4" } };
        long[] observed = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expected = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 100 };

        HashSet<Object> cPop = new HashSet<Object>(); // {0,1,2,3,4}
        for (int i = 0; i < 5; i++) {
            cPop.add(Integer.toString(i));
        }

        Object[] sets = new Object[10]; // 2-sets from 5
        for (int i = 0; i < 10; i++) {
            HashSet<Object> hs = new HashSet<Object>();
            hs.add(c[i][0]);
            hs.add(c[i][1]);
            sets[i] = hs;
        }

        for (int i = 0; i < 1000; i++) {
            Object[] cSamp = randomData.nextSample(cPop, 2);
            observed[findSample(sets, cSamp)]++;
        }

        /*
         * Use ChiSquare dist with df = 10-1 = 9, alpha = .001 Change to 21.67
         * for alpha = .01
         */
        Assert.assertTrue("chi-square test -- will fail about 1 in 1000 times",
                UnitTestUtils.chiSquare(expected, observed) < 27.88);

        // Make sure sample of size = size of collection returns same collection
        HashSet<Object> hs = new HashSet<Object>();
        hs.add("one");
        Object[] one = randomData.nextSample(hs, 1);
        String oneString = (String) one[0];
        if ((one.length != 1) || !oneString.equals("one")) {
            Assert.fail("bad sample for set size = 1, sample size = 1");
        }

        // Make sure we fail for sample size > collection size
        try {
            one = randomData.nextSample(hs, 2);
            Assert.fail("sample size > set size, expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }

        // Make sure we fail for empty collection
        try {
            hs = new HashSet<Object>();
            one = randomData.nextSample(hs, 0);
            Assert.fail("n = k = 0, expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @SuppressWarnings("unchecked")
    private int findSample(Object[] u, Object[] samp) {
        for (int i = 0; i < u.length; i++) {
            HashSet<Object> set = (HashSet<Object>) u[i];
            HashSet<Object> sampSet = new HashSet<Object>();
            for (int j = 0; j < samp.length; j++) {
                sampSet.add(samp[j]);
            }
            if (set.equals(sampSet)) {
                return i;
            }
        }
        Assert.fail("sample not found:{" + samp[0] + "," + samp[1] + "}");
        return -1;
    }

    /** tests for nextPermutation */
    @Test
    public void testNextPermutation() {
        int[][] p = { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 },
                { 2, 0, 1 }, { 2, 1, 0 } };
        long[] observed = { 0, 0, 0, 0, 0, 0 };
        double[] expected = { 100, 100, 100, 100, 100, 100 };

        for (int i = 0; i < 600; i++) {
            int[] perm = randomData.nextPermutation(3, 3);
            observed[findPerm(p, perm)]++;
        }

        String[] labels = {"{0, 1, 2}", "{ 0, 2, 1 }", "{ 1, 0, 2 }",
                "{ 1, 2, 0 }", "{ 2, 0, 1 }", "{ 2, 1, 0 }"};
        UnitTestUtils.assertChiSquareAccept(labels, expected, observed, 0.001);

        // Check size = 1 boundary case
        int[] perm = randomData.nextPermutation(1, 1);
        if ((perm.length != 1) || (perm[0] != 0)) {
            Assert.fail("bad permutation for n = 1, sample k = 1");

            // Make sure we fail for k size > n
            try {
                perm = randomData.nextPermutation(2, 3);
                Assert.fail("permutation k > n, expecting MathIllegalArgumentException");
            } catch (MathIllegalArgumentException ex) {
                // ignored
            }

            // Make sure we fail for n = 0
            try {
                perm = randomData.nextPermutation(0, 0);
                Assert.fail("permutation k = n = 0, expecting MathIllegalArgumentException");
            } catch (MathIllegalArgumentException ex) {
                // ignored
            }

            // Make sure we fail for k < n < 0
            try {
                perm = randomData.nextPermutation(-1, -3);
                Assert.fail("permutation k < n < 0, expecting MathIllegalArgumentException");
            } catch (MathIllegalArgumentException ex) {
                // ignored
            }

        }
    }

    private int findPerm(int[][] p, int[] samp) {
        for (int i = 0; i < p.length; i++) {
            boolean good = true;
            for (int j = 0; j < samp.length; j++) {
                if (samp[j] != p[i][j]) {
                    good = false;
                }
            }
            if (good) {
                return i;
            }
        }
        Assert.fail("permutation not found");
        return -1;
    }

    @Test
    public void testNextBeta() {
        double[] quartiles = UnitTestUtils.getDistributionQuartiles(new BetaDistribution(2,5));
        long[] counts = new long[4];
        randomData.setSeed(1000);
        for (int i = 0; i < 1000; i++) {
            double value = randomData.nextBeta(2, 5);
            UnitTestUtils.updateCounts(value, counts, quartiles);
        }
        UnitTestUtils.assertChiSquareAccept(expected, counts, 0.001);
    }

    @Test
    public void testNextGamma() {
        double[] quartiles;
        long[] counts;

        // Tests shape > 1, one case in the rejection sampling
        quartiles = UnitTestUtils.getDistributionQuartiles(new GammaDistribution(4, 2));
        counts = new long[4];
        randomData.setSeed(1000);
        for (int i = 0; i < 1000; i++) {
            double value = randomData.nextGamma(4, 2);
            UnitTestUtils.updateCounts(value, counts, quartiles);
        }
        UnitTestUtils.assertChiSquareAccept(expected, counts, 0.001);

        // Tests shape <= 1, another case in the rejection sampling
        quartiles = UnitTestUtils.getDistributionQuartiles(new GammaDistribution(0.3, 3));
        counts = new long[4];
        randomData.setSeed(1000);
        for (int i = 0; i < 1000; i++) {
            double value = randomData.nextGamma(0.3, 3);
            UnitTestUtils.updateCounts(value, counts, quartiles);
        }
        UnitTestUtils.assertChiSquareAccept(expected, counts, 0.001);
    }

    @Test
    public void testNextGamma2() {
        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(1000);
        final int sampleSize = 1000;
        final double alpha = 0.001;
        GammaDistribution dist = new GammaDistribution(3, 1);
        double[] values = randomDataGenerator.nextDeviates(dist, sampleSize);
        UnitTestUtils.assertGTest(dist, values, alpha);
        dist = new GammaDistribution(.4, 2);
        values = randomDataGenerator.nextDeviates(dist, sampleSize);
        UnitTestUtils.assertGTest(dist, values, alpha);
    }

    @Test
    public void testNextBeta2() {
        final double[] alphaBetas = {0.1, 1, 10, 100, 1000};
        final RandomGenerator random = new Well1024a(0x7829862c82fec2dal);
        final RandomDataGenerator randomDataGenerator = RandomDataGenerator.of(random);
        final int sampleSize = 1000;
        final double alphaCrit = 0.001;
        for (final double alpha : alphaBetas) {
            for (final double beta : alphaBetas) {
                final BetaDistribution betaDistribution = new BetaDistribution(alpha, beta);
                final double[] values = randomDataGenerator.nextDeviates(betaDistribution, sampleSize);
                UnitTestUtils.assertGTest(betaDistribution, values, alphaCrit);
            }
        }
    }

    @Test
    public void testNextZipf() {
        int sampleSize = 1000;

        int[] numPointsValues = {
            2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100
        };
        double[] exponentValues = {
            1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 2e-1, 5e-1,
            1. - 1e-9, 1.0, 1. + 1e-9, 1.1, 1.2, 1.3, 1.5, 1.6, 1.7, 1.8, 2.0,
            2.5, 3.0, 4., 5., 6., 7., 8., 9., 10., 20., 30., 100., 150.
        };

        for (int numPoints : numPointsValues) {
            for (double exponent : exponentValues) {
                double weightSum = 0.;
                double[] weights = new double[numPoints];
                for (int i = numPoints; i>=1; i-=1) {
                    weights[i-1] = Math.pow(i, -exponent);
                    weightSum += weights[i-1];
                }
                // use fixed seed, the test is expected to fail for more than 50% of all seeds because each test case can fail
                // with probability 0.001, the chance that all test cases do not fail is 0.999^(32*22) = 0.49442874426
                ZipfDistribution distribution = new ZipfDistribution(numPoints, exponent);
                randomData.setSeed(1001);
                double[] expectedCounts = new double[numPoints];
                long[] observedCounts = new long[numPoints];
                for (int i = 0; i < numPoints; i++) {
                    expectedCounts[i] = sampleSize * (weights[i]/weightSum);
                }
                int[] sample = randomData.nextDeviates(distribution,sampleSize);
                for (int s : sample) {
                    observedCounts[s-1]++;
                }
                UnitTestUtils.assertChiSquareAccept(expectedCounts, observedCounts, 0.001);
            }
        }
    }

    @Test
    public void testNextSampleWithReplacement() {
        final int sampleSize = 1000;
        final double[] weights = {1, 2, 3, 4};
        final int[] sample = randomData.nextSampleWithReplacement(sampleSize, weights);
        final double[] expected = {sampleSize/10d, sampleSize/5d, 3*sampleSize/10d, 2*sampleSize/5d};
        final long[] observed = {0, 0, 0, 0};
        for (int i = 0; i < sampleSize; i++) {
            observed[sample[i]]++;
        }
        UnitTestUtils.assertChiSquareAccept(new String[] {"0", "1", "2","3"}, expected, observed, 0.01);
    }

    @Test
    public void testNextSampleWithReplacementPointMass() {
        final int sampleSize = 2;
        double[] weights = {1};
        final int[] expected = new int[] {0, 0};
        UnitTestUtils.assertEquals(expected, randomData.nextSampleWithReplacement(sampleSize, weights));
        weights = new double[] {1, 0};
        UnitTestUtils.assertEquals(expected, randomData.nextSampleWithReplacement(sampleSize, weights));
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testNextSampleWithReplacementAllZeroWeights() {
        final double[] weights = {0, 0, 0};
        randomData.nextSampleWithReplacement(1, weights);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testNextSampleWithReplacementNegativeWeights() {
        final double[] weights = {-1, 1, 0};
        randomData.nextSampleWithReplacement(1, weights);
    }

    @Test
    public void testNextSampleWithReplacement0SampleSize() {
        final double[] weights = {1, 0};
        final int[] expected = {};
        UnitTestUtils.assertEquals(expected, randomData.nextSampleWithReplacement(0, weights));
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testNextSampleWithReplacementNegativeSampleSize() {
        final double[] weights = {1, 0};
        randomData.nextSampleWithReplacement(-1, weights);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testNextSampleWithReplacementNaNWeights() {
        final double[] weights = {1, Double.NaN};
        randomData.nextSampleWithReplacement(0, weights);
    }

    @Test
    public void testNextDeviateEnumeratedIntegerDistribution() {
        final int sampleSize = 1000;
        final int[] data = new int[] {0, 1, 1, 2, 2, 2};
        final EnumeratedIntegerDistribution dist = new EnumeratedIntegerDistribution(data);
        final int[] sample = randomData.nextDeviates(dist, sampleSize);
        final double[] expected = {sampleSize/6d, sampleSize/3d, sampleSize/2d};
        final long[] observed = {0, 0, 0};
        for (int i = 0; i < sampleSize; i++) {
            observed[sample[i]]++;
        }
        UnitTestUtils.assertChiSquareAccept(new String[] {"0", "1", "2"}, expected, observed, 0.01);
    }

    @Test
    public void testNextDeviateEnumeratedRealDistribution() {
        final int sampleSize = 1000;
        final double[] data = new double[] {0, 1, 1, 2, 2, 2};
        final EnumeratedRealDistribution dist = new EnumeratedRealDistribution(data);
        final double[] sample = randomData.nextDeviates(dist, sampleSize);
        final double[] expected = {sampleSize/6d, sampleSize/3d, sampleSize/2d};
        final long[] observed = {0, 0, 0};
        for (int i = 0; i < sampleSize; i++) {
            observed[(int)sample[i]]++;
        }
        UnitTestUtils.assertChiSquareAccept(new String[] {"0", "1", "2"}, expected, observed, 0.01);
    }

}
