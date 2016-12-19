/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.stat.descriptive.rank;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.distribution.continuous.ExponentialDistribution;
import org.hipparchus.distribution.continuous.GammaDistribution;
import org.hipparchus.distribution.continuous.LogNormalDistribution;
import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.MersenneTwister;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well19937c;
import org.hipparchus.stat.descriptive.StorelessUnivariateStatistic;
import org.hipparchus.stat.descriptive.StorelessUnivariateStatisticAbstractTest;
import org.hipparchus.util.FastMath;
import org.junit.Test;

/**
 * Test cases for the {@link PSquarePercentile} class which naturally extends
 * {@link StorelessUnivariateStatisticAbstractTest}.
 */
public class RandomPercentileTest extends
        StorelessUnivariateStatisticAbstractTest {

    protected double tolerance = 10E-12;

    private final RandomGenerator randomGenerator = new Well19937c(1000);

    @Override
    public RandomPercentile getUnivariateStatistic() {
        return new RandomPercentile();  // Median with default epsilon and PRNG
    }

    @Override
    public double expectedValue() {
        return this.median;
    }

    /**
     * Verifies that copied statistics remain equal to originals when
     * incremented the same way by making the copy after a majority of elements
     * are incremented
     */
    @Test
    public void testCopyConsistencyWithInitialMostElements() {

        StorelessUnivariateStatistic master = getUnivariateStatistic();
        StorelessUnivariateStatistic replica = null;

        // select a portion of testArray till 75 % of the length to load first
        long index = FastMath.round(0.75 * testArray.length);

        // Put first half in master and copy master to replica
        master.incrementAll(testArray, 0, (int) index);
        replica = master.copy();

        // Check same
        assertTrue(replica.equals(master));
        assertTrue(master.equals(replica));

        // Now add second part to both and check again
        master.incrementAll(testArray, (int) index, (int) (testArray.length - index));
        replica.incrementAll(testArray, (int) index, (int) (testArray.length - index));
        assertTrue(replica.equals(master));
        assertTrue(master.equals(replica));
    }

    /**
     * Verifies that copied statistics remain equal to originals when
     * incremented the same way by way of copying original after just a few
     * elements are incremented
     */
    @Test
    public void testCopyConsistencyWithInitialFirstFewElements() {

        StorelessUnivariateStatistic master = getUnivariateStatistic();
        StorelessUnivariateStatistic replica = null;

        // select a portion of testArray which is 10% of the length to load
        // first
        long index = FastMath.round(0.1 * testArray.length);

        // Put first half in master and copy master to replica
        master.incrementAll(testArray, 0, (int) index);
        replica = master.copy();

        // Check same
        assertTrue(replica.equals(master));
        assertTrue(master.equals(replica));
        // Now add second part to both and check again
        master.incrementAll(testArray, (int) index, (int) (testArray.length - index));
        replica.incrementAll(testArray, (int) index, (int) (testArray.length - index));
        assertTrue(master.equals(master));
        assertTrue(replica.equals(replica));
        assertTrue(replica.equals(master));
        assertTrue(master.equals(replica));
    }


    //@Test
    public void testPSquaredEqualsAndMin() {
        PSquarePercentile ptile = new PSquarePercentile(0);
        assertEquals(ptile, ptile);
        assertFalse(ptile.equals(null));
        assertFalse(ptile.equals(new String()));
        // Just to check if there is no data get result for zeroth and 100th
        // ptile returns NAN
        assertTrue(Double.isNaN(ptile.getResult()));
        assertTrue(Double.isNaN(new PSquarePercentile(100).getResult()));

        double[] d = new double[] { 1, 3, 2, 4, 9, 10, 11 };
        ptile.incrementAll(d);
        assertEquals(ptile, ptile);
        assertEquals(1d, ptile.getResult(), 1e-02);// this calls min
    }

    // @Test
    public void testString() {
        PSquarePercentile ptile = new PSquarePercentile(95);
        assertNotNull(ptile.toString());
        ptile.increment(1);
        ptile.increment(2);
        ptile.increment(3);
        assertNotNull(ptile.toString());
        assertEquals(expectedValue(), ptile.evaluate(testArray), getTolerance());
        assertNotNull(ptile.toString());
    }

    @Test
    public void testPercentileSmallSample() {
        double[] d = new double[] { 1, 3, 2, 4 };
        final RandomPercentile randomPercentile = new RandomPercentile();
        randomPercentile.incrementAll(d);
        Percentile p = new Percentile(30d);
        assertEquals(p.evaluate(d), randomPercentile.getResult(30d), 1.0e-5);
        p = new Percentile(25);
        assertEquals(p.evaluate(d), randomPercentile.getResult(25d), 1.0e-5);
        p = new Percentile(75);
        assertEquals(p.evaluate(d),randomPercentile.getResult(75d), 1.0e-5);
        p = new Percentile(50);
        assertEquals(p.evaluate(d),randomPercentile.getResult(50d), 1.0e-5);
    }

    @Test(expected = MathIllegalArgumentException.class)
    public void testNonPositiveEpsilon() {
        double[] d =
                new double[] { 95.1772, 95.1567, 95.1937, 95.1959, 95.1442,
                        95.0610, 95.1591, 95.1195, 95.1772, 95.0925, 95.1990,
                        95.1682 };
        RandomPercentile p = new RandomPercentile(0);
        p.evaluate(d, 0, d.length);
    }

    @Test
    public void testNISTExample() {
        double[] d =
                new double[] { 95.1772, 95.1567, 95.1937, 95.1959, 95.1442,
                        95.0610, 95.1591, 95.1195, 95.1772, 95.0925, 95.1990,
                        95.1682 };
        assertEquals(95.1981, new RandomPercentile().evaluate(90d, d), 1.0e-4);
        assertEquals(95.061, new RandomPercentile().evaluate(0d, d), 0);
        assertEquals(95.1990, new RandomPercentile().evaluate(100d, d, 0, d.length), 0);
    }

    @Test
    public void test5() {
        RandomPercentile percentile = new RandomPercentile();
        assertEquals(this.percentile5, percentile.evaluate(5, testArray), 0.0001);
    }

    @Test(expected = NullArgumentException.class)
    public void testNull() {
        PSquarePercentile percentile = new PSquarePercentile(50d);
        double[] nullArray = null;
        percentile.evaluate(nullArray);
    }

    @Test
    public void testEmpty() {
        PSquarePercentile percentile = new PSquarePercentile(50d);
        double[] emptyArray = new double[] {};
        assertTrue(Double.isNaN(percentile.evaluate(emptyArray)));
    }

    @Test
    public void testSingleton() {
        RandomPercentile percentile = new RandomPercentile();
        double[] singletonArray = new double[] { 1d };
        assertEquals(1d, percentile.evaluate(singletonArray), 0);
        assertEquals(1d, percentile.evaluate(singletonArray, 0, 1), 0);
        percentile = new RandomPercentile();
        assertEquals(1d, percentile.evaluate(5, singletonArray, 0, 1), 0);
        percentile = new RandomPercentile();
        assertEquals(1d, percentile.evaluate(100, singletonArray, 0, 1), 0);
        percentile = new RandomPercentile();
        assertTrue(Double.isNaN(percentile.evaluate(100, singletonArray, 0, 0)));
    }

    @Test
    public void testSpecialValues() {
        RandomPercentile percentile = new RandomPercentile();
        double[] specialValues = new double[] { 0d, 1d, 2d, 3d, 4d, Double.NaN };
        assertEquals(2d, percentile.evaluate(specialValues), 0);
        specialValues =
            new double[] { Double.NEGATIVE_INFINITY, 1d, 2d, 3d, Double.NaN, Double.POSITIVE_INFINITY };
        assertEquals(2d, percentile.evaluate(specialValues), 0);
        specialValues =
            new double[] { 1d, 1d, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY };
        assertTrue(Double.isInfinite(percentile.evaluate(specialValues)));
        specialValues = new double[] { 1d, 1d, Double.NaN, Double.NaN };
        assertFalse(Double.isNaN(percentile.evaluate(specialValues)));
        specialValues =
            new double[] { 1d, 1d, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY };
        percentile = new RandomPercentile();
        assertTrue(Double.isNaN(percentile.evaluate(specialValues)));
    }

    @Test
    public void testArrayExample() {
        assertEquals(this.percentile95, new RandomPercentile().evaluate(95d,testArray), getTolerance());
    }

    @Test
    public void testReduceSmallDataSet() {
        final RandomDataGenerator random = new RandomDataGenerator(1000);
        final long n = 1000;
        double[] combined = new double[10000];
        int i = 0;
        final List<RandomPercentile> aggregates = new ArrayList<RandomPercentile>();
        for (int j = 0; j < 10; j++) {
            final RandomPercentile randomPercentile = new RandomPercentile();
            for (int k = 0; k < n; k++) {
                final double value = random.nextGaussian();
                randomPercentile.accept(value);
                combined[i++] = value;
            }
            aggregates.add(randomPercentile);
        }
        // Check some quantiles
        final Percentile master = new Percentile();
        final RandomPercentile randomMaster = new RandomPercentile();
        for (int l = 0; l < 5; l++) {
            final double percentile = l * 15 + 1;
            assertEquals(master.evaluate(combined, percentile),
                    randomMaster.reduce(percentile, aggregates), Double.MIN_VALUE);
        }
    }

    @Test
    public void testReduceLargeDataSet() {
        final RandomDataGenerator random = new RandomDataGenerator(1000);
        final long n = 1000000;
        final RandomGenerator randomGenerator = new RandomDataGenerator(1000);
        final RandomPercentile randomMaster = new RandomPercentile(randomGenerator);
        final PSquarePercentile pSquare = new PSquarePercentile(1);
        final List<RandomPercentile> aggregates = new ArrayList<RandomPercentile>();
        for (int j = 0; j < 5; j++) {
            final RandomPercentile randomPercentile = new RandomPercentile(randomGenerator);
            for (int k = 0; k < n; k++) {
                final double value = random.nextGaussian();
                randomPercentile.accept(value);
                randomMaster.accept(value);
                pSquare.increment(value);
            }
            aggregates.add(randomPercentile);
        }
        // Check some quantiles
        for (int l = 0; l < 5; l++) {
            final double percentile = l * 13 + 1;
            assertEquals("percentile = " + percentile, randomMaster.getResult(percentile),
                    randomMaster.reduce(percentile, aggregates), 5E-3);
        }
    }

    private Double[] randomTestData(int factor, int values) {
        Double[] test = new Double[values];
        for (int i = 0; i < test.length; i++) {
            test[i] = Math.abs(randomGenerator.nextDouble() * factor);
        }
        return test;
    }

    @Test
    public void testAccept() {
        final RandomPercentile randomPercentile = new RandomPercentile();
        assertTrue(Double.isNaN(randomPercentile.getResult()));
        Double[] test = randomTestData(100, 10000);

        for (Double value : test) {
            randomPercentile.increment(value);
            assertTrue(randomPercentile.getResult() >= 0);
        }
    }

    private void assertValues(Double a, Double b, double delta) {
        if (Double.isNaN(a)) {
            assertTrue("" + b + " is not NaN.", Double.isNaN(a));
        } else if (Double.isInfinite(a)) {
            assertTrue(a.equals(b));
        } else {
            double max = FastMath.max(a, b);
            double percentage = FastMath.abs(a - b) / max;
            double deviation = delta;
            assertTrue(String.format("Deviated = %f and is beyond %f as a=%f,  b=%f",
                                     percentage, deviation, a, b), percentage < deviation);
        }
    }

    /**
     * Checks to make sure that the actual quantile position (normalized rank) of value
     * is within tolerance of quantile
     *
     * @param data data array
     * @param value value to test
     * @param quantile purported quantile
     * @param tolerance max difference between actual quantile of value and quantile
     */
    private void checkQuantileError(double[] data, double value, double quantile, double tolerance,
                                    double referenceValue) {
         final double n = (double) data.length;
         int nLess = 0;
         for (double val : data) {
             if (val < value) {
                 nLess++;
             }
         }
         if (Double.isNaN(referenceValue)) {
             assertTrue(Double.isNaN(value));
         } else if (Double.isInfinite(value)) {
             assertTrue(value == referenceValue);
         } else {
             assertTrue("Quantile error exceeded: value returned = " + value +
                        " Reference value = " + referenceValue +
                        " quantile = " + quantile + " n = " + n +
                        " error = " + (quantile - (double) nLess / n),
                        FastMath.abs(quantile - (double) nLess / n) < tolerance);
         }
    }


    private void doCalculatePercentile(Double percentile, Number[] test) {
        doCalculatePercentile(percentile, test, Double.MAX_VALUE);
    }


    private void doCalculatePercentile(Double percentile, Number[] test, double delta) {
        RandomPercentile random = new RandomPercentile(new Well19937c(200));
        for (Number value : test) {
            random.increment(value.doubleValue());
        }

        Percentile p2 = new Percentile(percentile * 100);

        double[] dall = new double[test.length];
        for (int i = 0; i < test.length; i++) {
            dall[i] = test[i].doubleValue();
        }

        Double referenceValue = p2.evaluate(dall);
        assertValues(random.getResult(percentile), referenceValue, delta);
    }

    private void doCalculatePercentile(double percentile, double[] test, double delta) {
        RandomPercentile randomEstimated = new RandomPercentile(new Well19937c(200));
        for (double value : test) {
            randomEstimated.increment(value);
        }

        Percentile p2 = new Percentile(percentile < 1 ? percentile * 100 : percentile);
        /*
         * double[] dall = new double[test.length]; for (int i = 0; i <
         * test.length; i++) dall[i] = test[i];
         */
        Double referenceValue = p2.evaluate(test);
        if (test.length < LARGE) {
            assertValues(randomEstimated.getResult(percentile), referenceValue, delta);
        } else {
            checkQuantileError(test,randomEstimated.getResult(percentile), percentile / 100, delta, referenceValue);
        }
    }

    @Test
    public void testCannedDataSet() {
        Integer[] seedInput =
                new Integer[] { 283, 285, 298, 304, 310, 31, 319, 32, 33, 339,
                        342, 348, 350, 354, 354, 357, 36, 36, 369, 37, 37, 375,
                        378, 383, 390, 396, 405, 408, 41, 414, 419, 416, 42,
                        420, 430, 430, 432, 444, 447, 447, 449, 45, 451, 456,
                        468, 470, 471, 474, 600, 695, 70, 83, 97, 109, 113, 128 };
        Integer[] input = new Integer[seedInput.length * 100];
        for (int i = 0; i < input.length; i++) {
            input[i] = seedInput[i % seedInput.length] + i;
        }
        doCalculatePercentile(0.50d, input);
        doCalculatePercentile(0.95d, input);
    }

    @Test
    public void test99Percentile() {
        Double[] test = randomTestData(100, 10000);
        doCalculatePercentile(0.99d, test);
    }

    @Test
    public void test90Percentile() {
        Double[] test = randomTestData(100, 10000);
        doCalculatePercentile(0.90d, test);
    }

    @Test
    public void test20Percentile() {
        Double[] test = randomTestData(100, 100000);
        doCalculatePercentile(0.20d, test);
    }

    @Test
    public void test5Percentile() {
        Double[] test = randomTestData(50, 990000);
        doCalculatePercentile(0.50d, test);
    }

    @Test
    public void test99PercentileHighValues() {
        Double[] test = randomTestData(100000, 10000);
        doCalculatePercentile(0.99d, test);
    }

    @Test
    public void test90PercentileHighValues() {
        Double[] test = randomTestData(100000, 100000);
        doCalculatePercentile(0.90d, test);
    }

    @Test
    public void test20PercentileHighValues() {
        Double[] test = randomTestData(100000, 100000);
        doCalculatePercentile(0.20d, test);
    }

    @Test
    public void test5PercentileHighValues() {
        Double[] test = randomTestData(100000, 100000);
        doCalculatePercentile(0.05d, test);
    }

    @Test
    public void test0PercentileValuesWithFewerThan5Values() {
        double[] test = { 1d, 2d, 3d, 4d };
        RandomPercentile p = new RandomPercentile();
        assertEquals(1d, p.evaluate(0d,test), 0);
    }


    final int TINY = 10, SMALL = 50, NOMINAL = 100, MEDIUM = 500,
              STANDARD = 1000, BIG = 10000, VERY_BIG = 50000, LARGE = 1000000,
              VERY_LARGE = 10000000;

    private void doDistributionTest(RealDistribution distribution) {
        double data[];

        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(100);
        data = randomDataGenerator.nextDeviates(distribution, LARGE);
        doCalculatePercentile(50, data, 0.0005);
        doCalculatePercentile(95, data, 0.0005);

        data = randomDataGenerator.nextDeviates(distribution, VERY_BIG);
        doCalculatePercentile(50, data, 0.0001);
        doCalculatePercentile(95, data, 0.0001);

        data = randomDataGenerator.nextDeviates(distribution, BIG);
        doCalculatePercentile(50, data, 0.0001);
        doCalculatePercentile(95, data, 0.0001);

        data = randomDataGenerator.nextDeviates(distribution, STANDARD);
        doCalculatePercentile(50, data, 0.0001);
        doCalculatePercentile(95, data, 0.0001);

        data = randomDataGenerator.nextDeviates(distribution, MEDIUM);
        doCalculatePercentile(50, data, 0.0001);
        doCalculatePercentile(95, data, 0.0001);

        data = randomDataGenerator.nextDeviates(distribution, NOMINAL);
        doCalculatePercentile(50, data, 0.0001);
        doCalculatePercentile(95, data, 0.0001);

        data = randomDataGenerator.nextDeviates(distribution, SMALL);
        doCalculatePercentile(50, data, 0.0001);
        doCalculatePercentile(95, data, 0.0001);

        data = randomDataGenerator.nextDeviates(distribution, TINY);
        doCalculatePercentile(50, data, 0.0001);
        doCalculatePercentile(95, data, 0.0001);
    }

    /**
     * Test Various Dist
     */
    @Test
    public void testDistribution() {
        doDistributionTest(new NormalDistribution(4000, 50));
        doDistributionTest(new LogNormalDistribution(4000, 50));
        doDistributionTest(new ExponentialDistribution(4000));
        doDistributionTest(new GammaDistribution(5d,1d));
    }

    /**
     * Large sample streaming tests.
     */
    @Test
    public void testDistributionStreaming() {
        checkQuartiles(new NormalDistribution(), 5000000, 5E-3);
        checkQuartiles(new ExponentialDistribution(1), 100000, 1E-3);
        checkQuartiles(new GammaDistribution(4d,2d), 100000, 1E-3);
    }

    /**
     * Verify the quantile-accuracy contract using a large sample from dist.
     *
     * @param dist distribution to generate sample data from
     * @param sampleSize size of the generated sample
     * @param tolerance normalized rank error tolerance (epsilon in contract)
     */
    private void checkQuartiles(RealDistribution dist, int sampleSize, double tolerance) {
        final long seed = 1000;
        RandomDataGenerator randomDataGenerator = RandomDataGenerator.of(new MersenneTwister(seed));
        final RandomPercentile randomPercentile = new RandomPercentile(RandomPercentile.DEFAULT_EPSILON,
                                                                       randomDataGenerator);
        for (int i = 0; i < sampleSize; i++) {
            randomPercentile.increment(randomDataGenerator.nextDeviate(dist));
        }
        final double q1 = randomPercentile.getResult(25);
        final double q2 = randomPercentile.getResult();
        final double q3 = randomPercentile.getResult(75);

        // Respin the sample and (brutally) count to compute actual ranks
        randomDataGenerator.setSeed(seed);
        double v;
        double ct1 = 0;
        double ct2 = 0;
        double ct3 = 0;
        for (int i = 0; i < sampleSize; i++) {
            v = randomDataGenerator.nextDeviate(dist);
            if (v < q1) {
                ct1++;
                ct2++;
                ct3++;
            } else if (v < q2) {
                ct2++;
                ct3++;
            } else if (v < q3) {
                ct3++;
            }
        }
        assertEquals(0.25, ct1/sampleSize, tolerance);
        assertEquals(0.5, ct2/sampleSize, tolerance);
        assertEquals(0.75, ct3/sampleSize, tolerance);
    }
}
