/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.stat.descriptive.rank;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.SplittableRandom;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.distribution.continuous.ExponentialDistribution;
import org.hipparchus.distribution.continuous.GammaDistribution;
import org.hipparchus.distribution.continuous.LogNormalDistribution;
import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.random.MersenneTwister;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well19937c;
import org.hipparchus.stat.StatUtils;
import org.hipparchus.stat.descriptive.StorelessUnivariateStatistic;
import org.hipparchus.stat.descriptive.StorelessUnivariateStatisticAbstractTest;
import org.hipparchus.stat.inference.AlternativeHypothesis;
import org.hipparchus.stat.inference.BinomialTest;
import org.hipparchus.util.FastMath;
import org.junit.Test;

/**
 * Test cases for the {@link RandomPercentileTest} class.
 * Some tests are adapted from PSquarePercentileTest.
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
    public void testAggregateSmallSamplesA() {
        doTestAggregateSmallSamples(0.5, 11.0);
    }

    @Test
    public void testAggregateSmallSamplesB() {
        doTestAggregateSmallSamples(0.1, 7.0);
    }

    @Test
    public void testAggregateSmallSamplesC() {
        doTestAggregateSmallSamples(0.01, 7.0);
    }

    private void doTestAggregateSmallSamples(double epsilon, double expected) {
        SplittableRandom seeds = new SplittableRandom(0x218560e08c8df220l);

        RandomPercentile rp1   = new RandomPercentile(epsilon, new Well19937c(seeds.nextLong()));
        double[]         data1 = new double[] { 3, 5, -2, 7, 14, 6 };
        for (double d : data1 ) {
            rp1.accept(d);
        }
        RandomPercentile rp2   = new RandomPercentile(epsilon, new Well19937c(seeds.nextLong()));
        double[]         data2 = new double[] { 9, 12, 15 };
        for (double d : data2 ) {
            rp2.accept(d);
        }

        rp1.aggregate(rp2);

        assertEquals(expected, rp1.getResult(), 1.0e-10);

    }

    @Test
    public void testBufferConsumeLevel0() throws Exception {
        final int len = testArray.length;
        final double[] sorted = Arrays.copyOf(testArray, len);
        Arrays.sort(sorted);
        BufferMock buffer = new BufferMock(len, 0, randomGenerator);
        for (int i = 0; i < len; i++) {
            buffer.consume(testArray[i]);
        }
        UnitTestUtils.assertEquals(sorted, buffer.getData(), Double.MIN_VALUE);
        assertFalse(buffer.hasCapacity());
        assertEquals(StatUtils.min(testArray),buffer.min(), Double.MIN_VALUE);
        assertEquals(StatUtils.max(testArray),buffer.max(), Double.MIN_VALUE);
    }

    @Test
    public void testBufferSampling() throws Exception {
        // Capacity = 10, level = 2 - should take 1 out of each block of 4 values from the stream
        BufferMock buffer = new BufferMock(10, 2, randomGenerator);
        for (int i = 0; i < 40; i++) {
            buffer.consume(i);
        }
        // Buffer should be full, including one value from each consecutive sequence of 4
        assertFalse(buffer.hasCapacity());
        final double[] data = buffer.getData();
        for (int i = 0; i < 10; i++) {
           assertTrue(data[i] < 4 * (i + 1));
           assertTrue(data[i] >= 4 * i);
           // Rank of 4n should be n for n = 1, ..., 10.
           assertEquals(i + 1, buffer.rankOf((i + 1) * 4));
           // Check boundary ranks
           assertEquals(0, buffer.rankOf(-1));
           assertEquals(10, buffer.rankOf(100));
        }
    }

    @Test
    public void testBufferMergeWith() throws Exception {
        // Create 2 level 0 buffers of size 20
        BufferMock buffer1 = new BufferMock(20, 0, randomGenerator);
        BufferMock buffer2 = new BufferMock(20, 0, randomGenerator);

        // fill buffer1 with 0, ..., 19 and buffer2 with 20, ..., 39
        for (int i = 0; i < 20; i++) {
            buffer1.consume(i);
            buffer2.consume(i + 20);
        }
        assertEquals(0, buffer1.getLevel());
        assertEquals(0, buffer2.getLevel());

        // Merge 1 with 2
        buffer1.mergeWith(buffer2.getInstance()); // Need to pass the real thing for buffer2

        // Both should now have level 1, buffer1 should be full, buffer2 should be free
        assertEquals(1,buffer1.getLevel());
        assertFalse(buffer1.hasCapacity());
        assertEquals(1,buffer2.getLevel());
        assertTrue(buffer2.hasCapacity());

        // Check the contents of buffer1 - should be the merged contents of both
        final double[] data = buffer1.getData();
        int nSmall = 0;
        for (double value : data) {
            assertTrue(value >=0 && value < 40);
            if (value < 20) {
                nSmall++;
            }
        }

        // Verify merge selection was close to fair
        BinomialTest bTest = new BinomialTest();
        assertFalse(bTest.binomialTest(20, nSmall, 0.5, AlternativeHypothesis.TWO_SIDED, 0.01));

        // Check sort on merged buffer
        buffer1.checkSorted();
    }

    @Test
    public void testBufferMergeInto() throws Exception {
        BufferMock buffer1 = new BufferMock(20, 0, randomGenerator);
        BufferMock buffer2 = new BufferMock(20, 2, randomGenerator);

        // fill buffer1 with 0, ..., 19
        for (int i = 0; i < 20; i++) {
            buffer1.consume(i);
        }

        // fill buffer 2 - level 2 means it will take 80 values to fill it
        for (int i = 0; i < 80; i++) {
            buffer2.consume(i + 20);
        }
        assertFalse(buffer1.hasCapacity());
        assertFalse(buffer2.hasCapacity());
        buffer1.mergeInto(buffer2.getInstance());

        // levels should be unchanged
        assertEquals(0, buffer1.getLevel());
        assertEquals(2, buffer2.getLevel());

        // buffer2 should have about 1/4 values under 20
        final double[] data = buffer2.getData();
        int nSmall = 0;
        for (double value : data) {
            assertTrue(value >=0 && value < 100);
            if (value < 20) {
                nSmall++;
            }
        }
        BinomialTest bTest = new BinomialTest();
        assertFalse(bTest.binomialTest(20, nSmall, 0.25, AlternativeHypothesis.TWO_SIDED, 0.01));

        // Check sort on merged buffer
        buffer2.checkSorted();
    }

    static class BufferMock {
        private Object instance;
        private Method consumeMethod;
        private Method getDataMethod;
        private Method hasCapacityMethod;
        private Method minMethod;
        private Method maxMethod;
        private Method rankOfMethod;
        private Method getLevelMethod;
        private Method mergeWithMethod;
        private Method mergeIntoMethod;
        @SuppressWarnings({ "unchecked", "rawtypes" })
        public BufferMock(int size, int level, RandomGenerator randomGenerator) throws Exception {
            final Class[] classes = RandomPercentile.class.getDeclaredClasses();
            for (Class cls : classes) {
                if (cls.getName().endsWith("$Buffer")) {
                   Constructor constructor = cls.getDeclaredConstructor(Integer.TYPE, Integer.TYPE, RandomGenerator.class);
                   instance = constructor.newInstance(size, level, randomGenerator);
                   consumeMethod = cls.getDeclaredMethod("consume", Double.TYPE);
                   getDataMethod = cls.getDeclaredMethod("getData");
                   hasCapacityMethod = cls.getDeclaredMethod("hasCapacity");
                   minMethod = cls.getDeclaredMethod("min");
                   maxMethod = cls.getDeclaredMethod("max");
                   rankOfMethod = cls.getDeclaredMethod("rankOf", Double.TYPE);
                   getLevelMethod = cls.getDeclaredMethod("getLevel");
                   mergeWithMethod = cls.getDeclaredMethod("mergeWith", cls);
                   mergeIntoMethod = cls.getDeclaredMethod("mergeInto", cls);
                }
            }
        }

        public void consume(double value) throws Exception {
            consumeMethod.invoke(instance, value);
        }

        public boolean hasCapacity() throws Exception {
            return (boolean) hasCapacityMethod.invoke(instance);
        }

        public double[] getData() throws Exception {
            return (double[]) getDataMethod.invoke(instance);
        }

        public double min() throws Exception {
            return (double) minMethod.invoke(instance);
        }

        public double max() throws Exception {
            return (double) maxMethod.invoke(instance);
        }

        public int rankOf(double value) throws Exception {
            return (int) rankOfMethod.invoke(instance, value);
        }

        public int getLevel() throws Exception {
            return (int) getLevelMethod.invoke(instance);
        }

        public void mergeWith(Object other) throws Exception {
            mergeWithMethod.invoke(instance, other);
        }

        public void mergeInto(Object other) throws Exception {
            mergeIntoMethod.invoke(instance, other);
        }

        public Object getInstance() {
            return instance;
        }

        public void checkSorted() throws Exception {
            final double[] data = getData();
            final double[] copy = Arrays.copyOf(data, data.length);
            Arrays.sort(copy);
            UnitTestUtils.assertEquals(data, copy, Double.MIN_VALUE);
        }

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

    // @Test
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
    
    @Test
    public void testMaxValuesRetained() {
        assertEquals(546795, RandomPercentile.maxValuesRetained(RandomPercentile.DEFAULT_EPSILON));
        assertEquals(34727, RandomPercentile.maxValuesRetained(1E-3));
        assertEquals(2064, RandomPercentile.maxValuesRetained(1E-2));
    }
    
    @Test(expected = MathIllegalArgumentException.class)
    public void testMaxValuesRetained0Epsilon() {
        RandomPercentile.maxValuesRetained(0);
    }
    
    @Test(expected = MathIllegalArgumentException.class)
    public void testMaxValuesRetained1Epsilon() {
        RandomPercentile.maxValuesRetained(1);
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
        checkQuartiles(new NormalDistribution(), 1000000, 5E-4);
        checkQuartiles(new ExponentialDistribution(1), 600000, 5E-4);
        checkQuartiles(new GammaDistribution(4d,2d), 600000, 5E-4);
    }

    /**
     * Verify no sequential bias even when buffer size is small.
     */
    @Test
    public void testSequentialData() {
        final long seed = 1000;
        double epsilon = 1e-4;
        for (int j = 0; j < 3; j++) {
            epsilon *= 10;
            final RandomPercentile randomPercentile = new RandomPercentile(epsilon,
                    new MersenneTwister(seed));
            final int n = 5000000;
            for (int i = 1; i <= n; i++) {
                randomPercentile.accept(i);
            }
            for (int i = 1; i < 5; i++) {
                assertEquals(0.2 * i, randomPercentile.getResult(i * 20) / n, 2 * epsilon);
            }
        }
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
                randomGenerator);
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
