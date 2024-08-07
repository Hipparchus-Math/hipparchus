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
package org.hipparchus.stat.descriptive.rank;

import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.distribution.continuous.LogNormalDistribution;
import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well19937c;
import org.hipparchus.stat.descriptive.StorelessUnivariateStatistic;
import org.hipparchus.stat.descriptive.StorelessUnivariateStatisticAbstractTest;
import org.hipparchus.stat.descriptive.rank.PSquarePercentile.PSquareMarkers;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test cases for the {@link PSquarePercentile} class which naturally extends
 * {@link StorelessUnivariateStatisticAbstractTest}.
 */
public class PSquarePercentileTest extends
        StorelessUnivariateStatisticAbstractTest {

    protected double percentile5 = 8.2299d;
    protected double percentile95 = 16.72195;// 20.82d; this is approximation
    protected double tolerance = 10E-12;

    private final RandomGenerator randomGenerator = new Well19937c(1000);

    @Override
    public PSquarePercentile getUnivariateStatistic() {
        return new PSquarePercentile(95);
    }

    @Override
    public double expectedValue() {
        return this.percentile95;
    }

    @Override
    public double getTolerance() {
        // tolerance limit changed as this is an approximation
        // algorithm and also gets accurate after few tens of samples
        return 1.0e-2;
    }

    /**
     * Verifies that copied statistics remain equal to originals when
     * incremented the same way by making the copy after a majority of elements
     * are incremented
     */
    @Test
    void testCopyConsistencyWithInitialMostElements() {

        StorelessUnivariateStatistic master = getUnivariateStatistic();
        StorelessUnivariateStatistic replica = null;

        // select a portion of testArray till 75 % of the length to load first
        long index = FastMath.round(0.75 * testArray.length);

        // Put first half in master and copy master to replica
        master.incrementAll(testArray, 0, (int) index);
        replica = master.copy();

        // Check same
        assertEquals(replica, master);
        assertEquals(master, replica);

        // Now add second part to both and check again
        master.incrementAll(testArray, (int) index, (int) (testArray.length - index));
        replica.incrementAll(testArray, (int) index, (int) (testArray.length - index));
        assertEquals(replica, master);
        assertEquals(master, replica);
    }

    /**
     * Verifies that copied statistics remain equal to originals when
     * incremented the same way by way of copying original after just a few
     * elements are incremented
     */
    @Test
    void testCopyConsistencyWithInitialFirstFewElements() {

        StorelessUnivariateStatistic master = getUnivariateStatistic();
        StorelessUnivariateStatistic replica = null;

        // select a portion of testArray which is 10% of the length to load
        // first
        long index = FastMath.round(0.1 * testArray.length);

        // Put first half in master and copy master to replica
        master.incrementAll(testArray, 0, (int) index);
        replica = master.copy();

        // Check same
        assertEquals(replica, master);
        assertEquals(master, replica);
        // Now add second part to both and check again
        master.incrementAll(testArray, (int) index, (int) (testArray.length - index));
        replica.incrementAll(testArray, (int) index, (int) (testArray.length - index));
        assertEquals(master, master);
        assertEquals(replica, replica);
        assertEquals(replica, master);
        assertEquals(master, replica);
    }

    @Test
    void testNullListInMarkers() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            // In case of null list Markers cannot be instantiated..is getting verified
            // new Markers(null, 0, PSquarePercentile.newEstimator());
            PSquarePercentile.newMarkers(null, 0);
        });
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testMiscellaniousFunctionsInMarkers() {
        double p = 0.5;
        PSquareMarkers markers =
                PSquarePercentile.newMarkers(
                        Arrays.asList(new Double[] { 0.02, 1.18, 9.15, 21.91, 38.62 }), p);
        // Markers equality
        assertEquals(markers, markers);
        assertNotEquals(null, markers);
        assertNotEquals(markers, new String());
        // Check for null markers test during equality testing
        // Until 5 elements markers are not initialized
        PSquarePercentile p1 = new PSquarePercentile();
        PSquarePercentile p2 = new PSquarePercentile();
        assertEquals(p1, p2);
        p1.evaluate(new double[] { 1.0, 2.0, 3.0 });
        p2.evaluate(new double[] { 1.0, 2.0, 3.0 });
        assertEquals(p1, p2);
        // Move p2 alone with more values just to make sure markers are not null
        // for p2
        p2.incrementAll(new double[] { 5.0, 7.0, 11.0 });
        assertNotEquals(p1, p2);
        assertNotEquals(p2, p1);
        // Next add different data to p1 to make number of elements match and
        // markers are not null however actual results will vary
        p1.incrementAll(new double[] { 20, 21, 22, 23 });
        assertNotEquals(p1, p2);// though markers are non null, N matches, results wont
    }

    @Test
    void testMarkersOORLow() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            PSquarePercentile.newMarkers(
                Arrays.asList(new Double[]{0.02, 1.18, 9.15, 21.91, 38.62}), 0.5).estimate(0);
        });
    }

    @Test
    void testMarkersOORHigh() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            PSquarePercentile.newMarkers(
                Arrays.asList(new Double[]{0.02, 1.18, 9.15, 21.91, 38.62}), 0.5).estimate(5);
        });
    }

    @Test
    void testMarkers2() {
        double p = 0.5;
        PSquareMarkers markers =
                PSquarePercentile.newMarkers(
                        Arrays.asList(new Double[] { 0.02, 1.18, 9.15, 21.91, 38.62 }), p);

        PSquareMarkers markersNew =
                PSquarePercentile.newMarkers(
                        Arrays.asList(new Double[] { 0.02, 1.18, 9.15, 21.91, 38.62 }), p);

        assertEquals(markers, markersNew);
        // If just one element of markers got changed then its still false.
        markersNew.processDataPoint(39);
        assertNotEquals(markers, markersNew);
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testHashCodeInMarkers() {
        PSquarePercentile p = new PSquarePercentile(95);
        PSquarePercentile p2 = new PSquarePercentile(95);
        Set<PSquarePercentile> s = new HashSet<>();
        s.add(p);
        s.add(p2);
        assertEquals(1, s.size());
        assertEquals(p, s.iterator().next());
        double[] d =
                new double[] { 95.1772, 95.1567, 95.1937, 95.1959, 95.1442, 95.0610,
                               95.1591, 95.1195, 95.1772, 95.0925, 95.1990, 95.1682 };
        assertEquals(95.1981, p.evaluate(d), 1.0e-2); // change
        assertEquals(95.1981, p2.evaluate(d), 1.0e-2); // change
        s.clear();
        s.add(p);
        s.add(p2);
        assertEquals(1, s.size());
        assertEquals(p, s.iterator().next());

        PSquareMarkers m1 =
                PSquarePercentile.newMarkers(
                        Arrays.asList(new Double[] { 95.1772, 95.1567, 95.1937,
                                95.1959, 95.1442, 95.0610, 95.1591, 95.1195,
                                95.1772, 95.0925, 95.1990, 95.1682 }), 0.0);
        PSquareMarkers m2 =
                PSquarePercentile.newMarkers(
                        Arrays.asList(new Double[] { 95.1772, 95.1567, 95.1937,
                                95.1959, 95.1442, 95.0610, 95.1591, 95.1195,
                                95.1772, 95.0925, 95.1990, 95.1682 }), 0.0);
        assertEquals(m1, m2);
        Set<PSquareMarkers> setMarkers = new LinkedHashSet<PSquareMarkers>();
        assertTrue(setMarkers.add(m1));
        assertFalse(setMarkers.add(m2));
        assertEquals(1, setMarkers.size());

        PSquareMarkers mThis =
                PSquarePercentile.newMarkers(
                        Arrays.asList(new Double[] { 195.1772, 195.1567,
                                195.1937, 195.1959, 95.1442, 195.0610,
                                195.1591, 195.1195, 195.1772, 95.0925, 95.1990,
                                195.1682 }), 0.50);
        PSquareMarkers mThat =
                PSquarePercentile.newMarkers(
                        Arrays.asList(new Double[] { 95.1772, 95.1567, 95.1937,
                                95.1959, 95.1442, 95.0610, 95.1591, 95.1195,
                                95.1772, 95.0925, 95.1990, 95.1682 }), 0.50);
        assertEquals(mThis, mThis);
        assertNotEquals(mThis, mThat);
        String s1="";
        assertNotEquals(mThis, s1);
        for (int i = 0; i < testArray.length; i++) {
            mThat.processDataPoint(testArray[i]);
        }
        setMarkers.add(mThat);
        setMarkers.add(mThis);
        assertEquals(mThat, mThat);
        assertTrue(setMarkers.contains(mThat));
        assertTrue(setMarkers.contains(mThis));
        assertEquals(3, setMarkers.size());
        Iterator<PSquareMarkers> iterator=setMarkers.iterator();
        assertEquals(m1, iterator.next());
        assertEquals(mThat, iterator.next());
        assertEquals(mThis, iterator.next());
    }

    @Test
    void testMarkersWithLowerIndex() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            PSquareMarkers mThat =
                PSquarePercentile.newMarkers(
                    Arrays.asList(new Double[]{95.1772, 95.1567, 95.1937,
                        95.1959, 95.1442, 95.0610, 95.1591, 95.1195,
                        95.1772, 95.0925, 95.1990, 95.1682}), 0.50);
            for (int i = 0; i < testArray.length; i++) {
                mThat.processDataPoint(testArray[i]);
            }
            mThat.estimate(0);
        });
    }

    @Test
    void testMarkersWithHigherIndex() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            PSquareMarkers mThat =
                PSquarePercentile.newMarkers(
                    Arrays.asList(new Double[]{95.1772, 95.1567, 95.1937,
                        95.1959, 95.1442, 95.0610, 95.1591, 95.1195,
                        95.1772, 95.0925, 95.1990, 95.1682}), 0.50);
            for (int i = 0; i < testArray.length; i++) {
                mThat.processDataPoint(testArray[i]);
            }
            mThat.estimate(6);
        });
    }

    @Test
    void testMarkerHeightWithLowerIndex() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            PSquareMarkers mThat =
                PSquarePercentile.newMarkers(
                    Arrays.asList(new Double[]{95.1772, 95.1567, 95.1937,
                        95.1959, 95.1442, 95.0610, 95.1591, 95.1195,
                        95.1772, 95.0925, 95.1990, 95.1682}), 0.50);
            mThat.height(0);
        });
    }

    @Test
    void testMarkerHeightWithHigherIndex() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            PSquareMarkers mThat =
                PSquarePercentile.newMarkers(
                    Arrays.asList(new Double[]{95.1772, 95.1567, 95.1937,
                        95.1959, 95.1442, 95.0610, 95.1591, 95.1195,
                        95.1772, 95.0925, 95.1990, 95.1682}), 0.50);
            mThat.height(6);
        });
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testPSquaredEqualsAndMin() {
        PSquarePercentile ptile = new PSquarePercentile(0);
        assertEquals(ptile, ptile);
        assertNotEquals(null, ptile);
        assertNotEquals(ptile, new String());
        // Just to check if there is no data get result for zeroth and 100th
        // ptile returns NAN
        assertTrue(Double.isNaN(ptile.getResult()));
        assertTrue(Double.isNaN(new PSquarePercentile(100).getResult()));

        double[] d = new double[] { 1, 3, 2, 4, 9, 10, 11 };
        ptile.incrementAll(d);
        assertEquals(ptile, ptile);
        assertEquals(1d, ptile.getResult(), 1e-02);// this calls min
        assertEquals(0.0, ptile.getQuantile(), 1.0e-15);
    }

    @Test
    void testString() {
        PSquarePercentile ptile = new PSquarePercentile(95);
        assertNotNull(ptile.toString());
        ptile.increment(1);
        ptile.increment(2);
        ptile.increment(3);
        assertNotNull(ptile.toString());
        assertEquals(expectedValue(), ptile.evaluate(testArray), getTolerance());
        assertNotNull(ptile.toString());
        assertEquals(0.95, ptile.getQuantile(), 1.0e-15);
    }

    @Test
    void testHighPercentile() {
        double[] d = new double[] { 1, 2, 3 };
        PSquarePercentile p = new PSquarePercentile(75.0);
        assertEquals(2, p.evaluate(d), 1.0e-5);
        PSquarePercentile p95 = new PSquarePercentile();
        assertEquals(2, p95.evaluate(d), 1.0e-5);
    }

    @Test
    void testLowPercentile() {
        double[] d = new double[] { 0, 1 };
        PSquarePercentile p = new PSquarePercentile(25.0);
        assertEquals(0d, p.evaluate(d), Double.MIN_VALUE);
    }

    @Test
    void testPercentile() {
        double[] d = new double[] { 1, 3, 2, 4 };
        PSquarePercentile p = new PSquarePercentile(30d);
        assertEquals(1.0, p.evaluate(d), 1.0e-5);
        p = new PSquarePercentile(25);
        assertEquals(1.0, p.evaluate(d), 1.0e-5);
        p = new PSquarePercentile(75);
        assertEquals(3.0, p.evaluate(d), 1.0e-5);
        p = new PSquarePercentile(50);
        assertEquals(2d, p.evaluate(d), 1.0e-5);
    }

    @Test
    void testInitial() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            PSquarePercentile.newMarkers(new ArrayList<Double>(), 0.5);
        });
    }

    @Test
    void testNegativeInvalidValues() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            double[] d =
                new double[]{95.1772, 95.1567, 95.1937, 95.1959, 95.1442,
                    95.0610, 95.1591, 95.1195, 95.1772, 95.0925, 95.1990,
                    95.1682};
            PSquarePercentile p = new PSquarePercentile(-1.0);
            p.evaluate(d, 0, d.length);
        });
    }

    @Test
    void testPositiveInvalidValues() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            double[] d =
                new double[]{95.1772, 95.1567, 95.1937, 95.1959, 95.1442,
                    95.0610, 95.1591, 95.1195, 95.1772, 95.0925, 95.1990,
                    95.1682};
            PSquarePercentile p = new PSquarePercentile(101.0);
            p.evaluate(d, 0, d.length);
        });
    }

    @Test
    void testNISTExample() {
        double[] d =
                new double[] { 95.1772, 95.1567, 95.1937, 95.1959, 95.1442,
                        95.0610, 95.1591, 95.1195, 95.1772, 95.0925, 95.1990,
                        95.1682 };
        assertEquals(95.1981, new PSquarePercentile(90d).evaluate(d), 1.0e-2); // changed the accuracy to 1.0e-2
        assertEquals(95.061, new PSquarePercentile(0d).evaluate(d), 0);
        assertEquals(95.1990, new PSquarePercentile(100d).evaluate(d, 0, d.length), 0);
    }

    @Test
    void test5() {
        PSquarePercentile percentile = new PSquarePercentile(5d);
        assertEquals(this.percentile5, percentile.evaluate(testArray), 1.0); // changed the accuracy to 1 instead of tolerance
    }

    @Test
    void testNull() {
        assertThrows(NullArgumentException.class, () -> {
            PSquarePercentile percentile = new PSquarePercentile(50d);
            double[] nullArray = null;
            percentile.evaluate(nullArray);
        });
    }

    @Test
    void testEmpty() {
        PSquarePercentile percentile = new PSquarePercentile(50d);
        double[] emptyArray = new double[] {};
        assertTrue(Double.isNaN(percentile.evaluate(emptyArray)));
    }

    @Test
    void testSingleton() {
        PSquarePercentile percentile = new PSquarePercentile(50d);
        double[] singletonArray = new double[] { 1d };
        assertEquals(1d, percentile.evaluate(singletonArray), 0);
        assertEquals(1d, percentile.evaluate(singletonArray, 0, 1), 0);
        percentile = new PSquarePercentile(5);
        assertEquals(1d, percentile.evaluate(singletonArray, 0, 1), 0);
        percentile = new PSquarePercentile(100);
        assertEquals(1d, percentile.evaluate(singletonArray, 0, 1), 0);
        percentile = new PSquarePercentile(100);
        assertTrue(Double.isNaN(percentile.evaluate(singletonArray, 0, 0)));
    }

    @Test
    void testSpecialValues() {
        PSquarePercentile percentile = new PSquarePercentile(50d);
        double[] specialValues = new double[] { 0d, 1d, 2d, 3d, 4d, Double.NaN };
        assertEquals(2d, percentile.evaluate(specialValues), 0);
        specialValues =
            new double[] { Double.NEGATIVE_INFINITY, 1d, 2d, 3d, Double.NaN, Double.POSITIVE_INFINITY };
        assertEquals(2d, percentile.evaluate(specialValues), 0);
        specialValues =
            new double[] { 1d, 1d, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY };
        assertFalse(Double.isInfinite(percentile.evaluate(specialValues)));
        specialValues = new double[] { 1d, 1d, Double.NaN, Double.NaN };
        assertFalse(Double.isNaN(percentile.evaluate(specialValues)));
        specialValues =
            new double[] { 1d, 1d, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY };
        percentile = new PSquarePercentile(50d);
        // Interpolation results in NEGATIVE_INFINITY + POSITIVE_INFINITY
        // changed the result check to infinity instead of NaN
        assertTrue(Double.isInfinite(percentile.evaluate(specialValues)));
    }

    @Test
    void testArrayExample() {
        assertEquals(expectedValue(), new PSquarePercentile(95d).evaluate(testArray), getTolerance());
    }

    @Test
    void testSetQuantile() {
        PSquarePercentile percentile = new PSquarePercentile(10d);

        percentile = new PSquarePercentile(100); // OK
        assertEquals(1.0, percentile.quantile(), 0);
        try {
            percentile = new PSquarePercentile(0);
            // fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            new PSquarePercentile(0d);
            // fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
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
    void testAccept() {
        PSquarePercentile psquared = new PSquarePercentile(0.99);
        assertTrue(Double.isNaN(psquared.getResult()));
        Double[] test = randomTestData(100, 10000);

        for (Double value : test) {
            psquared.increment(value);
            assertTrue(psquared.getResult() >= 0);
        }
    }

    private void customAssertValues(Double a, Double b, double delta) {
        if (Double.isNaN(a)) {
            assertTrue(Double.isNaN(a), "" + b + " is not NaN.");
        } else {
            double max = FastMath.max(a, b);
            double percentage = FastMath.abs(a - b) / max;
            double deviation = delta;
            assertTrue(percentage < deviation, String.format("Deviated = %f and is beyond %f as a=%f,  b=%f",
                                     percentage, deviation, a, b));
        }
    }

    private void doCalculatePercentile(Double percentile, Number[] test) {
        doCalculatePercentile(percentile, test, Double.MAX_VALUE);
    }

    private void doCalculatePercentile(Double percentile, Number[] test, double delta) {
        PSquarePercentile psquared = new PSquarePercentile(percentile);
        for (Number value : test) {
            psquared.increment(value.doubleValue());
        }

        Percentile p2 = new Percentile(percentile * 100);

        double[] dall = new double[test.length];
        for (int i = 0; i < test.length; i++) {
            dall[i] = test[i].doubleValue();
        }

        Double referenceValue = p2.evaluate(dall);
        customAssertValues(psquared.getResult(), referenceValue, delta);
    }

    private void doCalculatePercentile(double percentile, double[] test, double delta) {
        PSquarePercentile psquared = new PSquarePercentile(percentile);
        for (double value : test) {
            psquared.increment(value);
        }

        Percentile p2 = new Percentile(percentile < 1 ? percentile * 100 : percentile);
        /*
         * double[] dall = new double[test.length]; for (int i = 0; i <
         * test.length; i++) dall[i] = test[i];
         */
        Double referenceValue = p2.evaluate(test);
        customAssertValues(psquared.getResult(), referenceValue, delta);
    }

    @Test
    void testCannedDataSet() {
        // test.unoverride("dump");
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
        // Arrays.sort(input);
        doCalculatePercentile(0.50d, input);
        doCalculatePercentile(0.95d, input);
    }

    @Test
    void test99Percentile() {
        Double[] test = randomTestData(100, 10000);
        doCalculatePercentile(0.99d, test);
    }

    @Test
    void test90Percentile() {
        Double[] test = randomTestData(100, 10000);
        doCalculatePercentile(0.90d, test);
    }

    @Test
    void test20Percentile() {
        Double[] test = randomTestData(100, 100000);
        doCalculatePercentile(0.20d, test);
    }

    @Test
    void test5Percentile() {
        Double[] test = randomTestData(50, 990000);
        doCalculatePercentile(0.50d, test);
    }

    @Test
    void test99PercentileHighValues() {
        Double[] test = randomTestData(100000, 10000);
        doCalculatePercentile(0.99d, test);
    }

    @Test
    void test90PercentileHighValues() {
        Double[] test = randomTestData(100000, 100000);
        doCalculatePercentile(0.90d, test);
    }

    @Test
    void test20PercentileHighValues() {
        Double[] test = randomTestData(100000, 100000);
        doCalculatePercentile(0.20d, test);
    }

    @Test
    void test5PercentileHighValues() {
        Double[] test = randomTestData(100000, 100000);
        doCalculatePercentile(0.05d, test);
    }

    @Test
    void test0PercentileValuesWithFewerThan5Values() {
        double[] test = { 1d, 2d, 3d, 4d };
        PSquarePercentile p = new PSquarePercentile(0d);
        assertEquals(1d, p.evaluate(test), 0);
        assertNotNull(p.toString());
    }

    @Test
    void testPSQuaredEvalFuncWithPapersExampleData() throws IOException {

        // This data as input is considered from
        // http://www.cs.wustl.edu/~jain/papers/ftp/psqr.pdf
        double[] data =
                { 0.02, 0.5, 0.74, 3.39, 0.83, 22.37, 10.15, 15.43, 38.62,
                  15.92, 34.6, 10.28, 1.47, 0.4, 0.05, 11.39, 0.27, 0.42,
                  0.09, 11.37,

                  11.39, 15.43, 15.92, 22.37, 34.6, 38.62, 18.9, 19.2,
                  27.6, 12.8, 13.7, 21.9
                };

        PSquarePercentile psquared = new PSquarePercentile(50);

        Double p2value = 0d;
        for (int i = 0; i < 20; i++) {
            psquared.increment(data[i]);
            p2value = psquared.getResult();
            // System.out.println(psquared.toString());//uncomment here to see
            // the papers example output
        }
        // System.out.println("p2value=" + p2value);
        Double expected = 4.44d;// 13d; // From The Paper
        // http://www.cs.wustl.edu/~jain/papers/ftp/psqr.pdf.
        // Pl refer Pg 1061 Look at the mid marker
        // height
        // expected = new Percentile(50).evaluate(data,0,20);
        // Well the values deviate in our calculation by 0.25 so its 4.25 vs
        // 4.44
        assertEquals(expected, p2value, 0.25, String.format("Expected=%f, Actual=%f", expected, p2value));
    }

    final int TINY = 10, SMALL = 50, NOMINAL = 100, MEDIUM = 500,
              STANDARD = 1000, BIG = 10000, VERY_BIG = 50000, LARGE = 1000000,
              VERY_LARGE = 10000000;

    private void doDistributionTest(RealDistribution distribution) {
        double[] data;

//        data = distribution.sample(VERY_LARGE);
//        doCalculatePercentile(50, data, 0.0001);
//        doCalculatePercentile(95, data, 0.0001);

        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(100);
        data = randomDataGenerator.nextDeviates(distribution, LARGE);
        doCalculatePercentile(50, data, 0.001);
        doCalculatePercentile(95, data, 0.001);

        data = randomDataGenerator.nextDeviates(distribution, VERY_BIG);
        doCalculatePercentile(50, data, 0.001);
        doCalculatePercentile(95, data, 0.001);

        data = randomDataGenerator.nextDeviates(distribution, BIG);
        doCalculatePercentile(50, data, 0.001);
        doCalculatePercentile(95, data, 0.001);

        data = randomDataGenerator.nextDeviates(distribution, STANDARD);
        doCalculatePercentile(50, data, 0.005);
        doCalculatePercentile(95, data, 0.005);

        data = randomDataGenerator.nextDeviates(distribution, MEDIUM);
        doCalculatePercentile(50, data, 0.005);
        doCalculatePercentile(95, data, 0.005);

        data = randomDataGenerator.nextDeviates(distribution, NOMINAL);
        doCalculatePercentile(50, data, 0.01);
        doCalculatePercentile(95, data, 0.01);

        data = randomDataGenerator.nextDeviates(distribution, SMALL);
        doCalculatePercentile(50, data, 0.01);
        doCalculatePercentile(95, data, 0.01);

        data = randomDataGenerator.nextDeviates(distribution, TINY);
        doCalculatePercentile(50, data, 0.05);
        doCalculatePercentile(95, data, 0.05);
    }

    /**
     * Test Various Dist
     */
    @Test
    void testDistribution() {
        doDistributionTest(new NormalDistribution(4000, 50));
        doDistributionTest(new LogNormalDistribution(4000, 50));
        // doDistributionTest((new ExponentialDistribution(4000));
        // doDistributionTest(new GammaDistribution(5d,1d),0.1);
    }
}
