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

import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.stat.descriptive.UnivariateStatistic;
import org.hipparchus.stat.descriptive.UnivariateStatisticAbstractTest;
import org.hipparchus.stat.descriptive.rank.Percentile.EstimationType;
import org.hipparchus.stat.ranking.NaNStrategy;
import org.hipparchus.util.KthSelector;
import org.hipparchus.util.PivotingStrategy;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the {@link Percentile} class.
 */
public class PercentileTest extends UnivariateStatisticAbstractTest{

    private double quantile;

    /**
     * {@link Percentile.EstimationType type}
     * of estimation to be used while calling {@link #getUnivariateStatistic()}
     */
    private Percentile.EstimationType type;

    /**
     * {@link NaNStrategy}
     * of estimation to be used while calling {@link #getUnivariateStatistic()}
     */
    private NaNStrategy nanStrategy;

    /**
     * kth selector
     */
    private KthSelector kthSelector;

    /**
     * A default percentile to be used for {@link #getUnivariateStatistic()}
     */
    protected final double DEFAULT_PERCENTILE = 95d;

    /**
     * Before method to ensure defaults retained
     */
    @BeforeEach
    void setup() {
        quantile         = 95.0;
        type             = Percentile.EstimationType.LEGACY;
        nanStrategy      = NaNStrategy.REMOVED;
        kthSelector      = new KthSelector(PivotingStrategy.MEDIAN_OF_3);
    }

    private void reset(final double p, final Percentile.EstimationType type) {
        this.quantile = p;
        this.type     = type;
        nanStrategy   = (type == Percentile.EstimationType.LEGACY) ? NaNStrategy.FIXED : NaNStrategy.REMOVED;
    }

    @Override
    public Percentile getUnivariateStatistic() {
        return new Percentile(quantile)
                .withEstimationType(type)
                .withNaNStrategy(nanStrategy)
                .withKthSelector(kthSelector);
    }

    @Override
    public double expectedValue() {
        return this.percentile95;
    }

    @Test
    void testHighPercentile(){
        final double[] d = new double[]{1, 2, 3};
        final Percentile p = new Percentile(75);
        assertEquals(3.0, p.evaluate(d), 1.0e-5);
    }

    @Test
    void testLowPercentile() {
        final double[] d = new double[] {0, 1};
        final Percentile p = new Percentile(25);
        assertEquals(0d, p.evaluate(d), Double.MIN_VALUE);
    }

    @Test
    void testPercentile() {
        final double[] d = new double[] {1, 3, 2, 4};
        final Percentile p = new Percentile(30);
        assertEquals(1.5, p.evaluate(d), 1.0e-5);
        p.setQuantile(25);
        assertEquals(1.25, p.evaluate(d), 1.0e-5);
        p.setQuantile(75);
        assertEquals(3.75, p.evaluate(d), 1.0e-5);
        p.setQuantile(50);
        assertEquals(2.5, p.evaluate(d), 1.0e-5);

        // invalid percentiles
        try {
            p.evaluate(d, 0, d.length, -1.0);
            fail();
        } catch (final MathIllegalArgumentException ex) {
            // success
        }
        try {
            p.evaluate(d, 0, d.length, 101.0);
            fail();
        } catch (final MathIllegalArgumentException ex) {
            // success
        }
    }

    @Test
    void testNISTExample() {
        final double[] d = new double[] {95.1772, 95.1567, 95.1937, 95.1959,
                95.1442, 95.0610,  95.1591, 95.1195, 95.1772, 95.0925, 95.1990, 95.1682
        };
        final Percentile p = new Percentile(90);
        assertEquals(95.1981, p.evaluate(d), 1.0e-4);
        assertEquals(95.1990, p.evaluate(d,0,d.length, 100d), 0);
    }

    @Test
    void test5() {
        final Percentile percentile = new Percentile(5);
        assertEquals(this.percentile5, percentile.evaluate(testArray), getTolerance());
    }

    @Test
    void testNullEmpty() {
        final Percentile percentile = new Percentile(50);
        final double[] nullArray = null;
        final double[] emptyArray = new double[] {};
        try {
            percentile.evaluate(nullArray);
            fail("Expecting NullArgumentException for null array");
        } catch (final NullArgumentException ex) {
            // expected
        }
        assertTrue(Double.isNaN(percentile.evaluate(emptyArray)));
    }

    @Test
    void testSingleton() {
        final Percentile percentile = new Percentile(50);
        final double[] singletonArray = new double[] {1d};
        assertEquals(1d, percentile.evaluate(singletonArray), 0);
        assertEquals(1d, percentile.evaluate(singletonArray, 0, 1), 0);
        assertEquals(1d, percentile.evaluate(singletonArray, 0, 1, 5), 0);
        assertEquals(1d, percentile.evaluate(singletonArray, 0, 1, 100), 0);
        assertTrue(Double.isNaN(percentile.evaluate(singletonArray, 0, 0)));
    }

    @Test
    void testSpecialValues() {
        final Percentile percentile = new Percentile(50);
        double[] specialValues = new double[] {0d, 1d, 2d, 3d, 4d,  Double.NaN};
        assertEquals(/*2.5d*/2d, percentile.evaluate(specialValues), 0);
        specialValues = new double[] { Double.NEGATIVE_INFINITY, 1d, 2d, 3d, Double.NaN, Double.POSITIVE_INFINITY };
        assertEquals(/*2.5d*/2d, percentile.evaluate(specialValues), 0);
        specialValues = new double[] {1d, 1d, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        assertTrue(Double.isInfinite(percentile.evaluate(specialValues)));
        specialValues = new double[] {1d, 1d, Double.NaN, Double.NaN};
        assertFalse(Double.isNaN(percentile.evaluate(specialValues)));
        assertEquals(1d, percentile.evaluate(specialValues));
        specialValues = new double[] {1d, 1d, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
        // Interpolation results in NEGATIVE_INFINITY + POSITIVE_INFINITY
        assertTrue(Double.isNaN(percentile.evaluate(specialValues)));
    }

    @Test
    void testSetQuantile() {
        final Percentile percentile = new Percentile(10);
        percentile.setQuantile(100); // OK
        assertEquals(100, percentile.getQuantile(), 0);
        try {
            percentile.setQuantile(0);
            fail("Expecting MathIllegalArgumentException");
        } catch (final MathIllegalArgumentException ex) {
            // expected
        }
        try {
            new Percentile(0);
            fail("Expecting MathIllegalArgumentException");
        } catch (final MathIllegalArgumentException ex) {
            // expected
        }
    }

    //Below tests are basically to run for all estimation types.
    /**
     * While {@link #testHighPercentile()} checks only for the existing
     * implementation; this method verifies for all the types including Percentile.Type.CM Percentile.Type.
     */
    @Test
    void testAllTechniquesHighPercentile() {
        final double[] d = new double[] { 1, 2, 3 };
        testAssertMappedValues(d, new Object[][] { { Percentile.EstimationType.LEGACY, 3d }, { Percentile.EstimationType.R_1, 3d },
                { Percentile.EstimationType.R_2, 3d }, { Percentile.EstimationType.R_3, 2d }, { Percentile.EstimationType.R_4, 2.25 }, { Percentile.EstimationType.R_5, 2.75 },
                { Percentile.EstimationType.R_6, 3d }, { Percentile.EstimationType.R_7, 2.5 },{ Percentile.EstimationType.R_8, 2.83333 }, {Percentile.EstimationType.R_9,2.81250} },
                75d, 1.0e-5);
    }

    @Test
    void testAllTechniquesLowPercentile() {
        final double[] d = new double[] { 0, 1 };
        testAssertMappedValues(d, new Object[][] { { Percentile.EstimationType.LEGACY, 0d }, { Percentile.EstimationType.R_1, 0d },
                { Percentile.EstimationType.R_2, 0d }, { Percentile.EstimationType.R_3, 0d }, { Percentile.EstimationType.R_4, 0d }, {Percentile.EstimationType.R_5, 0d}, {Percentile.EstimationType.R_6, 0d},
                { Percentile.EstimationType.R_7, 0.25 }, { Percentile.EstimationType.R_8, 0d }, {Percentile.EstimationType.R_9, 0d} },
                25d, Double.MIN_VALUE);
    }

    public void checkAllTechniquesPercentile() {
        final double[] d = new double[] { 1, 3, 2, 4 };

        testAssertMappedValues(d, new Object[][] { { Percentile.EstimationType.LEGACY, 1.5d },
                { Percentile.EstimationType.R_1, 2d }, { Percentile.EstimationType.R_2, 2d }, { Percentile.EstimationType.R_3, 1d }, { Percentile.EstimationType.R_4, 1.2 }, {Percentile.EstimationType.R_5, 1.7},
                { Percentile.EstimationType.R_6, 1.5 },{ Percentile.EstimationType.R_7, 1.9 }, { Percentile.EstimationType.R_8, 1.63333 },{ Percentile.EstimationType.R_9, 1.65 } },
                30d, 1.0e-05);

        testAssertMappedValues(d, new Object[][] { { Percentile.EstimationType.LEGACY, 1.25d },
                { Percentile.EstimationType.R_1, 1d }, { Percentile.EstimationType.R_2, 1.5d }, { Percentile.EstimationType.R_3, 1d }, { Percentile.EstimationType.R_4, 1d }, {Percentile.EstimationType.R_5, 1.5},
                { Percentile.EstimationType.R_6, 1.25 },{ Percentile.EstimationType.R_7, 1.75 },
                { Percentile.EstimationType.R_8, 1.41667 }, { Percentile.EstimationType.R_9, 1.43750 } }, 25d, 1.0e-05);

        testAssertMappedValues(d, new Object[][] { { Percentile.EstimationType.LEGACY, 3.75d },
                { Percentile.EstimationType.R_1, 3d }, { Percentile.EstimationType.R_2, 3.5d }, { Percentile.EstimationType.R_3, 3d }, { Percentile.EstimationType.R_4, 3d },
                { Percentile.EstimationType.R_5, 3.5d },{ Percentile.EstimationType.R_6, 3.75d }, { Percentile.EstimationType.R_7, 3.25 },
                { Percentile.EstimationType.R_8, 3.58333 },{ Percentile.EstimationType.R_9, 3.56250} }, 75d, 1.0e-05);

        testAssertMappedValues(d, new Object[][] { { Percentile.EstimationType.LEGACY, 2.5d },
                { Percentile.EstimationType.R_1, 2d }, { Percentile.EstimationType.R_2, 2.5d }, { Percentile.EstimationType.R_3, 2d }, { Percentile.EstimationType.R_4, 2d },
                { Percentile.EstimationType.R_5, 2.5 },{ Percentile.EstimationType.R_6, 2.5 },{ Percentile.EstimationType.R_7, 2.5 },
                { Percentile.EstimationType.R_8, 2.5 },{ Percentile.EstimationType.R_9, 2.5 } }, 50d, 1.0e-05);

        // invalid percentiles
        for (final Percentile.EstimationType e : Percentile.EstimationType.values()) {
            try {
                reset(-1.0, e);
                getUnivariateStatistic().evaluate(d, 0, d.length);
                fail();
            } catch (final MathIllegalArgumentException ex) {
                // success
            }
        }

        for (final Percentile.EstimationType e : Percentile.EstimationType.values()) {
            try {
                reset(101.0, e);
                getUnivariateStatistic().evaluate(d, 0, d.length);
                fail();
            } catch (final MathIllegalArgumentException ex) {
                // success
            }
        }
    }

    @Test
    void testAllTechniquesPercentileUsingMedianOf3Pivoting() {
        kthSelector = new KthSelector(PivotingStrategy.MEDIAN_OF_3);
        checkAllTechniquesPercentile();
    }

    @Test
    void testAllTechniquesPercentileUsingCentralPivoting() {
        kthSelector = new KthSelector(PivotingStrategy.CENTRAL);
        checkAllTechniquesPercentile();
    }

    @Test
    void testAllTechniquesNISTExample() {
        final double[] d =
                new double[] { 95.1772, 95.1567, 95.1937, 95.1959, 95.1442, 95.0610,
                               95.1591, 95.1195, 95.1772, 95.0925, 95.1990, 95.1682 };

        testAssertMappedValues(d, new Object[][] { { Percentile.EstimationType.LEGACY, 95.1981 },
                { Percentile.EstimationType.R_1, 95.19590 }, { Percentile.EstimationType.R_2, 95.19590 }, { Percentile.EstimationType.R_3, 95.19590 },
                { Percentile.EstimationType.R_4, 95.19546 }, { Percentile.EstimationType.R_5, 95.19683 }, { Percentile.EstimationType.R_6, 95.19807 },
                { Percentile.EstimationType.R_7, 95.19568 }, { Percentile.EstimationType.R_8, 95.19724 }, { Percentile.EstimationType.R_9, 95.19714 } }, 90d,
                1.0e-04);

        for (final Percentile.EstimationType e : Percentile.EstimationType.values()) {
            reset(100.0, e);
            assertEquals(95.1990, getUnivariateStatistic().evaluate(d), 1.0e-4);
        }
    }

    @Test
    void testAllTechniques5() {
        reset(5, Percentile.EstimationType.LEGACY);
        final UnivariateStatistic percentile = getUnivariateStatistic();
        assertEquals(this.percentile5, percentile.evaluate(testArray), getTolerance());
        testAssertMappedValues(testArray,
                new Object[][] { { Percentile.EstimationType.LEGACY, percentile5 }, { Percentile.EstimationType.R_1, 8.8000 },
                        { Percentile.EstimationType.R_2, 8.8000 }, { Percentile.EstimationType.R_3, 8.2000 }, { Percentile.EstimationType.R_4, 8.2600 },
                        { Percentile.EstimationType.R_5, 8.5600 }, { Percentile.EstimationType.R_6, 8.2900 },
                        { Percentile.EstimationType.R_7, 8.8100 }, { Percentile.EstimationType.R_8, 8.4700 },
                        { Percentile.EstimationType.R_9, 8.4925 }}, 5d, getTolerance());
    }

    @Test
    void testAllTechniquesNullEmpty() {

        final double[] nullArray = null;
        final double[] emptyArray = new double[] {};
        for (final Percentile.EstimationType e : Percentile.EstimationType.values()) {
            reset (50, e);
            final UnivariateStatistic percentile = getUnivariateStatistic();
            try {
                percentile.evaluate(nullArray);
                fail("Expecting NullArgumentException "
                        + "for null array");
            } catch (final NullArgumentException ex) {
                // expected
            }
            assertTrue(Double.isNaN(percentile.evaluate(emptyArray)));
        }

    }

    @Test
    void testAllTechniquesSingleton() {
        final double[] singletonArray = new double[] { 1d };
        for (final Percentile.EstimationType e : Percentile.EstimationType.values()) {
            reset (50, e);
            final UnivariateStatistic percentile = getUnivariateStatistic();
            assertEquals(1d, percentile.evaluate(singletonArray), 0);
            assertEquals(1d, percentile.evaluate(singletonArray, 0, 1), 0);
            assertEquals(1d, new Percentile().evaluate(singletonArray, 0, 1, 5), 0);
            assertEquals(1d, new Percentile().evaluate(singletonArray, 0, 1, 100), 0);
            assertTrue(Double.isNaN(percentile.evaluate(singletonArray, 0, 0)));
        }
    }

    @Test
    void testAllTechniquesEmpty() {
        final double[] singletonArray = new double[] { };
        for (final Percentile.EstimationType e : Percentile.EstimationType.values()) {
            reset (50, e);
            final UnivariateStatistic percentile = getUnivariateStatistic();
            assertEquals(Double.NaN, percentile.evaluate(singletonArray), 0);
            assertEquals(Double.NaN, percentile.evaluate(singletonArray, 0, 0), 0);
            assertEquals(Double.NaN, new Percentile().evaluate(singletonArray, 0, 0, 5), 0);
            assertEquals(Double.NaN, new Percentile().evaluate(singletonArray, 0, 0, 100), 0);
            assertTrue(Double.isNaN(percentile.evaluate(singletonArray, 0, 0)));
        }
    }

    @Test
    void testReplaceNanInRange() {
        final double[] specialValues =
            new double[] { 0d, 1d, 2d, 3d, 4d, Double.NaN, Double.NaN, 5d, 7d, Double.NaN, 8d};
        assertEquals(/*Double.NaN*/3.5, new Percentile(50d).evaluate(specialValues), 0d);
        reset (50, Percentile.EstimationType.R_1);
        assertEquals(3d, getUnivariateStatistic().evaluate(specialValues), 0d);
        reset (50, Percentile.EstimationType.R_2);
        assertEquals(3.5d, getUnivariateStatistic().evaluate(specialValues), 0d);
        assertEquals(Double.POSITIVE_INFINITY,
                     new Percentile(70).withNaNStrategy(NaNStrategy.MAXIMAL).evaluate(specialValues),
                     0d);
    }

    @Test
    void testRemoveNan() {
        final double[] specialValues = new double[] { 0d, 1d, 2d, 3d, 4d, Double.NaN };
        final double[] expectedValues = new double[] { 0d, 1d, 2d, 3d, 4d };
        reset (50, Percentile.EstimationType.R_1);
        assertEquals(2.0, getUnivariateStatistic().evaluate(specialValues), 0d);
        assertEquals(2.0, getUnivariateStatistic().evaluate(expectedValues),0d);
        assertTrue(Double.isNaN(getUnivariateStatistic().evaluate(specialValues,5,1)));
        assertEquals(4d, getUnivariateStatistic().evaluate(specialValues, 4, 2), 0d);
        assertEquals(3d, getUnivariateStatistic().evaluate(specialValues,3,3),0d);
        reset(50, Percentile.EstimationType.R_2);
        assertEquals(3.5d, getUnivariateStatistic().evaluate(specialValues,3,3),0d);
    }

    @Test
    void testPercentileCopy() {
       reset(50d, Percentile.EstimationType.LEGACY);
       final Percentile original = getUnivariateStatistic();
       final Percentile copy = new Percentile(original);
       assertEquals(original.getNaNStrategy(),copy.getNaNStrategy());
       assertEquals(original.getQuantile(), copy.getQuantile(),0d);
       assertEquals(original.getEstimationType(),copy.getEstimationType());
       assertEquals(NaNStrategy.FIXED, original.getNaNStrategy());
    }

    @Test
    void testAllTechniquesSpecialValues() {
        reset(50d, Percentile.EstimationType.LEGACY);
        final UnivariateStatistic percentile = getUnivariateStatistic();
        double[] specialValues = new double[] { 0d, 1d, 2d, 3d, 4d, Double.NaN };
        assertEquals(2.5d, percentile.evaluate(specialValues), 0);

        testAssertMappedValues(specialValues, new Object[][] {
                { Percentile.EstimationType.LEGACY, 2.5d }, { Percentile.EstimationType.R_1, 2.0 }, { Percentile.EstimationType.R_2, 2.0 }, { Percentile.EstimationType.R_3, 1.0 },
                { Percentile.EstimationType.R_4, 1.5 }, { Percentile.EstimationType.R_5, 2.0 }, { Percentile.EstimationType.R_6, 2.0 },
                { Percentile.EstimationType.R_7, 2.0 }, { Percentile.EstimationType.R_8, 2.0 }, { Percentile.EstimationType.R_9, 2.0 }}, 50d, 0d);

        specialValues =
                new double[] { Double.NEGATIVE_INFINITY, 1d, 2d, 3d, Double.NaN, Double.POSITIVE_INFINITY };
        assertEquals(2.5d, percentile.evaluate(specialValues), 0);

        testAssertMappedValues(specialValues, new Object[][] {
                { Percentile.EstimationType.LEGACY, 2.5d }, { Percentile.EstimationType.R_1, 2.0 }, { Percentile.EstimationType.R_2, 2.0 }, { Percentile.EstimationType.R_3, 1.0 },
                { Percentile.EstimationType.R_4, 1.5 }, { Percentile.EstimationType.R_5, 2.0 }, { Percentile.EstimationType.R_7, 2.0 }, { Percentile.EstimationType.R_7, 2.0 },
                { Percentile.EstimationType.R_8, 2.0 }, { Percentile.EstimationType.R_9, 2.0 } }, 50d, 0d);

        specialValues =
                new double[] { 1d, 1d, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY };
        assertTrue(Double.isInfinite(percentile.evaluate(specialValues)));

        testAssertMappedValues(specialValues, new Object[][] {
                // This is one test not matching with R results.
                { Percentile.EstimationType.LEGACY, Double.POSITIVE_INFINITY },
                { Percentile.EstimationType.R_1,/* 1.0 */Double.NaN },
                { Percentile.EstimationType.R_2, /* Double.POSITIVE_INFINITY */Double.NaN },
                { Percentile.EstimationType.R_3, /* 1.0 */Double.NaN }, { Percentile.EstimationType.R_4, /* 1.0 */Double.NaN },
                { Percentile.EstimationType.R_5, Double.POSITIVE_INFINITY },
                { Percentile.EstimationType.R_6, Double.POSITIVE_INFINITY },
                { Percentile.EstimationType.R_7, Double.POSITIVE_INFINITY },
                { Percentile.EstimationType.R_8, Double.POSITIVE_INFINITY },
                { Percentile.EstimationType.R_9, Double.POSITIVE_INFINITY }, }, 50d, 0d);

        specialValues = new double[] { 1d, 1d, Double.NaN, Double.NaN };
        assertTrue(Double.isNaN(percentile.evaluate(specialValues)));
        testAssertMappedValues(specialValues, new Object[][] {
                { Percentile.EstimationType.LEGACY, Double.NaN }, { Percentile.EstimationType.R_1, 1.0 }, { Percentile.EstimationType.R_2, 1.0 }, { Percentile.EstimationType.R_3, 1.0 },
                { Percentile.EstimationType.R_4, 1.0 }, { Percentile.EstimationType.R_5, 1.0 },{ Percentile.EstimationType.R_6, 1.0 },{ Percentile.EstimationType.R_7, 1.0 },
                { Percentile.EstimationType.R_8, 1.0 }, { Percentile.EstimationType.R_9, 1.0 },}, 50d, 0d);

        specialValues =
                new double[] { 1d, 1d, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY };

        testAssertMappedValues(specialValues, new Object[][] {
                { Percentile.EstimationType.LEGACY, Double.NaN }, { Percentile.EstimationType.R_1, Double.NaN },
                { Percentile.EstimationType.R_2, Double.NaN }, { Percentile.EstimationType.R_3, Double.NaN }, { Percentile.EstimationType.R_4, Double.NaN },
                { Percentile.EstimationType.R_5, Double.NaN }, { Percentile.EstimationType.R_6, Double.NaN },
                { Percentile.EstimationType.R_7, Double.NaN }, { Percentile.EstimationType.R_8, Double.NaN }, { Percentile.EstimationType.R_9, Double.NaN }
                }, 50d, 0d);
    }

    @Test
    void testAllTechniquesSetQuantile() {
        for (final Percentile.EstimationType e : Percentile.EstimationType.values()) {
            reset(10, e);
            final Percentile percentile = getUnivariateStatistic();
            percentile.setQuantile(100); // OK
            assertEquals(100, percentile.getQuantile(), 0);
            try {
                percentile.setQuantile(0);
                fail("Expecting MathIllegalArgumentException");
            } catch (final MathIllegalArgumentException ex) {
                // expected
            }
            try {
                new Percentile(0);
                fail("Expecting MathIllegalArgumentException");
            } catch (final MathIllegalArgumentException ex) {
                // expected
            }
        }
    }

    @Test
    void testAllTechniquesEvaluateArraySegmentWeighted() {
        for (final Percentile.EstimationType e : Percentile.EstimationType.values()) {
            reset(quantile, e);
            testEvaluateArraySegmentWeighted();
        }
    }

    @Test
    void testAllTechniquesEvaluateArraySegment() {
        for (final Percentile.EstimationType e : Percentile.EstimationType.values()) {
            reset(quantile, e);
            testEvaluateArraySegment();
        }
    }

    @Test
    void testAllTechniquesWeightedConsistency() {
        for (final Percentile.EstimationType e : Percentile.EstimationType.values()) {
            reset(quantile, e);
            testWeightedConsistency();
        }
    }

    @Test
    void testAllTechniquesEvaluation() {

        testAssertMappedValues(testArray, new Object[][] { { Percentile.EstimationType.LEGACY, 20.820 },
                { Percentile.EstimationType.R_1, 19.800 }, { Percentile.EstimationType.R_2, 19.800 }, { Percentile.EstimationType.R_3, 19.800 },
                { Percentile.EstimationType.R_4, 19.310 }, { Percentile.EstimationType.R_5, 20.280 }, { Percentile.EstimationType.R_6, 20.820 },
                { Percentile.EstimationType.R_7, 19.555 }, { Percentile.EstimationType.R_8, 20.460 },{ Percentile.EstimationType.R_9, 20.415} },
                DEFAULT_PERCENTILE, tolerance);
    }

    @Test
    void testPercentileWithTechnique() {
        reset (50, Percentile.EstimationType.LEGACY);;
        final Percentile p = getUnivariateStatistic();
        assertEquals(Percentile.EstimationType.LEGACY, p.getEstimationType());
        assertNotEquals(Percentile.EstimationType.R_1, p.getEstimationType());
    }

    static final int TINY = 10, SMALL = 50, NOMINAL = 100, MEDIUM = 500,
                     STANDARD = 1000, BIG = 10000, VERY_BIG = 50000, LARGE = 1000000,
                     VERY_LARGE = 10000000;
    static final int[] sampleSizes= { TINY, SMALL, NOMINAL, MEDIUM, STANDARD, BIG };

    @Test
    void testStoredVsDirect() {
        final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(100);
        final NormalDistribution normalDistribution = new NormalDistribution(4000, 50);
        for (final int sampleSize:sampleSizes) {
            final double[] data = randomDataGenerator.nextDeviates(normalDistribution, sampleSize);
            for (final double p:new double[] {50d,95d}) {
                for (final Percentile.EstimationType e : Percentile.EstimationType.values()) {
                    reset(p, e);
                    final Percentile pStoredData = getUnivariateStatistic();
                    pStoredData.setData(data);
                    final double storedDataResult=pStoredData.evaluate();
                    pStoredData.setData(null);
                    final Percentile pDirect = getUnivariateStatistic();
                    assertEquals(storedDataResult,
                                 pDirect.evaluate(data), 0d, "Sample="+sampleSize+",P="+p+" e="+e);
                }
            }
        }
    }

    @Test
    void testPercentileWithDataRef() {
        reset(50.0, Percentile.EstimationType.R_7);
        final Percentile p = getUnivariateStatistic();
        p.setData(testArray);
        assertEquals(Percentile.EstimationType.R_7, p.getEstimationType());
        assertNotEquals(Percentile.EstimationType.R_1, p.getEstimationType());
        assertEquals(12d, p.evaluate(), 0d);
        assertEquals(12.16d, p.evaluate(60d), 0d);
    }

    @Test
    void testNullEstimation() {
        assertThrows(NullArgumentException.class, () -> {
            type = null;
            getUnivariateStatistic();
        });
    }

    @Test
    void testAllEstimationTechniquesOnlyLimits() {
        final int N=testArray.length;

        final double[] input = testArray.clone();
        Arrays.sort(input);
        final double min = input[0];
        final double max=input[input.length-1];
        //limits may be ducked by 0.01 to induce the condition of p<pMin
        final Object[][] map =
                new Object[][] { { Percentile.EstimationType.LEGACY, 0d, 1d }, { Percentile.EstimationType.R_1, 0d, 1d },
                        { Percentile.EstimationType.R_2, 0d,1d }, { Percentile.EstimationType.R_3, 0.5/N,1d },
                        { Percentile.EstimationType.R_4, 1d/N-0.001,1d },
                        { Percentile.EstimationType.R_5, 0.5/N-0.001,(N-0.5)/N}, { Percentile.EstimationType.R_6, 0.99d/(N+1),
                            1.01d*N/(N+1)},
                        { Percentile.EstimationType.R_7, 0d,1d}, { Percentile.EstimationType.R_8, 1.99d/3/(N+1d/3),
                            (N-1d/3)/(N+1d/3)},
                        { Percentile.EstimationType.R_9, 4.99d/8/(N+0.25), (N-3d/8)/(N+0.25)} };

        for(final Object[] arr:map) {
            final Percentile.EstimationType t= (Percentile.EstimationType) arr[0];
            double pMin=(Double)arr[1];
            final double pMax=(Double)arr[2];
            assertEquals(0d, t.index(pMin, N),0d,"Type:"+t);
            assertEquals(N, t.index(pMax, N),0.5d,"Type:"+t);
            pMin=pMin==0d?pMin+0.01:pMin;
            testAssertMappedValues(testArray, new Object[][] { { t, min }}, pMin, 0.01);
            testAssertMappedValues(testArray, new Object[][] { { t, max }}, pMax * 100, tolerance);
        }
    }

    @Test
    void testAllEstimationTechniquesOnly() {
        assertEquals("Legacy Hipparchus",Percentile.EstimationType.LEGACY.getName());
        final Object[][] map =
                new Object[][] { { Percentile.EstimationType.LEGACY, 20.82 }, { Percentile.EstimationType.R_1, 19.8 },
                        { Percentile.EstimationType.R_2, 19.8 }, { Percentile.EstimationType.R_3, 19.8 }, { Percentile.EstimationType.R_4, 19.310 },
                        { Percentile.EstimationType.R_5, 20.280}, { Percentile.EstimationType.R_6, 20.820},
                        { Percentile.EstimationType.R_7, 19.555 }, { Percentile.EstimationType.R_8, 20.460 },{Percentile.EstimationType.R_9,20.415} };
        try {
            Percentile.EstimationType.LEGACY.evaluate(testArray, -1d, new KthSelector(PivotingStrategy.MEDIAN_OF_3));
        } catch (final MathIllegalArgumentException oore) {
        }
        try {
            Percentile.EstimationType.LEGACY.evaluate(testArray, 101d, new KthSelector());
        } catch (final MathIllegalArgumentException oore) {
        }
        try {
            Percentile.EstimationType.LEGACY.evaluate(testArray, 50d, new KthSelector());
        } catch(final MathIllegalArgumentException oore) {
        }
        for (final Object[] o : map) {
            final Percentile.EstimationType e = (Percentile.EstimationType) o[0];
            final double expected = (Double) o[1];
            final double result = e.evaluate(testArray, DEFAULT_PERCENTILE, new KthSelector());
            assertEquals(expected, result, tolerance, "expected[" + e + "] = " + expected +
                    " but was = " + result);
        }
    }

    @Test
    void testAllEstimationTechniquesOnlyForAllPivotingStrategies() {

        assertEquals("Legacy Hipparchus",Percentile.EstimationType.LEGACY.getName());

        for (final PivotingStrategy strategy : PivotingStrategy.values()) {
            kthSelector = new KthSelector(strategy);
            testAllEstimationTechniquesOnly();
        }
    }

    @Test
    void testAllEstimationTechniquesOnlyForExtremeIndexes() {
        final double MAX=100;
        final Object[][] map =
                new Object[][] { { Percentile.EstimationType.LEGACY, 0d, MAX}, { Percentile.EstimationType.R_1, 0d,MAX+0.5 },
                { Percentile.EstimationType.R_2, 0d,MAX}, { Percentile.EstimationType.R_3, 0d,MAX }, { Percentile.EstimationType.R_4, 0d,MAX },
                { Percentile.EstimationType.R_5, 0d,MAX }, { Percentile.EstimationType.R_6, 0d,MAX },
                { Percentile.EstimationType.R_7, 0d,MAX }, { Percentile.EstimationType.R_8, 0d,MAX }, { Percentile.EstimationType.R_9, 0d,MAX }  };
        for (final Object[] o : map) {
            final Percentile.EstimationType e = (Percentile.EstimationType) o[0];
                assertEquals(((Double)o[1]).doubleValue(), e.index(0d, (int)MAX),0d);
                assertEquals(((Double)o[2]).doubleValue(), e.index(1.0, (int)MAX),0d,"Enum:"+e);
            }
    }

    @Test
    void testAllEstimationTechniquesOnlyForNullsAndOOR() {

        final Object[][] map =
                new Object[][] { { Percentile.EstimationType.LEGACY, 20.82 }, { Percentile.EstimationType.R_1, 19.8 },
                        { Percentile.EstimationType.R_2, 19.8 }, { Percentile.EstimationType.R_3, 19.8 }, { Percentile.EstimationType.R_4, 19.310 },
                        { Percentile.EstimationType.R_5, 20.280}, { Percentile.EstimationType.R_6, 20.820},
                        { Percentile.EstimationType.R_7, 19.555 }, { Percentile.EstimationType.R_8, 20.460 },{ Percentile.EstimationType.R_9, 20.415 } };
        for (final Object[] o : map) {
            final Percentile.EstimationType e = (Percentile.EstimationType) o[0];
            try {
                e.evaluate(null, DEFAULT_PERCENTILE, new KthSelector());
                fail("Expecting NullArgumentException");
            } catch (final NullArgumentException nae) {
                // expected
            }
            try {
                e.evaluate(testArray, 120, new KthSelector());
                fail("Expecting MathIllegalArgumentException");
            } catch (final MathIllegalArgumentException oore) {
                // expected
            }
        }
    }

    /**
     * Simple test assertion utility method assuming {@link NaNStrategy default}
     * nan handling strategy specific to each {@link EstimationType type}
     *
     * @param data input data
     * @param map of expected result against a {@link EstimationType}
     * @param p the quantile to compute for
     * @param tolerance the tolerance of difference allowed
     */
    protected void testAssertMappedValues(final double[] data, final Object[][] map,
            final Double p, final Double tolerance) {
        for (final Object[] o : map) {
            final Percentile.EstimationType e = (Percentile.EstimationType) o[0];
            final double expected = (Double) o[1];
            try {
                reset(p, e);
                final double result = getUnivariateStatistic().evaluate(data);
                assertEquals(expected, result, tolerance, "expected[" + e + "] = " + expected +
                             " but was = " + result);
            } catch(final Exception ex) {
                fail("Exception occured for estimation type "+e+":"+
                     ex.getLocalizedMessage());
            }
        }
    }

    // Some NaNStrategy specific testing
    @Test
    void testNanStrategySpecific() {
        double[] specialValues = new double[] { 0d, 1d, 2d, 3d, 4d, Double.NaN };
        assertTrue(Double.isNaN(new Percentile(50d).withEstimationType(Percentile.EstimationType.LEGACY).withNaNStrategy(NaNStrategy.MAXIMAL).evaluate(specialValues, 3, 3)));
        assertEquals(2d,new Percentile(50d).withEstimationType(Percentile.EstimationType.R_1).withNaNStrategy(NaNStrategy.REMOVED).evaluate(specialValues),0d);
        assertEquals(Double.NaN,new Percentile(50d).withEstimationType(Percentile.EstimationType.R_5).withNaNStrategy(NaNStrategy.REMOVED).evaluate(new double[] {Double.NaN,Double.NaN,Double.NaN}),0d);
        assertEquals(50d,new Percentile(50d).withEstimationType(Percentile.EstimationType.R_7).withNaNStrategy(NaNStrategy.MINIMAL).evaluate(new double[] {50d,50d,50d},1,2),0d);

        specialValues = new double[] { 0d, 1d, 2d, 3d, 4d, Double.NaN, Double.NaN };
        assertEquals(3.5,new Percentile().evaluate(specialValues, 3, 4),0d);
        assertEquals(4d,new Percentile().evaluate(specialValues, 4, 3),0d);
        assertTrue(Double.isNaN(new Percentile().evaluate(specialValues, 5, 2)));

        specialValues = new double[] { 0d, 1d, 2d, 3d, 4d, Double.NaN, Double.NaN, 5d, 6d };
        assertEquals(4.5,new Percentile().evaluate(specialValues, 3, 6),0d);
        assertEquals(5d,new Percentile().evaluate(specialValues, 4, 5),0d);
        assertTrue(Double.isNaN(new Percentile().evaluate(specialValues, 5, 2)));
        assertTrue(Double.isNaN(new Percentile().evaluate(specialValues, 5, 1)));
        assertEquals(5.5,new Percentile().evaluate(specialValues, 5, 4),0d);
    }

    // Some NaNStrategy specific testing
    @Test
    void testNanStrategyFailed() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            double[] specialValues =
                new double[]{0d, 1d, 2d, 3d, 4d, Double.NaN};
            new Percentile(50d).
                withEstimationType(Percentile.EstimationType.R_9).
                withNaNStrategy(NaNStrategy.FAILED).
                evaluate(specialValues);
        });
    }

    @Test
    void testAllTechniquesSpecialValuesWithNaNStrategy() {
        double[] specialValues =
                new double[] { 0d, 1d, 2d, 3d, 4d, Double.NaN };
        try {
            new Percentile(50d).withEstimationType(Percentile.EstimationType.LEGACY).withNaNStrategy(null);
            fail("Expecting NullArgumentArgumentException "
                    + "for null Nan Strategy");
        } catch (NullArgumentException ex) {
            // expected
        }
        //This is as per each type's default NaNStrategy
        testAssertMappedValues(specialValues, new Object[][] {
                { Percentile.EstimationType.LEGACY, 2.5d }, { Percentile.EstimationType.R_1, 2.0 }, { Percentile.EstimationType.R_2, 2.0 }, { Percentile.EstimationType.R_3, 1.0 },
                { Percentile.EstimationType.R_4, 1.5 }, { Percentile.EstimationType.R_5, 2.0 }, { Percentile.EstimationType.R_6, 2.0 },
                { Percentile.EstimationType.R_7, 2.0 }, { Percentile.EstimationType.R_8, 2.0 }, { Percentile.EstimationType.R_9, 2.0 }}, 50d, 0d);

        //This is as per MAXIMAL and hence the values tend a +0.5 upward
        testAssertMappedValues(specialValues, new Object[][] {
                { Percentile.EstimationType.LEGACY, 2.5d }, { Percentile.EstimationType.R_1, 2.0 }, { Percentile.EstimationType.R_2, 2.5 }, { Percentile.EstimationType.R_3, 2.0 },
                { Percentile.EstimationType.R_4, 2.0 }, { Percentile.EstimationType.R_5, 2.5 }, { Percentile.EstimationType.R_6, 2.5 },
                { Percentile.EstimationType.R_7, 2.5 }, { Percentile.EstimationType.R_8, 2.5 }, { Percentile.EstimationType.R_9, 2.5 }}, 50d, 0d,
                NaNStrategy.MAXIMAL);

        //This is as per MINIMAL and hence the values tend a -0.5 downward
        testAssertMappedValues(specialValues, new Object[][] {
                { Percentile.EstimationType.LEGACY, 1.5d }, { Percentile.EstimationType.R_1, 1.0 }, { Percentile.EstimationType.R_2, 1.5 }, { Percentile.EstimationType.R_3, 1.0 },
                { Percentile.EstimationType.R_4, 1.0 }, { Percentile.EstimationType.R_5, 1.5 }, { Percentile.EstimationType.R_6, 1.5 },
                { Percentile.EstimationType.R_7, 1.5 }, { Percentile.EstimationType.R_8, 1.5 }, { Percentile.EstimationType.R_9, 1.5 }}, 50d, 0d,
                NaNStrategy.MINIMAL);

        //This is as per REMOVED as here Percentile.Type.CM changed its value from default
        //while rest of Estimation types were anyways defaulted to REMOVED
        testAssertMappedValues(specialValues, new Object[][] {
                { Percentile.EstimationType.LEGACY, 2.0 }, { Percentile.EstimationType.R_1, 2.0 }, { Percentile.EstimationType.R_2, 2.0 }, { Percentile.EstimationType.R_3, 1.0 },
                { Percentile.EstimationType.R_4, 1.5 }, { Percentile.EstimationType.R_5, 2.0 }, { Percentile.EstimationType.R_6, 2.0 },
                { Percentile.EstimationType.R_7, 2.0 }, { Percentile.EstimationType.R_8, 2.0 }, { Percentile.EstimationType.R_9, 2.0 }}, 50d, 0d,
                NaNStrategy.REMOVED);
    }

    /**
     * Simple test assertion utility method
     *
     * @param data input data
     * @param map of expected result against a {@link EstimationType}
     * @param p the quantile to compute for
     * @param tolerance the tolerance of difference allowed
     * @param nanStrategy NaNStrategy to be passed
     */
    protected void testAssertMappedValues(double[] data, Object[][] map,
                                          Double p, Double tolerance, NaNStrategy nanStrategy) {
        for (Object[] o : map) {
            Percentile.EstimationType e = (Percentile.EstimationType) o[0];
            double expected = (Double) o[1];
            try {
                double result = new Percentile(p).withEstimationType(e).withNaNStrategy(nanStrategy).evaluate(data);
                assertEquals(expected, result, tolerance, "expected[" + e + "] = " + expected + " but was = " + result);
            } catch(Exception ex) {
                fail("Exception occured for estimation type " + e + ":" + ex.getLocalizedMessage());
            }
        }
    }

}
