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
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.hipparchus.TestUtils;
import org.hipparchus.stat.StatUtils;
import org.hipparchus.stat.descriptive.moment.GeometricMean;
import org.hipparchus.stat.descriptive.moment.Mean;
import org.hipparchus.stat.descriptive.moment.Variance;
import org.hipparchus.stat.descriptive.summary.Sum;
import org.hipparchus.util.FastMath;
import org.junit.Test;

/**
 * Test cases for the {@link SummaryStatistics} class.
 */
public class SummaryStatisticsTest {

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

    protected SummaryStatistics createSummaryStatistics() {
        return SummaryStatistics.create();
    }

    protected SummaryStatistics.Builder createBuilder() {
        return SummaryStatistics.builder();
    }

    /** test stats */
    @Test
    public void testStats() {
        SummaryStatistics u = createSummaryStatistics();
        assertEquals("total count",0,u.getN(),tolerance);
        u.addValue(one);
        u.addValue(twoF);
        u.addValue(twoL);
        u.addValue(three);
        assertEquals("N",n,u.getN(),tolerance);
        assertEquals("sum",sum,u.getSum(),tolerance);
        assertEquals("sumsq",sumSq,u.getSumOfSquares(),tolerance);
        assertEquals("var",var,u.getVariance(),tolerance);
        assertEquals("population var",popVar,u.getPopulationVariance(),tolerance);
        assertEquals("std",std,u.getStandardDeviation(),tolerance);
        assertEquals("mean",mean,u.getMean(),tolerance);
        assertEquals("min",min,u.getMin(),tolerance);
        assertEquals("max",max,u.getMax(),tolerance);
        u.clear();
        assertEquals("total count",0,u.getN(),tolerance);
    }

    @Test
    public void testN0andN1Conditions() {
        SummaryStatistics u = createSummaryStatistics();
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
        SummaryStatistics u = createSummaryStatistics();
        u.addValue( 1.0 );
        u.addValue( 2.0 );
        u.addValue( 3.0 );
        u.addValue( 4.0 );

        assertEquals( "Geometric mean not expected", 2.213364, u.getGeometricMean(), 0.00001 );
    }

    @Test
    public void testNaNContracts() {
        SummaryStatistics u = createSummaryStatistics();
        assertTrue("mean not NaN",Double.isNaN(u.getMean()));
        assertTrue("min not NaN",Double.isNaN(u.getMin()));
        assertTrue("std dev not NaN",Double.isNaN(u.getStandardDeviation()));
        assertTrue("var not NaN",Double.isNaN(u.getVariance()));
        assertTrue("geom mean not NaN",Double.isNaN(u.getGeometricMean()));

        u.addValue(1.0);

        assertEquals("mean not expected", 1.0, u.getMean(), Double.MIN_VALUE);
        assertEquals("variance not expected", 0.0, u.getVariance(), Double.MIN_VALUE);
        assertEquals("geometric mean not expected", 1.0, u.getGeometricMean(), Double.MIN_VALUE);

        u.addValue(-1.0);

        assertTrue("geom mean not NaN",Double.isNaN(u.getGeometricMean()));

        u.addValue(0.0);

        assertTrue("geom mean not NaN",Double.isNaN(u.getGeometricMean()));

        //FiXME: test all other NaN contract specs
    }

    @Test
    public void testGetSummary() {
        SummaryStatistics u = createSummaryStatistics();
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
        SummaryStatistics u = createSummaryStatistics();
        // Empty test
        TestUtils.checkSerializedEquality(u);
        SummaryStatistics s = (SummaryStatistics) TestUtils.serializeAndRecover(u);
        StatisticalSummary summary = s.getSummary();
        verifySummary(u, summary);

        // Add some data
        u.addValue(2d);
        u.addValue(1d);
        u.addValue(3d);
        u.addValue(4d);
        u.addValue(5d);

        // Test again
        TestUtils.checkSerializedEquality(u);
        s = (SummaryStatistics) TestUtils.serializeAndRecover(u);
        summary = s.getSummary();
        verifySummary(u, summary);
    }

    @Test
    public void testEqualsAndHashCode() {
        SummaryStatistics u = createSummaryStatistics();
        SummaryStatistics t = null;
        int emptyHash = u.hashCode();
        assertTrue("reflexive", u.equals(u));
        assertFalse("non-null compared to null", u.equals(t));
        assertFalse("wrong type", u.equals(Double.valueOf(0)));
        t = createSummaryStatistics();
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
        SummaryStatistics u = createSummaryStatistics();
        u.addValue(2d);
        u.addValue(1d);
        u.addValue(3d);
        u.addValue(4d);
        SummaryStatistics v = u.copy();
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

    private void verifySummary(SummaryStatistics u, StatisticalSummary s) {
        assertEquals("N",s.getN(),u.getN());
        TestUtils.assertEquals("sum",s.getSum(),u.getSum(),tolerance);
        TestUtils.assertEquals("var",s.getVariance(),u.getVariance(),tolerance);
        TestUtils.assertEquals("std",s.getStandardDeviation(),u.getStandardDeviation(),tolerance);
        TestUtils.assertEquals("mean",s.getMean(),u.getMean(),tolerance);
        TestUtils.assertEquals("min",s.getMin(),u.getMin(),tolerance);
        TestUtils.assertEquals("max",s.getMax(),u.getMax(),tolerance);
    }

    @Test
    public void testSetterInjection() {
        SummaryStatistics u =
                createBuilder()
                    .withMeanImpl(new Sum())
                    .withSumOfLogsImpl(new Sum())
                    .build();

        u.addValue(1);
        u.addValue(3);
        assertEquals(4, u.getMean(), 1E-14);
        assertEquals(4, u.getSumOfLogs(), 1E-14);
        assertEquals(StatUtils.geometricMean(new double[] { 1, 3}),
                     u.getGeometricMean(), 1E-14);
        u.clear();
        u.addValue(1);
        u.addValue(2);
        assertEquals(3, u.getMean(), 1E-14);
    }

    @Test
    public void testQuadraticMean() {
        final double[] values = { 1.2, 3.4, 5.6, 7.89 };
        final SummaryStatistics stats = createSummaryStatistics();

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

    /**
     * JIRA: MATH-691
     */
    @Test
    public void testOverrideVarianceWithMathClass() {
        double[] scores = {1, 2, 3, 4};
        // use "population variance"
        SummaryStatistics stats = createBuilder().withVarianceImpl(new Variance(false)).build();
        for(double i : scores) {
          stats.addValue(i);
        }
        assertEquals((new Variance(false)).evaluate(scores),stats.getVariance(), 0);
    }

    @Test
    public void testOverrideMeanWithMathClass() {
        double[] scores = {1, 2, 3, 4};
        SummaryStatistics stats = createBuilder().withMeanImpl(new Mean()).build();
        for(double i : scores) {
          stats.addValue(i);
        }
        assertEquals((new Mean()).evaluate(scores),stats.getMean(), 0);
    }

    @Test
    public void testOverrideGeoMeanWithMathClass() {
        double[] scores = {1, 2, 3, 4};
        SummaryStatistics stats = createBuilder().withGeometricMeanImpl(new GeometricMean()).build();
        for(double i : scores) {
          stats.addValue(i);
        }
        assertEquals((new GeometricMean()).evaluate(scores),stats.getGeometricMean(), 0);
    }

    @Test
    public void testToString() {
        SummaryStatistics u = createSummaryStatistics();
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
}
