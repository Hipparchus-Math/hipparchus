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
package org.hipparchus.stat.inference;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.stat.descriptive.StreamingStatistics;
import org.hipparchus.util.FastMath;
import org.junit.Test;


/**
 * Test cases for the InferenceTestUtils class.
 */
public class InferenceTestUtilsTest {

    private double[] classA = { 93.0, 103.0, 95.0, 101.0 };
    private double[] classB = { 99.0, 92.0, 102.0, 100.0, 102.0 };
    private double[] classC = { 110.0, 115.0, 111.0, 117.0, 128.0 };

    private List<double[]> classes = new ArrayList<double[]>();
    private OneWayAnova oneWayAnova = new OneWayAnova();


    @Test
    public void testChiSquare() {

        // Target values computed using R version 1.8.1
        // Some assembly required ;-)
        //      Use sum((obs - exp)^2/exp) for the chi-square statistic and
        //      1 - pchisq(sum((obs - exp)^2/exp), length(obs) - 1) for the p-value

        long[] observed = {10, 9, 11};
        double[] expected = {10, 10, 10};
        assertEquals("chi-square statistic", 0.2,
                     InferenceTestUtils.chiSquare(expected, observed), 10E-12);
        assertEquals("chi-square p-value", 0.904837418036,
                     InferenceTestUtils.chiSquareTest(expected, observed), 1E-10);

        long[] observed1 = { 500, 623, 72, 70, 31 };
        double[] expected1 = { 485, 541, 82, 61, 37 };
        assertEquals("chi-square test statistic", 9.023307936427388,
                     InferenceTestUtils.chiSquare(expected1, observed1), 1E-10);
        assertEquals("chi-square p-value", 0.06051952647453607,
                     InferenceTestUtils.chiSquareTest(expected1, observed1), 1E-9);
        assertTrue("chi-square test reject",
                   InferenceTestUtils.chiSquareTest(expected1, observed1, 0.07));
        assertTrue("chi-square test accept",
                   !InferenceTestUtils.chiSquareTest(expected1, observed1, 0.05));

        try {
            InferenceTestUtils.chiSquareTest(expected1, observed1, 95);
            fail("alpha out of range, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        long[] tooShortObs = { 0 };
        double[] tooShortEx = { 1 };
        try {
            InferenceTestUtils.chiSquare(tooShortEx, tooShortObs);
            fail("arguments too short, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        // unmatched arrays
        long[] unMatchedObs = { 0, 1, 2, 3 };
        double[] unMatchedEx = { 1, 1, 2 };
        try {
            InferenceTestUtils.chiSquare(unMatchedEx, unMatchedObs);
            fail("arrays have different lengths, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        // 0 expected count
        expected[0] = 0;
        try {
            InferenceTestUtils.chiSquareTest(expected, observed, .01);
            fail("bad expected count, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        // negative observed count
        expected[0] = 1;
        observed[0] = -1;
        try {
            InferenceTestUtils.chiSquareTest(expected, observed, .01);
            fail("bad expected count, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

    }

    @Test
    public void testChiSquareIndependence() {

        // Target values computed using R version 1.8.1

        long[][] counts = { {40, 22, 43}, {91, 21, 28}, {60, 10, 22}};
        assertEquals( "chi-square test statistic", 22.709027688, InferenceTestUtils.chiSquare(counts), 1E-9);
        assertEquals("chi-square p-value", 0.000144751460134, InferenceTestUtils.chiSquareTest(counts), 1E-9);
        assertTrue("chi-square test reject", InferenceTestUtils.chiSquareTest(counts, 0.0002));
        assertTrue("chi-square test accept", !InferenceTestUtils.chiSquareTest(counts, 0.0001));

        long[][] counts2 = {{10, 15}, {30, 40}, {60, 90} };
        assertEquals( "chi-square test statistic", 0.168965517241, InferenceTestUtils.chiSquare(counts2), 1E-9);
        assertEquals("chi-square p-value",0.918987499852, InferenceTestUtils.chiSquareTest(counts2), 1E-9);
        assertTrue("chi-square test accept", !InferenceTestUtils.chiSquareTest(counts2, 0.1));

        // ragged input array
        long[][] counts3 = { {40, 22, 43}, {91, 21, 28}, {60, 10}};
        try {
            InferenceTestUtils.chiSquare(counts3);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        // insufficient data
        long[][] counts4 = {{40, 22, 43}};
        try {
            InferenceTestUtils.chiSquare(counts4);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        long[][] counts5 = {{40}, {40}, {30}, {10}};
        try {
            InferenceTestUtils.chiSquare(counts5);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        // negative counts
        long[][] counts6 = {{10, -2}, {30, 40}, {60, 90} };
        try {
            InferenceTestUtils.chiSquare(counts6);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        // bad alpha
        try {
            InferenceTestUtils.chiSquareTest(counts, 0);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testChiSquareLargeTestStatistic() {
        double[] exp = new double[] {
                3389119.5, 649136.6, 285745.4, 25357364.76,
                11291189.78, 543628.0, 232921.0, 437665.75
        };

        long[] obs = new long[] { 2372383, 584222, 257170, 17750155, 7903832, 489265, 209628, 393899 };

        ChiSquareTest csti = new ChiSquareTest();
        double cst = csti.chiSquareTest(exp, obs);
        assertEquals("chi-square p-value", 0.0, cst, 1E-3);
        assertEquals("chi-square test statistic", 114875.90421929007,
                     InferenceTestUtils.chiSquare(exp, obs), 1E-9);
    }

    /** Contingency table containing zeros - PR # 32531 */
    @Test
    public void testChiSquareZeroCount() {
        // Target values computed using R version 1.8.1
        long[][] counts = { {40, 0, 4}, {91, 1, 2}, {60, 2, 0}};
        assertEquals( "chi-square test statistic", 9.67444662263, InferenceTestUtils.chiSquare(counts), 1E-9);
        assertEquals("chi-square p-value", 0.0462835770603, InferenceTestUtils.chiSquareTest(counts), 1E-9);
    }

    private double[] tooShortObs = { 1.0 };
    private double[] emptyObs = {};
    private StreamingStatistics emptyStats = new StreamingStatistics();

    @Test
    public void testOneSampleT() {
        double[] observed = {
            93.0, 103.0, 95.0, 101.0, 91.0, 105.0, 96.0, 94.0,
            101.0, 88.0, 98.0, 94.0, 101.0, 92.0, 95.0
        };
        double mu = 100.0;
        StreamingStatistics sampleStats = new StreamingStatistics();
        for (int i = 0; i < observed.length; i++) {
            sampleStats.addValue(observed[i]);
        }

        // Target comparison values computed using R version 1.8.1 (Linux version)
        assertEquals("t statistic",  -2.81976445346, InferenceTestUtils.t(mu, observed), 10E-10);
        assertEquals("t statistic",  -2.81976445346, InferenceTestUtils.t(mu, sampleStats), 10E-10);
        assertEquals("p value", 0.0136390585873, InferenceTestUtils.tTest(mu, observed), 10E-10);
        assertEquals("p value", 0.0136390585873, InferenceTestUtils.tTest(mu, sampleStats), 10E-10);

        try {
            InferenceTestUtils.t(mu, (double[]) null);
            fail("arguments too short, NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        try {
            InferenceTestUtils.t(mu, (StreamingStatistics) null);
            fail("arguments too short, NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        try {
            InferenceTestUtils.t(mu, emptyObs);
            fail("arguments too short, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            InferenceTestUtils.t(mu, emptyStats);
            fail("arguments too short, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            InferenceTestUtils.t(mu, tooShortObs);
            fail("insufficient data to compute t statistic, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            InferenceTestUtils.tTest(mu, tooShortObs);
            fail("insufficient data to perform t test, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            InferenceTestUtils.t(mu, (StreamingStatistics) null);
            fail("insufficient data to compute t statistic, NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }
        try {
            InferenceTestUtils.tTest(mu, (StreamingStatistics) null);
            fail("insufficient data to perform t test, NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testOneSampleTTest() {
        double[] oneSidedP = {
            2d, 0d, 6d, 6d, 3d, 3d, 2d, 3d, -6d, 6d,
            6d, 6d, 3d, 0d, 1d, 1d, 0d, 2d, 3d, 3d
        };
        StreamingStatistics oneSidedPStats = new StreamingStatistics();
        for (int i = 0; i < oneSidedP.length; i++) {
            oneSidedPStats.addValue(oneSidedP[i]);
        }
        // Target comparison values computed using R version 1.8.1 (Linux version)
        assertEquals("one sample t stat", 3.86485535541, InferenceTestUtils.t(0d, oneSidedP), 10E-10);
        assertEquals("one sample t stat", 3.86485535541, InferenceTestUtils.t(0d, oneSidedPStats),1E-10);
        assertEquals("one sample p value", 0.000521637019637, InferenceTestUtils.tTest(0d, oneSidedP) / 2d, 10E-10);
        assertEquals("one sample p value", 0.000521637019637, InferenceTestUtils.tTest(0d, oneSidedPStats) / 2d, 10E-5);
        assertTrue("one sample t-test reject", InferenceTestUtils.tTest(0d, oneSidedP, 0.01));
        assertTrue("one sample t-test reject", InferenceTestUtils.tTest(0d, oneSidedPStats, 0.01));
        assertTrue("one sample t-test accept", !InferenceTestUtils.tTest(0d, oneSidedP, 0.0001));
        assertTrue("one sample t-test accept", !InferenceTestUtils.tTest(0d, oneSidedPStats, 0.0001));

        try {
            InferenceTestUtils.tTest(0d, oneSidedP, 95);
            fail("alpha out of range, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            InferenceTestUtils.tTest(0d, oneSidedPStats, 95);
            fail("alpha out of range, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

    }

    @Test
    public void testTwoSampleTHeterscedastic() {
        double[] sample1 = { 7d, -4d, 18d, 17d, -3d, -5d, 1d, 10d, 11d, -2d };
        double[] sample2 = { -1d, 12d, -1d, -3d, 3d, -5d, 5d, 2d, -11d, -1d, -3d };
        StreamingStatistics sampleStats1 = new StreamingStatistics();
        for (int i = 0; i < sample1.length; i++) {
            sampleStats1.addValue(sample1[i]);
        }
        StreamingStatistics sampleStats2 = new StreamingStatistics();
        for (int i = 0; i < sample2.length; i++) {
            sampleStats2.addValue(sample2[i]);
        }

        // Target comparison values computed using R version 1.8.1 (Linux version)
        assertEquals("two sample heteroscedastic t stat", 1.60371728768,
                     InferenceTestUtils.t(sample1, sample2), 1E-10);
        assertEquals("two sample heteroscedastic t stat", 1.60371728768,
                     InferenceTestUtils.t(sampleStats1, sampleStats2), 1E-10);
        assertEquals("two sample heteroscedastic p value", 0.128839369622,
                     InferenceTestUtils.tTest(sample1, sample2), 1E-10);
        assertEquals("two sample heteroscedastic p value", 0.128839369622,
                     InferenceTestUtils.tTest(sampleStats1, sampleStats2), 1E-10);
        assertTrue("two sample heteroscedastic t-test reject",
                   InferenceTestUtils.tTest(sample1, sample2, 0.2));
        assertTrue("two sample heteroscedastic t-test reject",
                   InferenceTestUtils.tTest(sampleStats1, sampleStats2, 0.2));
        assertTrue("two sample heteroscedastic t-test accept",
                   !InferenceTestUtils.tTest(sample1, sample2, 0.1));
        assertTrue("two sample heteroscedastic t-test accept",
                   !InferenceTestUtils.tTest(sampleStats1, sampleStats2, 0.1));

        try {
            InferenceTestUtils.tTest(sample1, sample2, .95);
            fail("alpha out of range, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            InferenceTestUtils.tTest(sampleStats1, sampleStats2, .95);
            fail("alpha out of range, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            InferenceTestUtils.tTest(sample1, tooShortObs, .01);
            fail("insufficient data, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            InferenceTestUtils.tTest(sampleStats1, (StreamingStatistics) null, .01);
            fail("insufficient data, NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        try {
            InferenceTestUtils.tTest(sample1, tooShortObs);
            fail("insufficient data, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            InferenceTestUtils.tTest(sampleStats1, (StreamingStatistics) null);
            fail("insufficient data, NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        try {
            InferenceTestUtils.t(sample1, tooShortObs);
            fail("insufficient data, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            InferenceTestUtils.t(sampleStats1, (StreamingStatistics) null);
            fail("insufficient data, NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }
    }
    @Test
    public void testTwoSampleTHomoscedastic() {
        double[] sample1 ={2, 4, 6, 8, 10, 97};
        double[] sample2 = {4, 6, 8, 10, 16};
        StreamingStatistics sampleStats1 = new StreamingStatistics();
        for (int i = 0; i < sample1.length; i++) {
            sampleStats1.addValue(sample1[i]);
        }
        StreamingStatistics sampleStats2 = new StreamingStatistics();
        for (int i = 0; i < sample2.length; i++) {
            sampleStats2.addValue(sample2[i]);
        }

        // Target comparison values computed using R version 1.8.1 (Linux version)
        assertEquals("two sample homoscedastic t stat", 0.73096310086,
                     InferenceTestUtils.homoscedasticT(sample1, sample2), 10E-11);
        assertEquals("two sample homoscedastic p value", 0.4833963785,
                     InferenceTestUtils.homoscedasticTTest(sampleStats1, sampleStats2), 1E-10);
        assertTrue("two sample homoscedastic t-test reject",
                   InferenceTestUtils.homoscedasticTTest(sample1, sample2, 0.49));
        assertTrue("two sample homoscedastic t-test accept",
                   !InferenceTestUtils.homoscedasticTTest(sample1, sample2, 0.48));
    }

    @Test
    public void testSmallSamples() {
        double[] sample1 = {1d, 3d};
        double[] sample2 = {4d, 5d};

        // Target values computed using R, version 1.8.1 (linux version)
        assertEquals(-2.2360679775, InferenceTestUtils.t(sample1, sample2), 1E-10);
        assertEquals(0.198727388935, InferenceTestUtils.tTest(sample1, sample2), 1E-10);
    }

    @Test
    public void testPaired() {
        double[] sample1 = {1d, 3d, 5d, 7d};
        double[] sample2 = {0d, 6d, 11d, 2d};
        double[] sample3 = {5d, 7d, 8d, 10d};

        // Target values computed using R, version 1.8.1 (linux version)
        assertEquals(-0.3133, InferenceTestUtils.pairedT(sample1, sample2), 1E-4);
        assertEquals(0.774544295819, InferenceTestUtils.pairedTTest(sample1, sample2), 1E-10);
        assertEquals(0.001208, InferenceTestUtils.pairedTTest(sample1, sample3), 1E-6);
        assertFalse(InferenceTestUtils.pairedTTest(sample1, sample3, .001));
        assertTrue(InferenceTestUtils.pairedTTest(sample1, sample3, .002));
    }

    @Test
    public void testOneWayAnovaUtils() {
        classes.add(classA);
        classes.add(classB);
        classes.add(classC);
        assertEquals(oneWayAnova.anovaFValue(classes),
                     InferenceTestUtils.oneWayAnovaFValue(classes), 10E-12);
        assertEquals(oneWayAnova.anovaPValue(classes),
                     InferenceTestUtils.oneWayAnovaPValue(classes), 10E-12);
        assertEquals(oneWayAnova.anovaTest(classes, 0.01),
                     InferenceTestUtils.oneWayAnovaTest(classes, 0.01));
    }
    @Test
    public void testGTestGoodnesOfFit() throws Exception {
        double[] exp = new double[] { 0.54d, 0.40d, 0.05d, 0.01d };
        long[] obs = new long[] { 70, 79, 3, 4 };

        assertEquals("G test statistic", 13.144799, InferenceTestUtils.g(exp, obs), 1E-5);
        double p_gtgf = InferenceTestUtils.gTest(exp, obs);
        assertEquals("g-Test p-value", 0.004333, p_gtgf, 1E-5);
        assertTrue(InferenceTestUtils.gTest(exp, obs, 0.05));
    }

    @Test
    public void testGTestIndependance() throws Exception {
        long[] obs1 = new long[] { 268, 199, 42 };
        long[] obs2 = new long[] { 807, 759, 184 };

        double g = InferenceTestUtils.gDataSetsComparison(obs1, obs2);

        assertEquals("G test statistic", 7.3008170, g, 1E-4);
        double p_gti = InferenceTestUtils.gTestDataSetsComparison(obs1, obs2);

        assertEquals("g-Test p-value", 0.0259805, p_gti, 1E-4);
        assertTrue(InferenceTestUtils.gTestDataSetsComparison(obs1, obs2, 0.05));
    }

    @Test
    public void testRootLogLikelihood() {
        // positive where k11 is bigger than expected.
        assertTrue(InferenceTestUtils.rootLogLikelihoodRatio(904, 21060, 1144, 283012) > 0.0);

        // negative because k11 is lower than expected
        assertTrue(InferenceTestUtils.rootLogLikelihoodRatio(36, 21928, 60280, 623876) < 0.0);

        assertEquals(FastMath.sqrt(2.772589), InferenceTestUtils.rootLogLikelihoodRatio(1, 0, 0, 1), 0.000001);
        assertEquals(-FastMath.sqrt(2.772589), InferenceTestUtils.rootLogLikelihoodRatio(0, 1, 1, 0), 0.000001);
        assertEquals(FastMath.sqrt(27.72589), InferenceTestUtils.rootLogLikelihoodRatio(10, 0, 0, 10), 0.00001);

        assertEquals(FastMath.sqrt(39.33052), InferenceTestUtils.rootLogLikelihoodRatio(5, 1995, 0, 100000), 0.00001);
        assertEquals(-FastMath.sqrt(39.33052), InferenceTestUtils.rootLogLikelihoodRatio(0, 100000, 5, 1995), 0.00001);

        assertEquals(FastMath.sqrt(4730.737), InferenceTestUtils.rootLogLikelihoodRatio(1000, 1995, 1000, 100000), 0.001);
        assertEquals(-FastMath.sqrt(4730.737), InferenceTestUtils.rootLogLikelihoodRatio(1000, 100000, 1000, 1995), 0.001);

        assertEquals(FastMath.sqrt(5734.343), InferenceTestUtils.rootLogLikelihoodRatio(1000, 1000, 1000, 100000), 0.001);
        assertEquals(FastMath.sqrt(5714.932), InferenceTestUtils.rootLogLikelihoodRatio(1000, 1000, 1000, 99000), 0.001);
    }

    @Test
    public void testKSOneSample() throws Exception {
       final NormalDistribution unitNormal = new NormalDistribution(0d, 1d);
       final double[] sample = KolmogorovSmirnovTestTest.gaussian;
       final double tol = KolmogorovSmirnovTestTest.TOLERANCE;
       assertEquals(0.3172069207622391, InferenceTestUtils.kolmogorovSmirnovTest(unitNormal, sample), tol);
       assertEquals(0.0932947561266756, InferenceTestUtils.kolmogorovSmirnovStatistic(unitNormal, sample), tol);
    }

    @Test
    public void testKSTwoSample() throws Exception {
        final double tol = KolmogorovSmirnovTestTest.TOLERANCE;
        final double[] smallSample1 = { 6, 7, 9, 13, 19, 21, 22, 23, 24 };
        final double[] smallSample2 = { 10, 11, 12, 16, 20, 27, 28, 32, 44, 54 };

        assertEquals(0.105577085453247,
                     InferenceTestUtils.kolmogorovSmirnovTest(smallSample1, smallSample2, false), tol);
        final double d = InferenceTestUtils.kolmogorovSmirnovStatistic(smallSample1, smallSample2);
        assertEquals(0.5, d, tol);
        assertEquals(0.105577085453247,
                     InferenceTestUtils.exactP(d, smallSample1.length,smallSample2.length, false), tol);
    }
}
