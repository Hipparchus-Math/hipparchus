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


import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.stat.descriptive.StreamingStatistics;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

/**
 * Test cases for the TTestImpl class.
 *
 */
public class TTestTest {

    protected TTest testStatistic = new TTest();

    private double[] tooShortObs = { 1.0 };
    private double[] emptyObs = {};
    private StreamingStatistics emptyStats = new StreamingStatistics();
   StreamingStatistics tooShortStats = null;

    @BeforeEach
    public void setUp() {
        tooShortStats = new StreamingStatistics();
        tooShortStats.addValue(0d);
    }

    @Test
    public void testOneSampleT() {
        double[] observed =
            {93.0, 103.0, 95.0, 101.0, 91.0, 105.0, 96.0, 94.0, 101.0,  88.0, 98.0, 94.0, 101.0, 92.0, 95.0 };
        double mu = 100.0;
        StreamingStatistics sampleStats = null;
        sampleStats = new StreamingStatistics();
        for (int i = 0; i < observed.length; i++) {
            sampleStats.addValue(observed[i]);
        }

        // Target comparison values computed using R version 1.8.1 (Linux version)
        Assertions.assertEquals(-2.81976445346,
                testStatistic.t(mu, observed), 10E-10, "t statistic");
        Assertions.assertEquals(-2.81976445346,
                testStatistic.t(mu, sampleStats), 10E-10, "t statistic");
        Assertions.assertEquals(0.0136390585873,
                testStatistic.tTest(mu, observed), 10E-10, "p value");
        Assertions.assertEquals(0.0136390585873,
                testStatistic.tTest(mu, sampleStats), 10E-10, "p value");

        try {
            testStatistic.t(mu, (double[]) null);
            Assertions.fail("arguments too short, NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        try {
            testStatistic.t(mu, (StreamingStatistics) null);
            Assertions.fail("arguments too short, NullArgumentException expected");
        } catch (NullArgumentException ex) {
            // expected
        }

        try {
            testStatistic.t(mu, emptyObs);
            Assertions.fail("arguments too short, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.t(mu, emptyStats);
            Assertions.fail("arguments too short, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.t(mu, tooShortObs);
            Assertions.fail("insufficient data to compute t statistic, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            testStatistic.tTest(mu, tooShortObs);
            Assertions.fail("insufficient data to perform t test, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
           // expected
        }

        try {
            testStatistic.t(mu, tooShortStats);
            Assertions.fail("insufficient data to compute t statistic, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            testStatistic.tTest(mu, tooShortStats);
            Assertions.fail("insufficient data to perform t test, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testOneSampleTTest() {
        double[] oneSidedP =
            {2d, 0d, 6d, 6d, 3d, 3d, 2d, 3d, -6d, 6d, 6d, 6d, 3d, 0d, 1d, 1d, 0d, 2d, 3d, 3d };
        StreamingStatistics oneSidedPStats = new StreamingStatistics();
        for (int i = 0; i < oneSidedP.length; i++) {
            oneSidedPStats.addValue(oneSidedP[i]);
        }
        // Target comparison values computed using R version 1.8.1 (Linux version)
        Assertions.assertEquals(3.86485535541,
                testStatistic.t(0d, oneSidedP), 10E-10, "one sample t stat");
        Assertions.assertEquals(3.86485535541,
                testStatistic.t(0d, oneSidedPStats),1E-10,"one sample t stat");
        Assertions.assertEquals(0.000521637019637,
                testStatistic.tTest(0d, oneSidedP) / 2d, 10E-10, "one sample p value");
        Assertions.assertEquals(0.000521637019637,
                testStatistic.tTest(0d, oneSidedPStats) / 2d, 10E-5, "one sample p value");
        Assertions.assertTrue(testStatistic.tTest(0d, oneSidedP, 0.01), "one sample t-test reject");
        Assertions.assertTrue(testStatistic.tTest(0d, oneSidedPStats, 0.01), "one sample t-test reject");
        Assertions.assertFalse(testStatistic.tTest(0d, oneSidedP, 0.0001), "one sample t-test accept");
        Assertions.assertFalse(testStatistic.tTest(0d, oneSidedPStats, 0.0001), "one sample t-test accept");

        try {
            testStatistic.tTest(0d, oneSidedP, 95);
            Assertions.fail("alpha out of range, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.tTest(0d, oneSidedPStats, 95);
            Assertions.fail("alpha out of range, MathIllegalArgumentException expected");
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
        Assertions.assertEquals(1.60371728768,
                testStatistic.t(sample1, sample2), 1E-10, "two sample heteroscedastic t stat");
        Assertions.assertEquals(1.60371728768,
                testStatistic.t(sampleStats1, sampleStats2), 1E-10, "two sample heteroscedastic t stat");
        Assertions.assertEquals(0.128839369622,
                testStatistic.tTest(sample1, sample2), 1E-10, "two sample heteroscedastic p value");
        Assertions.assertEquals(0.128839369622,
                testStatistic.tTest(sampleStats1, sampleStats2), 1E-10, "two sample heteroscedastic p value");
        Assertions.assertTrue(testStatistic.tTest(sample1, sample2, 0.2),
                "two sample heteroscedastic t-test reject");
        Assertions.assertTrue(testStatistic.tTest(sampleStats1, sampleStats2, 0.2),
                "two sample heteroscedastic t-test reject");
        Assertions.assertFalse(testStatistic.tTest(sample1, sample2, 0.1), "two sample heteroscedastic t-test accept");
        Assertions.assertFalse(testStatistic.tTest(sampleStats1, sampleStats2, 0.1), "two sample heteroscedastic t-test accept");

        try {
            testStatistic.tTest(sample1, sample2, .95);
            Assertions.fail("alpha out of range, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.tTest(sampleStats1, sampleStats2, .95);
            Assertions.fail("alpha out of range, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.tTest(sample1, tooShortObs, .01);
            Assertions.fail("insufficient data, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.tTest(sampleStats1, tooShortStats, .01);
            Assertions.fail("insufficient data, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.tTest(sample1, tooShortObs);
            Assertions.fail("insufficient data, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
           // expected
        }

        try {
            testStatistic.tTest(sampleStats1, tooShortStats);
            Assertions.fail("insufficient data, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.t(sample1, tooShortObs);
            Assertions.fail("insufficient data, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        try {
            testStatistic.t(sampleStats1, tooShortStats);
            Assertions.fail("insufficient data, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
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
        Assertions.assertEquals(0.73096310086,
              testStatistic.homoscedasticT(sample1, sample2), 10E-11, "two sample homoscedastic t stat");
        Assertions.assertEquals(0.4833963785,
                testStatistic.homoscedasticTTest(sampleStats1, sampleStats2), 1E-10, "two sample homoscedastic p value");
        Assertions.assertTrue(testStatistic.homoscedasticTTest(sample1, sample2, 0.49),
                "two sample homoscedastic t-test reject");
        Assertions.assertFalse(testStatistic.homoscedasticTTest(sample1, sample2, 0.48), "two sample homoscedastic t-test accept");
    }

    @Test
    public void testSmallSamples() {
        double[] sample1 = {1d, 3d};
        double[] sample2 = {4d, 5d};

        // Target values computed using R, version 1.8.1 (linux version)
        Assertions.assertEquals(-2.2360679775, testStatistic.t(sample1, sample2),
                1E-10);
        Assertions.assertEquals(0.198727388935, testStatistic.tTest(sample1, sample2),
                1E-10);
    }

    @Test
    public void testPaired() {
        double[] sample1 = {1d, 3d, 5d, 7d};
        double[] sample2 = {0d, 6d, 11d, 2d};
        double[] sample3 = {5d, 7d, 8d, 10d};

        // Target values computed using R, version 1.8.1 (linux version)
        Assertions.assertEquals(-0.3133, testStatistic.pairedT(sample1, sample2), 1E-4);
        Assertions.assertEquals(0.774544295819, testStatistic.pairedTTest(sample1, sample2), 1E-10);
        Assertions.assertEquals(0.001208, testStatistic.pairedTTest(sample1, sample3), 1E-6);
        Assertions.assertFalse(testStatistic.pairedTTest(sample1, sample3, .001));
        Assertions.assertTrue(testStatistic.pairedTTest(sample1, sample3, .002));
    }
}
