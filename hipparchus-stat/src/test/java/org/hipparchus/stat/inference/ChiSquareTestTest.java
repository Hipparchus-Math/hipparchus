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
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;


/**
 * Test cases for the ChiSquareTestImpl class.
 *
 */

public class ChiSquareTestTest {

    protected ChiSquareTest testStatistic = new ChiSquareTest();

    @Test
    public void testChiSquare() {

        // Target values computed using R version 1.8.1
        // Some assembly required ;-)
        //      Use sum((obs - exp)^2/exp) for the chi-square statistic and
        //      1 - pchisq(sum((obs - exp)^2/exp), length(obs) - 1) for the p-value

        long[] observed = {10, 9, 11};
        double[] expected = {10, 10, 10};
        Assertions.assertEquals(0.2,  testStatistic.chiSquare(expected, observed), 10E-12, "chi-square statistic");
        Assertions.assertEquals(0.904837418036, testStatistic.chiSquareTest(expected, observed), 1E-10, "chi-square p-value");

        long[] observed1 = { 500, 623, 72, 70, 31 };
        double[] expected1 = { 485, 541, 82, 61, 37 };
        Assertions.assertEquals( 9.023307936427388, testStatistic.chiSquare(expected1, observed1), 1E-10, "chi-square test statistic");
        Assertions.assertEquals(0.06051952647453607, testStatistic.chiSquareTest(expected1, observed1), 1E-9, "chi-square p-value");
        Assertions.assertTrue(testStatistic.chiSquareTest(expected1, observed1, 0.08), "chi-square test reject");
        Assertions.assertFalse(testStatistic.chiSquareTest(expected1, observed1, 0.05), "chi-square test accept");

        try {
            testStatistic.chiSquareTest(expected1, observed1, 95);
            Assertions.fail("alpha out of range, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        long[] tooShortObs = { 0 };
        double[] tooShortEx = { 1 };
        try {
            testStatistic.chiSquare(tooShortEx, tooShortObs);
            Assertions.fail("arguments too short, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        // unmatched arrays
        long[] unMatchedObs = { 0, 1, 2, 3 };
        double[] unMatchedEx = { 1, 1, 2 };
        try {
            testStatistic.chiSquare(unMatchedEx, unMatchedObs);
            Assertions.fail("arrays have different lengths, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        // 0 expected count
        expected[0] = 0;
        try {
            testStatistic.chiSquareTest(expected, observed, .01);
            Assertions.fail("bad expected count, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        // negative observed count
        expected[0] = 1;
        observed[0] = -1;
        try {
            testStatistic.chiSquareTest(expected, observed, .01);
            Assertions.fail("bad expected count, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

    }

    @Test
    public void testChiSquareIndependence() {

        // Target values computed using R version 1.8.1

        long[][] counts = { {40, 22, 43}, {91, 21, 28}, {60, 10, 22}};
        Assertions.assertEquals( 22.709027688, testStatistic.chiSquare(counts), 1E-9, "chi-square test statistic");
        Assertions.assertEquals(0.000144751460134, testStatistic.chiSquareTest(counts), 1E-9, "chi-square p-value");
        Assertions.assertTrue(testStatistic.chiSquareTest(counts, 0.0002), "chi-square test reject");
        Assertions.assertFalse(testStatistic.chiSquareTest(counts, 0.0001), "chi-square test accept");

        long[][] counts2 = {{10, 15}, {30, 40}, {60, 90} };
        Assertions.assertEquals( 0.168965517241, testStatistic.chiSquare(counts2), 1E-9, "chi-square test statistic");
        Assertions.assertEquals(0.918987499852, testStatistic.chiSquareTest(counts2), 1E-9, "chi-square p-value");
        Assertions.assertFalse(testStatistic.chiSquareTest(counts2, 0.1), "chi-square test accept");

        // ragged input array
        long[][] counts3 = { {40, 22, 43}, {91, 21, 28}, {60, 10}};
        try {
            testStatistic.chiSquare(counts3);
            Assertions.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        // insufficient data
        long[][] counts4 = {{40, 22, 43}};
        try {
            testStatistic.chiSquare(counts4);
            Assertions.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        long[][] counts5 = {{40}, {40}, {30}, {10}};
        try {
            testStatistic.chiSquare(counts5);
            Assertions.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        // negative counts
        long[][] counts6 = {{10, -2}, {30, 40}, {60, 90} };
        try {
            testStatistic.chiSquare(counts6);
            Assertions.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        // bad alpha
        try {
            testStatistic.chiSquareTest(counts, 0);
            Assertions.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testChiSquareLargeTestStatistic() {
        double[] exp = new double[] {
            3389119.5, 649136.6, 285745.4, 25357364.76, 11291189.78, 543628.0,
            232921.0, 437665.75
        };

        long[] obs = new long[] {
            2372383, 584222, 257170, 17750155, 7903832, 489265, 209628, 393899
        };
        org.hipparchus.stat.inference.ChiSquareTest csti =
            new org.hipparchus.stat.inference.ChiSquareTest();
        double cst = csti.chiSquareTest(exp, obs);
        Assertions.assertEquals(0.0, cst, 1E-3, "chi-square p-value");
        Assertions.assertEquals( 114875.90421929007, testStatistic.chiSquare(exp, obs), 1E-9, "chi-square test statistic");
    }

    /** Contingency table containing zeros - PR # 32531 */
    @Test
    public void testChiSquareZeroCount() {
        // Target values computed using R version 1.8.1
        long[][] counts = { {40, 0, 4}, {91, 1, 2}, {60, 2, 0}};
        Assertions.assertEquals( 9.67444662263,
                testStatistic.chiSquare(counts), 1E-9, "chi-square test statistic");
        Assertions.assertEquals(0.0462835770603,
                testStatistic.chiSquareTest(counts), 1E-9, "chi-square p-value");
    }

    /** Target values verified using DATAPLOT version 2006.3 */
    @Test
    public void testChiSquareDataSetsComparisonEqualCounts()
        {
        long[] observed1 = {10, 12, 12, 10};
        long[] observed2 = {5, 15, 14, 10};
        Assertions.assertEquals(0.541096,
                testStatistic.chiSquareTestDataSetsComparison(
                observed1, observed2), 1E-6, "chi-square p value");
        Assertions.assertEquals(2.153846,
                testStatistic.chiSquareDataSetsComparison(
                observed1, observed2), 1E-6, "chi-square test statistic");
        Assertions.assertFalse(testStatistic.chiSquareTestDataSetsComparison(
                observed1, observed2, 0.4),
                "chi-square test result");
    }

    /** Target values verified using DATAPLOT version 2006.3 */
    @Test
    public void testChiSquareDataSetsComparisonUnEqualCounts()
        {
        long[] observed1 = {10, 12, 12, 10, 15};
        long[] observed2 = {15, 10, 10, 15, 5};
        Assertions.assertEquals(0.124115,
                testStatistic.chiSquareTestDataSetsComparison(
                observed1, observed2), 1E-6, "chi-square p value");
        Assertions.assertEquals(7.232189,
                testStatistic.chiSquareDataSetsComparison(
                observed1, observed2), 1E-6, "chi-square test statistic");
        Assertions.assertTrue(testStatistic.chiSquareTestDataSetsComparison(
                observed1, observed2, 0.13),
                "chi-square test result");
        Assertions.assertFalse(testStatistic.chiSquareTestDataSetsComparison(
                observed1, observed2, 0.12),
                "chi-square test result");
    }

    @Test
    public void testChiSquareDataSetsComparisonBadCounts()
        {
        long[] observed1 = {10, -1, 12, 10, 15};
        long[] observed2 = {15, 10, 10, 15, 5};
        try {
            testStatistic.chiSquareTestDataSetsComparison(
                    observed1, observed2);
            Assertions.fail("Expecting MathIllegalArgumentException - negative count");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        long[] observed3 = {10, 0, 12, 10, 15};
        long[] observed4 = {15, 0, 10, 15, 5};
        try {
            testStatistic.chiSquareTestDataSetsComparison(
                    observed3, observed4);
            Assertions.fail("Expecting MathIllegalArgumentException - double 0's");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        long[] observed5 = {10, 10, 12, 10, 15};
        long[] observed6 = {0, 0, 0, 0, 0};
        try {
            testStatistic.chiSquareTestDataSetsComparison(
                    observed5, observed6);
            Assertions.fail("Expecting MathIllegalArgumentException - vanishing counts");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }
}
