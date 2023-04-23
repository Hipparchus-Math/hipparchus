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
package org.hipparchus.stat.ranking;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.random.JDKRandomGenerator;
import org.hipparchus.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;


/**
 * Test cases for NaturalRanking class
 *
 */
public class NaturalRankingTest {

    private final double[] exampleData = { 20, 17, 30, 42.3, 17, 50,
            Double.NaN, Double.NEGATIVE_INFINITY, 17 };
    private final double[] tiesFirst = { 0, 0, 2, 1, 4 };
    private final double[] tiesLast = { 4, 4, 1, 0 };
    private final double[] multipleNaNs = { 0, 1, Double.NaN, Double.NaN };
    private final double[] multipleTies = { 3, 2, 5, 5, 6, 6, 1 };
    private final double[] allSame = { 0, 0, 0, 0 };

    @Test
    public void testDefault() { // Ties averaged, NaNs failed
        NaturalRanking ranking = new NaturalRanking();
        double[] ranks;

        try {
            ranks = ranking.rank(exampleData);
            Assert.fail("expected MathIllegalArgumentException due to NaNStrategy.FAILED");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        ranks = ranking.rank(tiesFirst);
        double[] correctRanks = new double[] { 1.5, 1.5, 4, 3, 5 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesLast);
        correctRanks = new double[] { 3.5, 3.5, 2, 1 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);

        try {
            ranks = ranking.rank(multipleNaNs);
            Assert.fail("expected MathIllegalArgumentException due to NaNStrategy.FAILED");
        } catch (MathIllegalArgumentException e) {
            // expected
        }

        ranks = ranking.rank(multipleTies);
        correctRanks = new double[] { 3, 2, 4.5, 4.5, 6.5, 6.5, 1 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(allSame);
        correctRanks = new double[] { 2.5, 2.5, 2.5, 2.5 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
    }

    @Test
    public void testNaNsMaximalTiesMinimum() {
        NaturalRanking ranking = new NaturalRanking(NaNStrategy.MAXIMAL, TiesStrategy.MINIMUM);
        double[] ranks = ranking.rank(exampleData);
        double[] correctRanks = { 5, 2, 6, 7, 2, 8, 9, 1, 2 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesFirst);
        correctRanks = new double[] { 1, 1, 4, 3, 5 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesLast);
        correctRanks = new double[] { 3, 3, 2, 1 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(multipleNaNs);
        correctRanks = new double[] { 1, 2, 3, 3 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(multipleTies);
        correctRanks = new double[] { 3, 2, 4, 4, 6, 6, 1 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(allSame);
        correctRanks = new double[] { 1, 1, 1, 1 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
    }

    @Test
    public void testNaNsRemovedTiesSequential() {
        NaturalRanking ranking = new NaturalRanking(NaNStrategy.REMOVED,
                TiesStrategy.SEQUENTIAL);
        double[] ranks = ranking.rank(exampleData);
        double[] correctRanks = { 5, 2, 6, 7, 3, 8, 1, 4 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesFirst);
        correctRanks = new double[] { 1, 2, 4, 3, 5 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesLast);
        correctRanks = new double[] { 3, 4, 2, 1 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(multipleNaNs);
        correctRanks = new double[] { 1, 2 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(multipleTies);
        correctRanks = new double[] { 3, 2, 4, 5, 6, 7, 1 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(allSame);
        correctRanks = new double[] { 1, 2, 3, 4 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
    }

    @Test
    public void testNaNsMinimalTiesMaximum() {
        NaturalRanking ranking = new NaturalRanking(NaNStrategy.MINIMAL,
                TiesStrategy.MAXIMUM);
        double[] ranks = ranking.rank(exampleData);
        double[] correctRanks = { 6, 5, 7, 8, 5, 9, 2, 2, 5 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesFirst);
        correctRanks = new double[] { 2, 2, 4, 3, 5 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesLast);
        correctRanks = new double[] { 4, 4, 2, 1 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(multipleNaNs);
        correctRanks = new double[] { 3, 4, 2, 2 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(multipleTies);
        correctRanks = new double[] { 3, 2, 5, 5, 7, 7, 1 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(allSame);
        correctRanks = new double[] { 4, 4, 4, 4 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
    }

    @Test
    public void testNaNsMinimalTiesAverage() {
        NaturalRanking ranking = new NaturalRanking(NaNStrategy.MINIMAL);
        double[] ranks = ranking.rank(exampleData);
        double[] correctRanks = { 6, 4, 7, 8, 4, 9, 1.5, 1.5, 4 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesFirst);
        correctRanks = new double[] { 1.5, 1.5, 4, 3, 5 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesLast);
        correctRanks = new double[] { 3.5, 3.5, 2, 1 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(multipleNaNs);
        correctRanks = new double[] { 3, 4, 1.5, 1.5 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(multipleTies);
        correctRanks = new double[] { 3, 2, 4.5, 4.5, 6.5, 6.5, 1 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(allSame);
        correctRanks = new double[] { 2.5, 2.5, 2.5, 2.5 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
    }

    @Test
    public void testNaNsFixedTiesRandom() {
        RandomGenerator randomGenerator = new JDKRandomGenerator();
        randomGenerator.setSeed(1000);
        NaturalRanking ranking = new NaturalRanking(NaNStrategy.FIXED,
                randomGenerator);
        double[] ranks = ranking.rank(exampleData);
        double[][] correctRanks = { {5}, {2, 3, 4}, {6}, {7}, {2, 3, 4}, {8}, {Double.NaN}, {1}, {2, 3, 4} };
        UnitTestUtils.assertContains(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesFirst);
        correctRanks = new double[][] { {1, 2}, {1, 2}, {4}, {3}, {5} };
        UnitTestUtils.assertContains(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesLast);
        correctRanks = new double[][] { {3, 4}, {3, 4}, {2}, {1} };
        UnitTestUtils.assertContains(correctRanks, ranks, 0d);
        ranks = ranking.rank(multipleNaNs);
        UnitTestUtils.assertEquals(new double[] { 1, 2, Double.NaN, Double.NaN }, ranks, 0d);
        ranks = ranking.rank(multipleTies);
        correctRanks = new double[][] { {3}, {2}, {4, 5}, {4, 5}, {6, 7}, {6, 7}, {1} };
        UnitTestUtils.assertContains(correctRanks, ranks, 0d);
        ranks = ranking.rank(allSame);
        correctRanks = new double[][] { {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4} };
        UnitTestUtils.assertContains(correctRanks, ranks, 0d);
    }

    @Test
    public void testNaNsAndInfs() {
        double[] data = { 0, Double.POSITIVE_INFINITY, Double.NaN,
                Double.NEGATIVE_INFINITY };
        NaturalRanking ranking = new NaturalRanking(NaNStrategy.MAXIMAL);
        double[] ranks = ranking.rank(data);
        double[] correctRanks = new double[] { 2, 3.5, 3.5, 1 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
        ranking = new NaturalRanking(NaNStrategy.MINIMAL);
        ranks = ranking.rank(data);
        correctRanks = new double[] { 3, 4, 1.5, 1.5 };
        UnitTestUtils.assertEquals(correctRanks, ranks, 0d);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testNaNsFailed() {
        double[] data = { 0, Double.POSITIVE_INFINITY, Double.NaN, Double.NEGATIVE_INFINITY };
        NaturalRanking ranking = new NaturalRanking(NaNStrategy.FAILED);
        ranking.rank(data);
    }

    @Test
    public void testNoNaNsFailed() {
        double[] data = { 1, 2, 3, 4 };
        NaturalRanking ranking = new NaturalRanking(NaNStrategy.FAILED);
        double[] ranks = ranking.rank(data);
        UnitTestUtils.assertEquals(data, ranks, 0d);
    }

}
