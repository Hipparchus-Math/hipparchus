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

import org.hipparchus.distribution.discrete.BinomialDistribution;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test cases for the BinomialTest class.
 */
public class BinomialTestTest {

    protected BinomialTest testStatistic = new BinomialTest();

    private static int successes = 51;
    private static int trials = 235;
    private static double probability = 1.0 / 6.0;

    @Test
    public void testBinomialTestPValues() {
        Assert.assertEquals(0.04375, testStatistic.binomialTest(
            trials, successes, probability, AlternativeHypothesis.TWO_SIDED), 1E-4);
        Assert.assertEquals(0.02654, testStatistic.binomialTest(
            trials, successes, probability, AlternativeHypothesis.GREATER_THAN), 1E-4);
        Assert.assertEquals(0.982, testStatistic.binomialTest(
            trials, successes, probability, AlternativeHypothesis.LESS_THAN), 1E-4);
    }

    @Test
    public void testBinomialTestExceptions() {
        try {
            testStatistic.binomialTest(10, -1, 0.5, AlternativeHypothesis.TWO_SIDED);
            Assert.fail("Expected not positive exception");
        } catch (MathIllegalArgumentException e) {
            // expected exception;
        }

        try {
            testStatistic.binomialTest(10, 11, 0.5, AlternativeHypothesis.TWO_SIDED);
            Assert.fail("Expected illegal argument exception");
        } catch (MathIllegalArgumentException e) {
            // expected exception;
        }
        try {
            testStatistic.binomialTest(10, 11, 0.5, null);
            Assert.fail("Expected illegal argument exception");
        } catch (MathIllegalArgumentException e) {
            // expected exception;
        }
    }

    @Test
    public void testBinomialTestAcceptReject() {
        double alpha05 = 0.05;
        double alpha01 = 0.01;

        Assert.assertTrue(testStatistic.binomialTest(trials, successes, probability, AlternativeHypothesis.TWO_SIDED, alpha05));
        Assert.assertTrue(testStatistic.binomialTest(trials, successes, probability, AlternativeHypothesis.GREATER_THAN, alpha05));
        Assert.assertFalse(testStatistic.binomialTest(trials, successes, probability, AlternativeHypothesis.LESS_THAN, alpha05));

        Assert.assertFalse(testStatistic.binomialTest(trials, successes, probability, AlternativeHypothesis.TWO_SIDED, alpha01));
        Assert.assertFalse(testStatistic.binomialTest(trials, successes, probability, AlternativeHypothesis.GREATER_THAN, alpha01));
        Assert.assertFalse(testStatistic.binomialTest(trials, successes, probability, AlternativeHypothesis.LESS_THAN, alpha05));
    }

    /**
     * All successes with p &gt;&gt; 0.5 - p-value picks up all mass points.
     */
    @Test
    public void testAllSuccessesTwoSidedHighP() {
        Assert.assertEquals(1d, testStatistic.binomialTest(200, 200, 0.9950429, AlternativeHypothesis.TWO_SIDED),
                            Double.MIN_VALUE);
    }

    /**
     * All successes with p = 0.5 - p-value is the sum of the two tails.
     */
    @Test
    public void testAllSuccessesTwoSidedEvenP() {
        Assert.assertEquals(2 * FastMath.pow(0.5, 5),
                            testStatistic.binomialTest(5, 5, 0.5,
                            AlternativeHypothesis.TWO_SIDED),
                            Double.MIN_VALUE);
    }

    /**
     * All successes with p = 0.5 - p-value is the sum of the two tails.
     */
    @Test
    public void testNoSuccessesTwoSidedEvenP() {
        Assert.assertEquals(2 * FastMath.pow(0.5, 5),
                            testStatistic.binomialTest(5, 0, 0.5,
                            AlternativeHypothesis.TWO_SIDED),
                            Double.MIN_VALUE);
    }

    /**
     * All successes with p &lt; 0.5 - p-value is 5 mass point.
     */
    @Test
    public void testAllSuccessesTwoSidedLowP() {
        final BinomialDistribution dist = new BinomialDistribution(5, 0.4);
        Assert.assertEquals(dist.probability(5),
                            testStatistic.binomialTest(5, 5, 0.4,
                            AlternativeHypothesis.TWO_SIDED),
                            Double.MIN_VALUE);
    }

    @Test
    /**
     * No successes, p > 0.5 - p-value is 0 mass point.
     */
    public void testNoSuccessesTwoSidedHighP() {
        final BinomialDistribution dist = new BinomialDistribution(5, 0.9);
        Assert.assertEquals(dist.probability(0),
                            testStatistic.binomialTest(5, 0, 0.9,
                            AlternativeHypothesis.TWO_SIDED),
                            Double.MIN_VALUE);
    }


    /**
     * In this case, the distribution looks like this:
     *    0: 0.32768
     *    1: 0.4096
     *    2: 0.2048
     *    3: 0.0512
     *    4: 0.0064
     *    5: 3.2E-4
     *  Algorithm picks up 5, 4, 3, 2 and then 0, so result is 1 - mass at 1.
     */
    @Test
    public void testNoSuccessesTwoSidedLowP() {
        final BinomialDistribution dist = new BinomialDistribution(5, 0.2);
        Assert.assertEquals(1 - dist.probability(1),
                            testStatistic.binomialTest(5, 0, 0.2,
                            AlternativeHypothesis.TWO_SIDED),
                            Double.MIN_VALUE);
    }

    /**
     * No successes has highest mass, so end up with everything here.
     */
    @Test
    public void testNoSuccessesTwoSidedVeryLowP() {
        Assert.assertEquals(1d,
                            testStatistic.binomialTest(5, 0,  0.001,
                            AlternativeHypothesis.TWO_SIDED),
                            Double.MIN_VALUE);
    }
}
