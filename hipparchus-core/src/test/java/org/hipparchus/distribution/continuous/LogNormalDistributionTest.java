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

package org.hipparchus.distribution.continuous;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test cases for {@link LogNormalDistribution}.
 */
public class LogNormalDistributionTest extends RealDistributionAbstractTest {

    //-------------- Implementations for abstract methods -----------------------

    /** Creates the default real distribution instance to use in tests. */
    @Override
    public LogNormalDistribution makeDistribution() {
        return new LogNormalDistribution(2.1, 1.4);
    }

    /** Creates the default cumulative probability distribution test input values */
    @Override
    public double[] makeCumulativeTestPoints() {
        // quantiles computed using R
        return new double[] { -2.226325228634938, -1.156887023657177,
                              -0.643949578356075, -0.2027950777320613,
                              0.305827808237559, 6.42632522863494,
                              5.35688702365718, 4.843949578356074,
                              4.40279507773206, 3.89417219176244 };
    }

    /** Creates the default cumulative probability density test expected values */
    @Override
    public double[] makeCumulativeTestValues() {
        return new double[] { 0, 0, 0, 0, 0.00948199951485, 0.432056525076,
                              0.381648158697, 0.354555726206, 0.329513316888,
                              0.298422824228 };
    }

    /** Creates the default probability density test expected values */
    @Override
    public double[] makeDensityTestValues() {
        return new double[] { 0, 0, 0, 0, 0.0594218160072, 0.0436977691036,
                              0.0508364857798, 0.054873528325, 0.0587182664085,
                              0.0636229042785 };
    }

    /**
     * Creates the default inverse cumulative probability distribution test
     * input values.
     */
    @Override
    public double[] makeInverseCumulativeTestPoints() {
        // Exclude the test points less than zero, as they have cumulative
        // probability of zero, meaning the inverse returns zero, and not the
        // points less than zero.
        double[] points = makeCumulativeTestValues();
        double[] points2 = new double[points.length - 4];
        System.arraycopy(points, 4, points2, 0, points2.length - 4);
        return points2;
        //return Arrays.copyOfRange(points, 4, points.length - 4);
    }

    /**
     * Creates the default inverse cumulative probability test expected
     * values.
     */
    @Override
    public double[] makeInverseCumulativeTestValues() {
        // Exclude the test points less than zero, as they have cumulative
        // probability of zero, meaning the inverse returns zero, and not the
        // points less than zero.
        double[] points = makeCumulativeTestPoints();
        double[] points2 = new double[points.length - 4];
        System.arraycopy(points, 4, points2, 0, points2.length - 4);
        return points2;
        //return Arrays.copyOfRange(points, 1, points.length - 4);
    }

    // --------------------- Override tolerance  --------------
    @BeforeEach
    @Override
    public void setUp() {
        super.setUp();
        setTolerance(1e-9);
    }

    //---------------------------- Additional test cases -------------------------

    private void verifyQuantiles() {
        LogNormalDistribution distribution = (LogNormalDistribution)getDistribution();
        double mu = distribution.getLocation();
        double sigma = distribution.getShape();
        setCumulativeTestPoints( new double[] { mu - 2 *sigma, mu - sigma,
                                                mu, mu + sigma, mu + 2 * sigma,
                                                mu + 3 * sigma,mu + 4 * sigma,
                                                mu + 5 * sigma });
        verifyCumulativeProbabilities();
    }

    @Test
    void testQuantiles() {
        setCumulativeTestValues(new double[] {0, 0.0396495152787,
                                              0.16601209243, 0.272533253269,
                                              0.357618409638, 0.426488363093,
                                              0.483255136841, 0.530823013877});
        setDensityTestValues(new double[] {0, 0.0873055825147, 0.0847676303432,
                                           0.0677935186237, 0.0544105523058,
                                           0.0444614628804, 0.0369750288945,
                                           0.0312206409653});
        verifyQuantiles();
        verifyDensities();

        setDistribution(new LogNormalDistribution(0, 1));
        setCumulativeTestValues(new double[] {0, 0, 0, 0.5, 0.755891404214,
                                              0.864031392359, 0.917171480998,
                                              0.946239689548});
        setDensityTestValues(new double[] {0, 0, 0, 0.398942280401,
                                           0.156874019279, 0.07272825614,
                                           0.0381534565119, 0.0218507148303});
        verifyQuantiles();
        verifyDensities();

        setDistribution(new LogNormalDistribution(0, 0.1));
        setCumulativeTestValues(new double[] {0, 0, 0, 1.28417563064e-117,
                                              1.39679883412e-58,
                                              1.09839325447e-33,
                                              2.52587961726e-20,
                                              2.0824223487e-12});
        setDensityTestValues(new double[] {0, 0, 0, 2.96247992535e-114,
                                           1.1283370232e-55, 4.43812313223e-31,
                                           5.85346445002e-18,
                                           2.9446618076e-10});
        verifyQuantiles();
        verifyDensities();
    }

    @Test
    void testInverseCumulativeProbabilityExtremes() {
        setInverseCumulativeTestPoints(new double[] {0, 1});
        setInverseCumulativeTestValues(
                new double[] {0, Double.POSITIVE_INFINITY});
        verifyInverseCumulativeProbabilities();
    }

    @Test
    void testGetLocation() {
        LogNormalDistribution distribution = (LogNormalDistribution)getDistribution();
        assertEquals(2.1, distribution.getLocation(), 0);
    }

    @Test
    void testGetShape() {
        LogNormalDistribution distribution = (LogNormalDistribution)getDistribution();
        assertEquals(1.4, distribution.getShape(), 0);
    }

    @Test
    void testPreconditions() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new LogNormalDistribution(1, 0);
        });
    }

    @Test
    void testDensity() {
        double [] x = new double[]{-2, -1, 0, 1, 2};
        // R 2.13: print(dlnorm(c(-2,-1,0,1,2)), digits=10)
        checkDensity(0, 1, x, new double[] { 0.0000000000, 0.0000000000,
                                             0.0000000000, 0.3989422804,
                                             0.1568740193 });
        // R 2.13: print(dlnorm(c(-2,-1,0,1,2), mean=1.1), digits=10)
        checkDensity(1.1, 1, x, new double[] { 0.0000000000, 0.0000000000,
                                               0.0000000000, 0.2178521770,
                                               0.1836267118});
    }

    private void checkDensity(double scale, double shape, double[] x,
        double[] expected) {
        LogNormalDistribution d = new LogNormalDistribution(scale, shape);
        for (int i = 0; i < x.length; i++) {
            assertEquals(expected[i], d.density(x[i]), 1e-9);
        }
    }

    /**
     * Check to make sure top-coding of extreme values works correctly.
     * Verifies fixes for JIRA MATH-167, MATH-414
     */
    @Test
    void testExtremeValues() {
        LogNormalDistribution d = new LogNormalDistribution(0, 1);
        for (int i = 0; i < 1e5; i++) { // make sure no convergence exception
            double upperTail = d.cumulativeProbability(i);
            if (i <= 72) { // make sure not top-coded
                assertTrue(upperTail < 1.0d);
            }
            else { // make sure top coding not reversed
                assertTrue(upperTail > 0.99999);
            }
        }

        assertEquals(1, d.cumulativeProbability(Double.MAX_VALUE), 0);
        assertEquals(0, d.cumulativeProbability(-Double.MAX_VALUE), 0);
        assertEquals(1, d.cumulativeProbability(Double.POSITIVE_INFINITY), 0);
        assertEquals(0, d.cumulativeProbability(Double.NEGATIVE_INFINITY), 0);
    }

    @Test
    void testMeanVariance() {
        final double tol = 1e-9;
        LogNormalDistribution dist;

        dist = new LogNormalDistribution(0, 1);
        assertEquals(1.6487212707001282, dist.getNumericalMean(), tol);
        assertEquals(
            4.670774270471604, dist.getNumericalVariance(), tol);

        dist = new LogNormalDistribution(2.2, 1.4);
        assertEquals(24.046753552064498, dist.getNumericalMean(), tol);
        assertEquals(
            3526.913651880464, dist.getNumericalVariance(), tol);

        dist = new LogNormalDistribution(-2000.9, 10.4);
        assertEquals(0.0, dist.getNumericalMean(), tol);
        assertEquals(0.0, dist.getNumericalVariance(), tol);
    }

    @Test
    void testTinyVariance() {
        LogNormalDistribution dist = new LogNormalDistribution(0, 1e-9);
        double t = dist.getNumericalVariance();
        assertEquals(1e-18, t, 1e-20);
    }

}
