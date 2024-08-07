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

import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test cases for {@link NormalDistribution}.
 */
public class NormalDistributionTest extends RealDistributionAbstractTest {

    //-------------- Implementations for abstract methods -----------------------

    /** Creates the default real distribution instance to use in tests. */
    @Override
    public NormalDistribution makeDistribution() {
        return new NormalDistribution(2.1, 1.4);
    }

    /** Creates the default cumulative probability distribution test input values */
    @Override
    public double[] makeCumulativeTestPoints() {
        // quantiles computed using R
        return new double[] {-2.226325228634938d, -1.156887023657177d, -0.643949578356075d, -0.2027950777320613d, 0.305827808237559d,
                6.42632522863494d, 5.35688702365718d, 4.843949578356074d, 4.40279507773206d, 3.89417219176244d};
    }

    /** Creates the default cumulative probability density test expected values */
    @Override
    public double[] makeCumulativeTestValues() {
        return new double[] {0.001d, 0.01d, 0.025d, 0.05d, 0.1d, 0.999d,
                0.990d, 0.975d, 0.950d, 0.900d};
    }

    /** Creates the default probability density test expected values */
    @Override
    public double[] makeDensityTestValues() {
        return new double[] {0.00240506434076, 0.0190372444310, 0.0417464784322, 0.0736683145538, 0.125355951380,
                0.00240506434076, 0.0190372444310, 0.0417464784322, 0.0736683145538, 0.125355951380};
    }

    // --------------------- Override tolerance  --------------
    protected double defaultTolerance = 1e-9;

    @BeforeEach
    @Override
    public void setUp() {
        super.setUp();
        setTolerance(defaultTolerance);
    }

    //---------------------------- Additional test cases -------------------------

    private void verifyQuantiles() {
        NormalDistribution distribution = (NormalDistribution) getDistribution();
        double mu = distribution.getMean();
        double sigma = distribution.getStandardDeviation();
        setCumulativeTestPoints( new double[] {mu - 2 *sigma, mu - sigma,
                mu, mu + sigma, mu + 2 * sigma,  mu + 3 * sigma, mu + 4 * sigma,
                mu + 5 * sigma});
        // Quantiles computed using R (same as Mathematica)
        setCumulativeTestValues(new double[] {0.02275013194817921, 0.158655253931457, 0.5, 0.841344746068543,
                0.977249868051821, 0.99865010196837, 0.999968328758167,  0.999999713348428});
        verifyCumulativeProbabilities();
    }

    @Test
    void testQuantiles() {
        setDensityTestValues(new double[] {0.0385649760808, 0.172836231799, 0.284958771715, 0.172836231799, 0.0385649760808,
                0.00316560600853, 9.55930184035e-05, 1.06194251052e-06});
        verifyQuantiles();
        verifyDensities();

        setDistribution(new NormalDistribution(0, 1));
        setDensityTestValues(new double[] {0.0539909665132, 0.241970724519, 0.398942280401, 0.241970724519, 0.0539909665132,
                0.00443184841194, 0.000133830225765, 1.48671951473e-06});
        verifyQuantiles();
        verifyDensities();

        setDistribution(new NormalDistribution(0, 0.1));
        setDensityTestValues(new double[] {0.539909665132, 2.41970724519, 3.98942280401, 2.41970724519,
                0.539909665132, 0.0443184841194, 0.00133830225765, 1.48671951473e-05});
        verifyQuantiles();
        verifyDensities();
    }

    @Test
    void testInverseCumulativeProbabilityExtremes() {
        setInverseCumulativeTestPoints(new double[] {0, 1});
        setInverseCumulativeTestValues(
                new double[] {Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY});
        verifyInverseCumulativeProbabilities();
    }

    // MATH-1257
    @Test
    void testCumulativeProbability() {
        final RealDistribution dist = new NormalDistribution(0, 1);
        double x = -10;
        double expected = 7.61985e-24;
        double v = dist.cumulativeProbability(x);
        double tol = 1e-5;
        assertEquals(1, v / expected, tol);
    }

    @Test
    void testGetMean() {
        NormalDistribution distribution = (NormalDistribution) getDistribution();
        assertEquals(2.1, distribution.getMean(), 0);
    }

    @Test
    void testGetStandardDeviation() {
        NormalDistribution distribution = (NormalDistribution) getDistribution();
        assertEquals(1.4, distribution.getStandardDeviation(), 0);
    }

    @Test
    void testPreconditions() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new NormalDistribution(1, 0);
        });
    }

    @Test
    void testDensity() {
        double [] x = new double[]{-2, -1, 0, 1, 2};
        // R 2.5: print(dnorm(c(-2,-1,0,1,2)), digits=10)
        checkDensity(0, 1, x, new double[]{0.05399096651, 0.24197072452, 0.39894228040, 0.24197072452, 0.05399096651});
        // R 2.5: print(dnorm(c(-2,-1,0,1,2), mean=1.1), digits=10)
        checkDensity(1.1, 1, x, new double[]{0.003266819056,0.043983595980,0.217852177033,0.396952547477,0.266085249899});
    }

    private void checkDensity(double mean, double sd, double[] x, double[] expected) {
        NormalDistribution d = new NormalDistribution(mean, sd);
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
        NormalDistribution distribution = new NormalDistribution(0, 1);
        for (int i = 0; i < 100; i++) { // make sure no convergence exception
            double lowerTail = distribution.cumulativeProbability(-i);
            double upperTail = distribution.cumulativeProbability(i);
            if (i < 9) { // make sure not top-coded
                // For i = 10, due to bad tail precision in erf (MATH-364), 1 is returned
                // TODO: once MATH-364 is resolved, replace 9 with 30
                assertTrue(lowerTail > 0.0d);
                assertTrue(upperTail < 1.0d);
            }
            else { // make sure top coding not reversed
                assertTrue(lowerTail < 0.00001);
                assertTrue(upperTail > 0.99999);
            }
        }

        assertEquals(1, distribution.cumulativeProbability(Double.MAX_VALUE), 0);
        assertEquals(0, distribution.cumulativeProbability(-Double.MAX_VALUE), 0);
        assertEquals(1, distribution.cumulativeProbability(Double.POSITIVE_INFINITY), 0);
        assertEquals(0, distribution.cumulativeProbability(Double.NEGATIVE_INFINITY), 0);
    }

    @Test
    void testMath280() {
        NormalDistribution normal = new NormalDistribution(0,1);
        double result = normal.inverseCumulativeProbability(0.9986501019683698);
        assertEquals(3.0, result, defaultTolerance);
        result = normal.inverseCumulativeProbability(0.841344746068543);
        assertEquals(1.0, result, defaultTolerance);
        result = normal.inverseCumulativeProbability(0.9999683287581673);
        assertEquals(4.0, result, defaultTolerance);
        result = normal.inverseCumulativeProbability(0.9772498680518209);
        assertEquals(2.0, result, defaultTolerance);
    }

    @Test
    void testMoments() {
        final double tol = 1e-9;
        NormalDistribution dist;

        dist = new NormalDistribution(0, 1);
        assertEquals(0, dist.getNumericalMean(), tol);
        assertEquals(1, dist.getNumericalVariance(), tol);

        dist = new NormalDistribution(2.2, 1.4);
        assertEquals(2.2, dist.getNumericalMean(), tol);
        assertEquals(dist.getNumericalVariance(), 1.4 * 1.4, tol);

        dist = new NormalDistribution(-2000.9, 10.4);
        assertEquals(dist.getNumericalMean(), -2000.9, tol);
        assertEquals(dist.getNumericalVariance(), 10.4 * 10.4, tol);
    }
}
