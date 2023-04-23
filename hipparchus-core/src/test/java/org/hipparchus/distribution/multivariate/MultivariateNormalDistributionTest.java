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

package org.hipparchus.distribution.multivariate;

import java.util.Random;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.random.Well19937c;
import org.hipparchus.util.Precision;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test cases for {@link MultivariateNormalDistribution}.
 */
public class MultivariateNormalDistributionTest {
    /**
     * Test the ability of the distribution to report its mean value parameter.
     */
    @Test
    public void testGetMean() {
        final double[] mu = { -1.5, 2 };
        final double[][] sigma = { { 2, -1.1 },
                                   { -1.1, 2 } };
        final MultivariateNormalDistribution d = new MultivariateNormalDistribution(mu, sigma);

        final double[] m = d.getMeans();
        for (int i = 0; i < m.length; i++) {
            Assert.assertEquals(mu[i], m[i], 0);
        }
    }

    /**
     * Test the ability of the distribution to report its covariance matrix parameter.
     */
    @Test
    public void testGetCovarianceMatrix() {
        final double[] mu = { -1.5, 2 };
        final double[][] sigma = { { 2, -1.1 },
                                   { -1.1, 2 } };
        final MultivariateNormalDistribution d = new MultivariateNormalDistribution(mu, sigma);

        final RealMatrix s = d.getCovariances();
        final int dim = d.getDimension();
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                Assert.assertEquals(sigma[i][j], s.getEntry(i, j), 0);
            }
        }
    }

    /**
     * Test the accuracy of sampling from the distribution.
     */
    @Test
    public void testSampling() {
        final double[] mu = { -1.5, 2 };
        final double[][] sigma = { { 2, -1.1 },
                                   { -1.1, 2 } };
        final MultivariateNormalDistribution d = new MultivariateNormalDistribution(mu, sigma);
        d.reseedRandomGenerator(50);

        final int n = 500000;

        final double[][] samples = d.sample(n);
        final int dim = d.getDimension();
        final double[] sampleMeans = new double[dim];

        for (int i = 0; i < samples.length; i++) {
            for (int j = 0; j < dim; j++) {
                sampleMeans[j] += samples[i][j];
            }
        }

        final double sampledValueTolerance = 1e-2;
        for (int j = 0; j < dim; j++) {
            sampleMeans[j] /= samples.length;
            Assert.assertEquals(mu[j], sampleMeans[j], sampledValueTolerance);
        }

        //final double[][] sampleSigma = new Covariance(samples).getCovarianceMatrix().getData();
        final RealMatrix sampleSigma = UnitTestUtils.covarianceMatrix(new Array2DRowRealMatrix(samples));
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                Assert.assertEquals(sigma[i][j], sampleSigma.getEntry(i, j), sampledValueTolerance);
            }
        }
    }

    /**
     * Test the accuracy of the distribution when calculating densities.
     */
    @Test
    public void testDensities() {
        final double[] mu = { -1.5, 2 };
        final double[][] sigma = { { 2, -1.1 },
                                   { -1.1, 2 } };
        final MultivariateNormalDistribution d = new MultivariateNormalDistribution(mu, sigma);

        final double[][] testValues = { { -1.5, 2 },
                                        { 4, 4 },
                                        { 1.5, -2 },
                                        { 0, 0 } };
        final double[] densities = new double[testValues.length];
        for (int i = 0; i < densities.length; i++) {
            densities[i] = d.density(testValues[i]);
        }

        // From dmvnorm function in R 2.15 CRAN package Mixtools v0.4.5
        final double[] correctDensities = { 0.09528357207691344,
                                            5.80932710124009e-09,
                                            0.001387448895173267,
                                            0.03309922090210541 };

        for (int i = 0; i < testValues.length; i++) {
            Assert.assertEquals(correctDensities[i], densities[i], 1e-16);
        }
    }

    /**
     * Test the accuracy of the distribution when calculating densities.
     */
    @Test
    public void testUnivariateDistribution() {
        final double[] mu = { -1.5 };
        final double[][] sigma = { { 1 } };

        final MultivariateNormalDistribution multi = new MultivariateNormalDistribution(mu, sigma);

        final NormalDistribution uni = new NormalDistribution(mu[0], sigma[0][0]);
        final Random rng = new Random();
        final int numCases = 100;
        final double tol = Math.ulp(1d);
        for (int i = 0; i < numCases; i++) {
            final double v = rng.nextDouble() * 10 - 5;
            Assert.assertEquals(uni.density(v), multi.density(new double[] { v }), tol);
        }
    }

    /**
     * Test getting/setting custom singularMatrixTolerance
     */
    @Test
    public void testGetSingularMatrixTolerance() {
        final double[] mu = { -1.5 };
        final double[][] sigma = { { 1 } };

        final double tolerance1 = 1e-2;
        final MultivariateNormalDistribution mvd1 = new MultivariateNormalDistribution(mu, sigma, tolerance1);
        Assert.assertEquals(tolerance1, mvd1.getSingularMatrixCheckTolerance(), Precision.EPSILON);

        final double tolerance2 = 1e-3;
        final MultivariateNormalDistribution mvd2 = new MultivariateNormalDistribution(mu, sigma, tolerance2);
        Assert.assertEquals(tolerance2, mvd2.getSingularMatrixCheckTolerance(), Precision.EPSILON);
    }

    @Test
    public void testNotPositiveDefinite() {
        try {
            new MultivariateNormalDistribution(new Well19937c(0x543l), new double[2],
                                               new double[][] { { -1.0, 0.0 }, { 0.0, -2.0 } });
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.NOT_POSITIVE_DEFINITE_MATRIX, miae.getSpecifier());
        }
    }

    @Test
    public void testStd() {
        MultivariateNormalDistribution d = new MultivariateNormalDistribution(new Well19937c(0x543l), new double[2],
                                                                              new double[][] { { 4.0, 0.0 }, { 0.0, 9.0 } });
        double[] s = d.getStandardDeviations();
        Assert.assertEquals(2, s.length);
        Assert.assertEquals(2.0, s[0], 1.0e-15);
        Assert.assertEquals(3.0, s[1], 1.0e-15);
    }

    @Test
    public void testWrongDensity() {
        try {
            MultivariateNormalDistribution d = new MultivariateNormalDistribution(new Well19937c(0x543l), new double[2],
                                                                                  new double[][] { { 4.0, 0.0 }, { 0.0, 4.0 } });
            d.density(new double[3]);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
        }
    }

    @Test
    public void testWrongArguments() {
        checkWrongArguments(new double[3], new double[6][6]);
        checkWrongArguments(new double[3], new double[3][6]);
    }

    private void checkWrongArguments(double[] means, double[][] covariances) {
        try {
            new MultivariateNormalDistribution(new Well19937c(0x543l), means, covariances);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
        }
    }
}
