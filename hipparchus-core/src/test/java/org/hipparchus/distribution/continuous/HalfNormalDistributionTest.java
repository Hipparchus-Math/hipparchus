/*
 * Licensed to the Hipparchus project under one or more
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

package org.hipparchus.distribution.continuous;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Test cases for {@link HalfNormalDistribution}.
 */
public class HalfNormalDistributionTest extends RealDistributionAbstractTest {

    //-------------- Implementations for abstract methods -----------------------

    /* Creates the default real distribution instance to use in tests. */
    /* The scale parameter is set to 1.0 to allow a comparison with the SciPy implementation */
    @Override
    public HalfNormalDistribution makeDistribution() {
        return new HalfNormalDistribution(1.0);
    }

    /** Creates the default cumulative probability distribution test input values */
    @Override
    public double[] makeCumulativeTestPoints() {
        // quantiles computed using SciPy
        return new double[] { 1.2533144654324165e-03, 1.2533469508069276e-02, 3.1337982021426479e-02,
                              6.2706777943213846e-02, 1.2566134685507416e-01, 3.2905267314919255e+00,
                              2.5758293035489004e+00, 2.2414027276049469e+00, 1.9599639845400540e+00,
                              1.6448536269514722e+00};
    }

    /** Creates the default cumulative probability density test expected values */
    @Override
    public double[] makeCumulativeTestValues() {
        return new double[] {0.001, 0.01, 0.025, 0.05, 0.1, 0.999, 0.990, 0.975, 0.950, 0.900};
    }

    /** Creates the default probability density test expected values */
    @Override
    public double[] makeDensityTestValues() {
        // Computed with SciPy
        return new double[] {0.7978839341457147, 0.7978218942756774, 0.7974928680849483,
                             0.7963174049782371, 0.7916097569752335, 0.0035543806938539,
                             0.0289194860538348, 0.0647168003177488, 0.1168901396100708,
                             0.2062712807507428};
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
        HalfNormalDistribution distribution = (HalfNormalDistribution) getDistribution();
        double mu = distribution.getMean();
        double sigma = distribution.getStandardDeviation();
        setCumulativeTestPoints( new double[] {mu - sigma, mu, mu + sigma, mu + 2 * sigma,
                                               mu + 3 * sigma, mu + 4 * sigma, mu + 5 * sigma});
        // Quantiles computed using SciPy
        setCumulativeTestValues(new double[] {0.1546652072068792, 0.5750625163166381,
                                              0.8386946523497458, 0.9548769004203032, 0.9908477853495347,
                                              0.9986686071007622, 0.9998621173822465});
        verifyCumulativeProbabilities();
    }

    @Test
    void testQuantiles() {
        verifyQuantiles();
    }


    @Test
    void testInverseCumulativeProbabilityExtremes() {
        setInverseCumulativeTestPoints(new double[] {0, 1});
        setInverseCumulativeTestValues(
                new double[] {0.0, Double.POSITIVE_INFINITY});
        verifyInverseCumulativeProbabilities();
    }

    @Test
    void testGetMean() {
        double mean;
        mean = FastMath.sqrt(2.0)/FastMath.sqrt(FastMath.PI);
        HalfNormalDistribution distribution = (HalfNormalDistribution) getDistribution();
        assertEquals(mean, distribution.getMean(), 0);
    }

    @Test
    void testGetStandardDeviation() {
        double sd;
        sd = FastMath.sqrt(1.0-2.0/FastMath.PI);
        HalfNormalDistribution distribution = (HalfNormalDistribution) getDistribution();
        assertEquals(sd, distribution.getStandardDeviation(), 0);
    }

    @Test
    void testPreconditions() {
        assertThrows(MathIllegalArgumentException.class, () -> new HalfNormalDistribution(0));
    }

    @Test
    void testDensity() {
        double [] x = new double[]{0.5, 1.0, 2.0, 4.0};
        checkDensity(1.0, x, new double[]{7.0413065352859905e-01, 4.8394144903828673e-01, 1.0798193302637613e-01,
                                          2.6766045152977074e-04});
    }

    private void checkDensity(double scale, double[] x, double[] expected) {
        HalfNormalDistribution d = new HalfNormalDistribution(scale);
        for (int i = 0; i < x.length; i++) {
            assertEquals(expected[i], d.density(x[i]), 1e-9);
        }
    }

    @Test
    void testMoments() {
        final double tol = 1e-9;
        double mean;
        double sd;
        HalfNormalDistribution dist;

        dist = new HalfNormalDistribution(1.0);
        mean = FastMath.sqrt(2.0)/FastMath.sqrt(FastMath.PI);
        sd = FastMath.sqrt(1.0-2.0/FastMath.PI);
        assertEquals(mean, dist.getNumericalMean(), tol);
        assertEquals(sd*sd, dist.getNumericalVariance(), tol);

        dist = new HalfNormalDistribution(1.4);
        mean = 1.4*FastMath.sqrt(2.0)/FastMath.sqrt(FastMath.PI);
        sd = 1.4*FastMath.sqrt(1.0-2.0/FastMath.PI);
        assertEquals(mean, dist.getNumericalMean(), tol);
        assertEquals(dist.getNumericalVariance(), sd * sd, tol);

        dist = new HalfNormalDistribution(10.4);
        mean = 10.4*FastMath.sqrt(2.0)/FastMath.sqrt(FastMath.PI);
        sd = 10.4*FastMath.sqrt(1.0-2.0/FastMath.PI);
        assertEquals(dist.getNumericalMean(), mean, tol);
        assertEquals(dist.getNumericalVariance(), sd * sd, tol);
    }

    @Test
    public void testNegative() {
        final HalfNormalDistribution dist = new HalfNormalDistribution(1.0);
        Assertions.assertEquals(0.0, dist.density(-0.5), 1.0e-15);
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, dist.logDensity(-0.5), 1.0);
        Assertions.assertEquals(0.0, dist.cumulativeProbability(-0.5), 1.0e-15);
        Assertions.assertEquals(0.0, dist.probability(-0.5, 1.0), 1.0e-15);
        Assertions.assertEquals(0.0, dist.probability(1.0, -0.5), 1.0e-15);
        Assertions.assertEquals(0.0, dist.probability(-0.5, -0.5), 1.0e-15);
    }

    @Test
    public void testInvariants() {
        final HalfNormalDistribution dist = new HalfNormalDistribution(1.0);
        Assertions.assertEquals(0.0, dist.getSupportLowerBound(), 1.0e-15);
        Assertions.assertEquals(Double.POSITIVE_INFINITY, dist.getSupportUpperBound(), 1.0);
        Assertions.assertTrue(dist.isSupportConnected());

    }

}

