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
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for InvGammaDistribution.
 */
public class InvGammaDistributionTest extends RealDistributionAbstractTest {

    //-------------- Implementations for abstract methods -----------------------

    /** Creates the default continuous distribution instance to use in tests. */
    @Override
    public InvGammaDistribution makeDistribution() {
        return new InvGammaDistribution(4.0, 2.0);
    }
    /* Reference values for the given distribution
    0.15311308632333426 0.001                0.067305117533028799
    0.19910170259800344 0.01                 0.36987481092972235
    0.22812110266103314 0.025000000000000001 0.67229999183907063
    0.25794281611455883 0.050000000000000003 1.002262087716445
    0.29936610417768517 0.10000000000000001  1.3916447406386885
    4.6668737274551964  0.999                0.00078473049996328028
    2.4293995643995729  0.98999999999999999  0.013833804188980926
    1.8350890379656435  0.97499999999999998  0.043088414094508649
    1.4637876535641741  0.94999999999999996  0.10120309476293882
    1.1462831783710461  0.90000000000000002  0.23537771521330381
    */
    /** Creates the default cumulative probability distribution test input values */
    @Override
    public double[] makeCumulativeTestPoints() {
        // Computed with Boost/1.74.0
        return new double[] {0.15311308632333426, 0.19910170259800344, 0.2281211026610331, 0.25794281611455883,
                             0.29936610417768517, 4.6668737274551964, 2.4293995643995729, 1.8350890379656435,
                             1.4637876535641741, 1.1462831783710461};
    }

    /** Creates the default cumulative probability density test expected values */
    @Override
    public double[] makeCumulativeTestValues() {
        return new double[] {0.001, 0.01, 0.025, 0.05, 0.1, 0.999, 0.990, 0.975, 0.950, 0.900};
    }

    /** Creates the default probability density test expected values */
    @Override
    public double[] makeDensityTestValues() {
        // Computed with Boost/1.74.0
        return new double[] {0.067305117533028799, 0.36987481092972235, 0.67229999183907063, 1.002262087716445,
                             1.3916447406386885, 0.00078473049996328028, 0.013833804188980926, 0.043088414094508649,
                             0.10120309476293882, 0.23537771521330381};
    }

    // --------------------- Override tolerance  --------------
    @BeforeEach
    @Override
    public void setUp() {
        super.setUp();
        setTolerance(1e-9);
    }

    //---------------------------- Additional test cases -------------------------
    @Test
    void testParameterAccessors() {
        InvGammaDistribution distribution = (InvGammaDistribution) getDistribution();
        assertEquals(4.0, distribution.getShape(), 0);
        assertEquals(2.0, distribution.getScale(), 0);
    }

    @Test
    void testPreconditions() {
        try {
            new InvGammaDistribution(0, 1);
            fail("Expecting MathIllegalArgumentException for alpha = 0");
        } catch (MathIllegalArgumentException ex) {
            // Expected.
        }
        try {
            new InvGammaDistribution(1, 0);
            fail("Expecting MathIllegalArgumentException for alpha = 0");
        } catch (MathIllegalArgumentException ex) {
            // Expected.
        }
    }

    @Test
    void testProbabilities() {
        // Computed with Boost 1.83.0
        testProbability(0.500,  4.0, 2.0, 0.43347012036670896);
        testProbability(1.000,  4.0, 2.0, 0.85712346049854704);
        testProbability(15.50,  4.0, 2.0, 0.99998958045696096);
    }

    private void testProbability(double x, double a, double b, double expected) {
        InvGammaDistribution distribution = new InvGammaDistribution( a, b );
        double actual = distribution.cumulativeProbability(x);
        assertEquals(expected, actual, 10e-4, "probability for " + x);
    }


    @Test
    void testDensity() {
        double[] x = new double[]{0.1, 0.5, 1, 2, 5};
        // Computed with Boost 1.83.0
        checkDensity(4.0, 2.0, x, new double[]{0.00054964096598361537, 1.5629345185053167, 0.36089408863096717, 0.030656620097620192, 0.00057200643928374557});
    }

    private void checkDensity(double alpha, double beta, double[] x, double[] expected) {
        InvGammaDistribution d = new InvGammaDistribution(alpha, beta);
        for (int i = 0; i < x.length; i++) {
            assertEquals(expected[i], d.density(x[i]), 1e-5);
        }
    }

    @Test
    void testInverseCumulativeProbabilityExtremes() {
        setInverseCumulativeTestPoints(new double[] {0, 1});
        setInverseCumulativeTestValues(new double[] {0, Double.POSITIVE_INFINITY});
        verifyInverseCumulativeProbabilities();
    }

    @Test
    void testMoments() {
        final double tol = 1.e-9;
        InvGammaDistribution dist;

        double shape = 3.0;
        double scale = 2.0;
        double one   = 1.0;
        dist = new InvGammaDistribution(shape,scale);
        assertEquals(one, dist.getNumericalMean(), tol);
        assertEquals(one, dist.getNumericalVariance(), tol);

        scale        = 4.0;
        double two   = 2.0;
        double four  = 4.0;
        dist = new InvGammaDistribution(shape,scale);
        assertEquals(two,dist.getNumericalMean(),tol);
        assertEquals(four,dist.getNumericalVariance(),tol);
    }

    @Test
    public void testNegative() {
        final InvGammaDistribution dist = new InvGammaDistribution(3.0, 2.0);
        Assertions.assertEquals(0.0, dist.density(-0.5), 1.0e-15);
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, dist.logDensity(-0.5), 1.0);
        Assertions.assertEquals(0.0, dist.cumulativeProbability(-0.5), 1.0e-15);
        Assertions.assertEquals(1.0, dist.getNumericalMean(), 1.0e-15);
        Assertions.assertTrue(Double.isNaN(new InvGammaDistribution(0.999, 2.0).getNumericalMean()));
        Assertions.assertTrue(Double.isNaN(new InvGammaDistribution(1.999, 2.0).getNumericalVariance()));
    }

    @Test
    public void testInvariants() {
        final InvGammaDistribution dist = new InvGammaDistribution(3.0, 2.0);
        Assertions.assertEquals(0.0, dist.getSupportLowerBound(), 1.0e-15);
        Assertions.assertEquals(Double.POSITIVE_INFINITY, dist.getSupportUpperBound(), 1.0);
        Assertions.assertTrue(dist.isSupportConnected());

    }

}

