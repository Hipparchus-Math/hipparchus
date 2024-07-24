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
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test cases for ExponentialDistribution.
 */
public class ExponentialDistributionTest extends RealDistributionAbstractTest {

    // --------------------- Override tolerance  --------------
    @BeforeEach
    @Override
    public void setUp() {
        super.setUp();
        setTolerance(1E-9);
    }

    //-------------- Implementations for abstract methods -----------------------

    /** Creates the default continuous distribution instance to use in tests. */
    @Override
    public ExponentialDistribution makeDistribution() {
        return new ExponentialDistribution(5.0);
    }

    /** Creates the default cumulative probability distribution test input values */
    @Override
    public double[] makeCumulativeTestPoints() {
        // quantiles computed using R version 2.9.2
        return new double[] {0.00500250166792, 0.0502516792675, 0.126589039921, 0.256466471938,
                0.526802578289, 34.5387763949, 23.0258509299, 18.4443972706, 14.9786613678, 11.5129254650};
    }

    /** Creates the default cumulative probability density test expected values */
    @Override
    public double[] makeCumulativeTestValues() {
        return new double[] {0.001, 0.01, 0.025, 0.05, 0.1, 0.999,
                0.990, 0.975, 0.950, 0.900};
    }

    /** Creates the default probability density test expected values */
    @Override
    public double[] makeDensityTestValues() {
        return new double[] {0.1998, 0.198, 0.195, 0.19, 0.18, 0.000200000000000,
                0.00200000000002, 0.00499999999997, 0.00999999999994, 0.0199999999999};
    }

    //------------ Additional tests -------------------------------------------

    @Test
    void testCumulativeProbabilityExtremes() {
        setCumulativeTestPoints(new double[] {-2, 0});
        setCumulativeTestValues(new double[] {0, 0});
        verifyCumulativeProbabilities();
    }

    @Test
    void testInverseCumulativeProbabilityExtremes() {
         setInverseCumulativeTestPoints(new double[] {0, 1});
         setInverseCumulativeTestValues(new double[] {0, Double.POSITIVE_INFINITY});
         verifyInverseCumulativeProbabilities();
    }

    @Test
    void testCumulativeProbability2() {
        double actual = getDistribution().probability(0.25, 0.75);
        assertEquals(0.0905214, actual, 10e-4);
    }

    @Test
    void testDensity() {
        ExponentialDistribution d1 = new ExponentialDistribution(1);
        assertTrue(Precision.equals(0.0, d1.density(-1e-9), 1));
        assertTrue(Precision.equals(1.0, d1.density(0.0), 1));
        assertTrue(Precision.equals(0.0, d1.density(1000.0), 1));
        assertTrue(Precision.equals(FastMath.exp(-1), d1.density(1.0), 1));
        assertTrue(Precision.equals(FastMath.exp(-2), d1.density(2.0), 1));

        ExponentialDistribution d2 = new ExponentialDistribution(3);
        assertTrue(Precision.equals(1/3.0, d2.density(0.0), 1));
        // computed using  print(dexp(1, rate=1/3), digits=10) in R 2.5
        assertEquals(0.2388437702, d2.density(1.0), 1e-8);

        // computed using  print(dexp(2, rate=1/3), digits=10) in R 2.5
        assertEquals(0.1711390397, d2.density(2.0), 1e-8);
    }

    @Test
    void testMeanAccessors() {
        ExponentialDistribution distribution = (ExponentialDistribution) getDistribution();
        assertEquals(5d, distribution.getMean(), Double.MIN_VALUE);
    }

    @Test
    void testPreconditions() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new ExponentialDistribution(0);
        });
    }

    @Test
    void testMoments() {
        final double tol = 1e-9;
        ExponentialDistribution dist;

        dist = new ExponentialDistribution(11d);
        assertEquals(11d, dist.getNumericalMean(), tol);
        assertEquals(dist.getNumericalVariance(), 11d * 11d, tol);

        dist = new ExponentialDistribution(10.5d);
        assertEquals(10.5d, dist.getNumericalMean(), tol);
        assertEquals(dist.getNumericalVariance(), 10.5d * 10.5d, tol);
    }
}
