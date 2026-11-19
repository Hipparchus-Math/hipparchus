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
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for HalfCauchyDistribution.
 */
public class HalfCauchyDistributionTest extends RealDistributionAbstractTest {

    // --------------------- Override tolerance  --------------
    protected double defaultTolerance = 1e-9;

    @BeforeEach
    @Override
    public void setUp() {
        super.setUp();
        setTolerance(defaultTolerance);
    }

    //-------------- Implementations for abstract methods -----------------------

    /* Creates the default continuous distribution instance to use in tests. */
    /* The scale parameter is set to 1.0 to allow a comparison with the SciPy implementation */
    @Override
    public HalfCauchyDistribution makeDistribution() {
        return new HalfCauchyDistribution(1.0);
    }

    /** Creates the default cumulative probability distribution test input values */
    @Override
    public double[] makeCumulativeTestPoints() {
        // quantiles computed using Scipy
        return new double[] {1.5707976187243667e-03, 1.5709255323664916e-02, 3.9290107007669640e-02,
                             7.8701706824618439e-02, 1.5838444032453627e-01, 6.3661924876873445e+02,
                             6.3656741162871697e+01, 2.5451699579357040e+01, 1.2706204736174696e+01,
                             6.3137515146750411e+00};
    }

    /** Creates the default cumulative probability density test expected values */
    @Override
    public double[] makeCumulativeTestValues() {
        return new double[] {0.001, 0.01, 0.025, 0.05, 0.1, 0.999, 0.990, 0.975, 0.950, 0.900};
    }

    /** Creates the default probability density test expected values */
    @Override
    public double[] makeDensityTestValues() {
        // Computed with Scipy
        return new double[] {6.3661820157254645e-01, 6.3646270565375884e-01, 6.3563852921903075e-01,
                             6.3270084946368599e-01, 6.2104057764005360e-01, 1.5707950348670535e-06,
                             1.5706671382255882e-04, 9.8124314855064783e-04, 3.9189229038953119e-03,
                             1.5579194727527890e-02 };
    }

    //---------------------------- Additional test cases -------------------------

    @Test
    void testInverseCumulativeProbabilityExtremes() {
        setInverseCumulativeTestPoints(new double[] {0.0, 1.0});
        setInverseCumulativeTestValues(
                new double[] {0.0, Double.POSITIVE_INFINITY});
        verifyInverseCumulativeProbabilities();
    }

    @Test
    void testScale() {
        HalfCauchyDistribution distribution;
        distribution = new HalfCauchyDistribution(2.1);
        assertEquals(2.1, distribution.getScale(), 0.0);
    }

    @Test
    void testPreconditions() {
        try {
            new HalfCauchyDistribution(0);
            fail("Cannot have zero scale");
        } catch (MathIllegalArgumentException ex) {
            // Expected.
        }
        try {
            new HalfCauchyDistribution(-1);
            fail("Cannot have negative scale");
        } catch (MathIllegalArgumentException ex) {
            // Expected.
        }
    }

    @Test
    void testMoments() {
        HalfCauchyDistribution dist;

        dist = new HalfCauchyDistribution(0.15);
        assertTrue(Double.isNaN(dist.getNumericalMean()));
        assertTrue(Double.isNaN(dist.getNumericalVariance()));

        dist = new HalfCauchyDistribution(2.12);
        assertTrue(Double.isNaN(dist.getNumericalMean()));
        assertTrue(Double.isNaN(dist.getNumericalVariance()));
    }

    @Test
    public void testNegative() {
        final HalfCauchyDistribution dist = new HalfCauchyDistribution(0.15);
        Assertions.assertEquals(0.0, dist.density(-0.5), 1.0e-15);
        Assertions.assertEquals(Double.NEGATIVE_INFINITY, dist.logDensity(-0.5), 1.0);
        Assertions.assertEquals(0.0, dist.cumulativeProbability(-0.5), 1.0e-15);
    }

    @Test
    public void testInvariants() {
        final HalfCauchyDistribution dist = new HalfCauchyDistribution(1.0);
        Assertions.assertEquals(0.0, dist.getSupportLowerBound(), 1.0e-15);
        Assertions.assertEquals(Double.POSITIVE_INFINITY, dist.getSupportUpperBound(), 1.0);
        Assertions.assertTrue(Double.isNaN(dist.getMedian()));
        Assertions.assertTrue(dist.isSupportConnected());

    }

}

