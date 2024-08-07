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
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for CauchyDistribution.
 */
public class CauchyDistributionTest extends RealDistributionAbstractTest {

    // --------------------- Override tolerance  --------------
    protected double defaultTolerance = 1e-9;

    @BeforeEach
    @Override
    public void setUp() {
        super.setUp();
        setTolerance(defaultTolerance);
    }

    //-------------- Implementations for abstract methods -----------------------

    /** Creates the default continuous distribution instance to use in tests. */
    @Override
    public CauchyDistribution makeDistribution() {
        return new CauchyDistribution(1.2, 2.1);
    }

    /** Creates the default cumulative probability distribution test input values */
    @Override
    public double[] makeCumulativeTestPoints() {
        // quantiles computed using R 2.9.2
        return new double[] {-667.24856187, -65.6230835029, -25.4830299460, -12.0588781808,
                -5.26313542807, 669.64856187, 68.0230835029, 27.8830299460, 14.4588781808, 7.66313542807};
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
        return new double[] {
            1.49599158008e-06, 0.000149550440335, 0.000933076881878,
            0.00370933207799,  0.0144742330437,   1.49599158008e-06,
            0.000149550440335, 0.000933076881878, 0.00370933207799,
            0.0144742330437
        };
    }

    //---------------------------- Additional test cases -------------------------

    @Test
    void testInverseCumulativeProbabilityExtremes() {
        setInverseCumulativeTestPoints(new double[] {0.0, 1.0});
        setInverseCumulativeTestValues(
                new double[] {Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY});
        verifyInverseCumulativeProbabilities();
    }

    @Test
    void testMedian() {
        CauchyDistribution distribution = (CauchyDistribution) getDistribution();
        assertEquals(1.2, distribution.getMedian(), 0.0);
    }

    @Test
    void testScale() {
        CauchyDistribution distribution = (CauchyDistribution) getDistribution();
        assertEquals(2.1, distribution.getScale(), 0.0);
    }

    @Test
    void testPreconditions() {
        try {
            new CauchyDistribution(0, 0);
            fail("Cannot have zero scale");
        } catch (MathIllegalArgumentException ex) {
            // Expected.
        }
        try {
            new CauchyDistribution(0, -1);
            fail("Cannot have negative scale");
        } catch (MathIllegalArgumentException ex) {
            // Expected.
        }
    }

    @Test
    void testMoments() {
        CauchyDistribution dist;

        dist = new CauchyDistribution(10.2, 0.15);
        assertTrue(Double.isNaN(dist.getNumericalMean()));
        assertTrue(Double.isNaN(dist.getNumericalVariance()));

        dist = new CauchyDistribution(23.12, 2.12);
        assertTrue(Double.isNaN(dist.getNumericalMean()));
        assertTrue(Double.isNaN(dist.getNumericalVariance()));
    }
}
