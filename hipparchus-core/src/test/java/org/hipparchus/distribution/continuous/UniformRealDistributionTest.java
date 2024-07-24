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

/**
 * Test cases for UniformRealDistribution.
 */
public class UniformRealDistributionTest extends RealDistributionAbstractTest {

    // --- Override tolerance -------------------------------------------------

    @BeforeEach
    @Override
    public void setUp() {
        super.setUp();
        setTolerance(1e-4);
    }

    //--- Implementations for abstract methods --------------------------------

    /** Creates the default uniform real distribution instance to use in tests. */
    @Override
    public UniformRealDistribution makeDistribution() {
        return new UniformRealDistribution(-0.5, 1.25);
    }

    /** Creates the default cumulative probability distribution test input values */
    @Override
    public double[] makeCumulativeTestPoints() {
        return new double[] {-0.5001, -0.5, -0.4999, -0.25, -0.0001, 0.0,
                             0.0001, 0.25, 1.0, 1.2499, 1.25, 1.2501};
    }

    /** Creates the default cumulative probability density test expected values */
    @Override
    public double[] makeCumulativeTestValues() {
        return new double[] {0.0, 0.0, 0.0001, 0.25/1.75, 0.4999/1.75,
                             0.5/1.75, 0.5001/1.75, 0.75/1.75, 1.5/1.75,
                             1.7499/1.75, 1.0, 1.0};
    }

    /** Creates the default probability density test expected values */
    @Override
    public double[] makeDensityTestValues() {
        double d = 1 / 1.75;
        return new double[] {0, d, d, d, d, d, d, d, d, d, d, 0};
    }

    //--- Additional test cases -----------------------------------------------

    /** Test lower bound getter. */
    @Test
    void testGetLowerBound() {
        UniformRealDistribution distribution = makeDistribution();
        assertEquals(-0.5, distribution.getSupportLowerBound(), 0);
    }

    /** Test upper bound getter. */
    @Test
    void testGetUpperBound() {
        UniformRealDistribution distribution = makeDistribution();
        assertEquals(1.25, distribution.getSupportUpperBound(), 0);
    }

    /** Test pre-condition for equal lower/upper bound. */
    @Test
    void testPreconditions1() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new UniformRealDistribution(0, 0);
        });
    }

    /** Test pre-condition for lower bound larger than upper bound. */
    @Test
    void testPreconditions2() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new UniformRealDistribution(1, 0);
        });
    }

    /** Test mean/variance. */
    @Test
    void testMeanVariance() {
        UniformRealDistribution dist;

        dist = new UniformRealDistribution(0, 1);
        assertEquals(0.5, dist.getNumericalMean(), 0);
        assertEquals(dist.getNumericalVariance(), 1/12.0, 0);

        dist = new UniformRealDistribution(-1.5, 0.6);
        assertEquals(dist.getNumericalMean(), -0.45, 0);
        assertEquals(0.3675, dist.getNumericalVariance(), 0);

        dist = new UniformRealDistribution(-0.5, 1.25);
        assertEquals(0.375, dist.getNumericalMean(), 0);
        assertEquals(0.2552083333333333, dist.getNumericalVariance(), 0);
    }

    /**
     * Check accuracy of analytical inverse CDF. Fails if a solver is used
     * with the default accuracy.
     */
    @Test
    void testInverseCumulativeDistribution() {
        UniformRealDistribution dist = new UniformRealDistribution(0, 1e-9);

        assertEquals(2.5e-10, dist.inverseCumulativeProbability(0.25), 0);
    }
}
