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

package org.hipparchus.distribution.discrete;

import org.hipparchus.distribution.IntegerDistribution;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test cases for UniformIntegerDistribution.
 */
public class UniformIntegerDistributionTest extends IntegerDistributionAbstractTest {

    // --- Override tolerance -------------------------------------------------

    @BeforeEach
    @Override
    public void setUp() {
        super.setUp();
        setTolerance(1e-9);
    }

    //--- Implementations for abstract methods --------------------------------

    /** Creates the default discrete distribution instance to use in tests. */
    @Override
    public IntegerDistribution makeDistribution() {
        return new UniformIntegerDistribution(-3, 5);
    }

    /** Creates the default probability density test input values. */
    @Override
    public int[] makeDensityTestPoints() {
        return new int[] {-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6};
    }

    /** Creates the default probability density test expected values. */
    @Override
    public double[] makeDensityTestValues() {
        double d = 1.0 / (5 - -3 + 1);
        return new double[] {0, d, d, d, d, d, d, d, d, d, 0};
    }

    /** Creates the default cumulative probability density test input values. */
    @Override
    public int[] makeCumulativeTestPoints() {
        return makeDensityTestPoints();
    }

    /** Creates the default cumulative probability density test expected values. */
    @Override
    public double[] makeCumulativeTestValues() {
        return new double[] {0, 1 / 9.0, 2 / 9.0, 3 / 9.0, 4 / 9.0, 5 / 9.0,
                             6 / 9.0, 7 / 9.0, 8 / 9.0, 1, 1};
    }

    /** Creates the default inverse cumulative probability test input values */
    @Override
    public double[] makeInverseCumulativeTestPoints() {
        return new double[] {0, 0.001, 0.010, 0.025, 0.050, 0.100, 0.200,
                             0.5, 0.999, 0.990, 0.975, 0.950, 0.900, 1};
    }

    /** Creates the default inverse cumulative probability density test expected values */
    @Override
    public int[] makeInverseCumulativeTestValues() {
        return new int[] {-3, -3, -3, -3, -3, -3, -2, 1, 5, 5, 5, 5, 5, 5};
    }

    //--- Additional test cases -----------------------------------------------

    /** Test mean/variance. */
    @Test
    void testMoments() {
        UniformIntegerDistribution dist;

        dist = new UniformIntegerDistribution(0, 5);
        assertEquals(2.5, dist.getNumericalMean(), 0);
        assertEquals(dist.getNumericalVariance(), 35 / 12.0, 0);

        dist = new UniformIntegerDistribution(0, 1);
        assertEquals(0.5, dist.getNumericalMean(), 0);
        assertEquals(dist.getNumericalVariance(), 3 / 12.0, 0);
    }

    // MATH-1141
    @Test
    void testPreconditionUpperBoundInclusive() {
        try {
            new UniformIntegerDistribution(1, 0);
        } catch (MathIllegalArgumentException e) {
            // Expected.
        }

        // Degenerate case is allowed.
        new UniformIntegerDistribution(0, 0);
    }
}
