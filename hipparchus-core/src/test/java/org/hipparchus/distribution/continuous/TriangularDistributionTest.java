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
 * Test cases for {@link TriangularDistribution}.
 */
public class TriangularDistributionTest extends RealDistributionAbstractTest {

    // --- Override tolerance -------------------------------------------------

    @BeforeEach
    @Override
    public void setUp() {
        super.setUp();
        setTolerance(1e-4);
    }

    //--- Implementations for abstract methods --------------------------------

    /**
     * Creates the default triangular distribution instance to use in tests.
     */
    @Override
    public TriangularDistribution makeDistribution() {
        // Left side 5 wide, right side 10 wide.
        return new TriangularDistribution(-3, 2, 12);
    }

    /**
     * Creates the default cumulative probability distribution test input
     * values.
     */
    @Override
    public double[] makeCumulativeTestPoints() {
        return new double[] { -3.0001,                 // below lower limit
                              -3.0,                    // at lower limit
                              -2.0, -1.0, 0.0, 1.0,    // on lower side
                              2.0,                     // at mode
                              3.0, 4.0, 10.0, 11.0,    // on upper side
                              12.0,                    // at upper limit
                              12.0001                  // above upper limit
                            };
    }

    /**
     * Creates the default cumulative probability density test expected values.
     */
    @Override
    public double[] makeCumulativeTestValues() {
        // Top at 2 / (b - a) = 2 / (12 - -3) = 2 / 15 = 7.5
        // Area left  = 7.5 * 5  * 0.5 = 18.75 (1/3 of the total area)
        // Area right = 7.5 * 10 * 0.5 = 37.5  (2/3 of the total area)
        // Area total = 18.75 + 37.5 = 56.25
        // Derivative left side = 7.5 / 5 = 1.5
        // Derivative right side = -7.5 / 10 = -0.75
        double third = 1 / 3.0;
        double left = 18.75;
        double area = 56.25;
        return new double[] { 0.0,
                              0.0,
                              0.75 / area, 3 / area, 6.75 / area, 12 / area,
                              third,
                              (left + 7.125) / area, (left + 13.5) / area,
                              (left + 36) / area, (left + 37.125) / area,
                              1.0,
                              1.0
                            };
    }

    /**
     * Creates the default inverse cumulative probability distribution test
     * input values.
     */
    @Override
    public double[] makeInverseCumulativeTestPoints() {
        // Exclude the points outside the limits, as they have cumulative
        // probability of zero and one, meaning the inverse returns the
        // limits and not the points outside the limits.
        double[] points = makeCumulativeTestValues();
        double[] points2 = new double[points.length-2];
        System.arraycopy(points, 1, points2, 0, points2.length);
        return points2;
        //return Arrays.copyOfRange(points, 1, points.length - 1);
    }

    /**
     * Creates the default inverse cumulative probability density test expected
     * values.
     */
    @Override
    public double[] makeInverseCumulativeTestValues() {
        // Exclude the points outside the limits, as they have cumulative
        // probability of zero and one, meaning the inverse returns the
        // limits and not the points outside the limits.
        double[] points = makeCumulativeTestPoints();
        double[] points2 = new double[points.length-2];
        System.arraycopy(points, 1, points2, 0, points2.length);
        return points2;
        //return Arrays.copyOfRange(points, 1, points.length - 1);
    }

    /** Creates the default probability density test expected values. */
    @Override
    public double[] makeDensityTestValues() {
        return new double[] { 0,
                              0,
                              2 / 75.0, 4 / 75.0, 6 / 75.0, 8 / 75.0,
                              10 / 75.0,
                              9 / 75.0, 8 / 75.0, 2 / 75.0, 1 / 75.0,
                              0,
                              0
                            };
    }

    //--- Additional test cases -----------------------------------------------

    /** Test lower bound getter. */
    @Test
    void testGetLowerBound() {
        TriangularDistribution distribution = makeDistribution();
        assertEquals(-3.0, distribution.getSupportLowerBound(), 0);
    }

    /** Test upper bound getter. */
    @Test
    void testGetUpperBound() {
        TriangularDistribution distribution = makeDistribution();
        assertEquals(12.0, distribution.getSupportUpperBound(), 0);
    }

    /** Test pre-condition for equal lower/upper limit. */
    @Test
    void testPreconditions1() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new TriangularDistribution(0, 0, 0);
        });
    }

    /** Test pre-condition for lower limit larger than upper limit. */
    @Test
    void testPreconditions2() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new TriangularDistribution(1, 1, 0);
        });
    }

    /** Test pre-condition for mode larger than upper limit. */
    @Test
    void testPreconditions3() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new TriangularDistribution(0, 2, 1);
        });
    }

    /** Test pre-condition for mode smaller than lower limit. */
    @Test
    void testPreconditions4() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new TriangularDistribution(2, 1, 3);
        });
    }

    /** Test mean/variance. */
    @Test
    void testMeanVariance() {
        TriangularDistribution dist;

        dist = new TriangularDistribution(0, 0.5, 1.0);
        assertEquals(0.5, dist.getNumericalMean(), 0);
        assertEquals(dist.getNumericalVariance(), 1 / 24.0, 0);

        dist = new TriangularDistribution(0, 1, 1);
        assertEquals(dist.getNumericalMean(), 2 / 3.0, 0);
        assertEquals(dist.getNumericalVariance(), 1 / 18.0, 0);

        dist = new TriangularDistribution(-3, 2, 12);
        assertEquals(dist.getNumericalMean(), 3 + (2 / 3.0), 0);
        assertEquals(dist.getNumericalVariance(), 175 / 18.0, 0);
    }
}
