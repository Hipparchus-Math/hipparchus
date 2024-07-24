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

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for FDistribution.
 */
public class FDistributionTest extends RealDistributionAbstractTest {

    //-------------- Implementations for abstract methods -----------------------

    /** Creates the default continuous distribution instance to use in tests. */
    @Override
    public FDistribution makeDistribution() {
        return new FDistribution(5.0, 6.0);
    }

    /** Creates the default cumulative probability distribution test input values */
    @Override
    public double[] makeCumulativeTestPoints() {
        // quantiles computed using R version 2.9.2
        return new double[] {0.0346808448626, 0.0937009113303, 0.143313661184, 0.202008445998, 0.293728320107,
                20.8026639595, 8.74589525602, 5.98756512605, 4.38737418741, 3.10751166664};
    }

    /** Creates the default cumulative probability density test expected values */
    @Override
    public double[] makeCumulativeTestValues() {
        return new double[] {0.001, 0.01, 0.025, 0.05, 0.1, 0.999, 0.990, 0.975, 0.950, 0.900};
    }

    /** Creates the default probability density test expected values */
    @Override
    public double[] makeDensityTestValues() {
        return new double[] {0.0689156576706, 0.236735653193, 0.364074131941, 0.481570789649, 0.595880479994,
                0.000133443915657, 0.00286681303403, 0.00969192007502, 0.0242883861471, 0.0605491314658};
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
    void testDfAccessors() {
        FDistribution dist = (FDistribution) getDistribution();
        assertEquals(5d, dist.getNumeratorDegreesOfFreedom(), Double.MIN_VALUE);
        assertEquals(6d, dist.getDenominatorDegreesOfFreedom(), Double.MIN_VALUE);
    }

    @Test
    void testPreconditions() {
        try {
            new FDistribution(0, 1);
            fail("Expecting MathIllegalArgumentException for df = 0");
        } catch (MathIllegalArgumentException ex) {
            // Expected.
        }
        try {
            new FDistribution(1, 0);
            fail("Expecting MathIllegalArgumentException for df = 0");
        } catch (MathIllegalArgumentException ex) {
            // Expected.
        }
    }

    @Test
    void testLargeDegreesOfFreedom() {
        FDistribution fd = new FDistribution(100000, 100000);
        double p = fd.cumulativeProbability(.999);
        double x = fd.inverseCumulativeProbability(p);
        assertEquals(.999, x, 1.0e-5);
    }

    @Test
    void testSmallDegreesOfFreedom() {
        FDistribution fd = new FDistribution(1, 1);
        double p = fd.cumulativeProbability(0.975);
        double x = fd.inverseCumulativeProbability(p);
        assertEquals(0.975, x, 1.0e-5);

        fd = new FDistribution(1, 2);
        p = fd.cumulativeProbability(0.975);
        x = fd.inverseCumulativeProbability(p);
        assertEquals(0.975, x, 1.0e-5);
    }

    @Test
    void testMoments() {
        final double tol = 1e-9;
        FDistribution dist;

        dist = new FDistribution(1, 2);
        assertTrue(Double.isNaN(dist.getNumericalMean()));
        assertTrue(Double.isNaN(dist.getNumericalVariance()));

        dist = new FDistribution(1, 3);
        assertEquals(dist.getNumericalMean(), 3d / (3d - 2d), tol);
        assertTrue(Double.isNaN(dist.getNumericalVariance()));

        dist = new FDistribution(1, 5);
        assertEquals(dist.getNumericalMean(), 5d / (5d - 2d), tol);
        assertEquals(dist.getNumericalVariance(), (2d * 5d * 5d * 4d) / 9d, tol);
    }

    @Test
    void testMath785() {
        // this test was failing due to inaccurate results from ContinuedFraction.

        assertDoesNotThrow(() -> {
            double prob = 0.01;
            FDistribution f = new FDistribution(200000, 200000);
            double result = f.inverseCumulativeProbability(prob);
            assertTrue(result < 1.0);
        }, "Failing to calculate inverse cumulative probability");
    }
}
