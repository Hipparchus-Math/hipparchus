/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
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

import static org.junit.Assert.assertEquals;

import org.hipparchus.distribution.RealDistribution;
import org.junit.Test;

/**
 * Test cases for ConstantRealDistribution.
 */
public class ConstantRealDistributionTest extends RealDistributionAbstractTest {

    // --- Override tolerance -------------------------------------------------

    @Override
    public void setUp() {
        super.setUp();
        setTolerance(0);
    }

    //--- Implementations for abstract methods --------------------------------

    /** Creates the default uniform real distribution instance to use in tests. */
    @Override
    public ConstantRealDistribution makeDistribution() {
        return new ConstantRealDistribution(1);
    }

    /** Creates the default cumulative probability distribution test input values */
    @Override
    public double[] makeCumulativeTestPoints() {
        return new double[] { 0, 0.5, 1 };
    }

    /** Creates the default cumulative probability distribution test expected values */
    @Override
    public double[] makeCumulativeTestValues() {
        return new double[] { 0, 0, 1 };
    }

    /** Creates the default probability density test expected values */
    @Override
    public double[] makeDensityTestValues() {
        return new double[] { 0, 0, 1 };
    }

    /** Override default test, verifying that inverse cum is constant */
    @Override
    @Test
    public void testInverseCumulativeProbabilities() {
        RealDistribution dist = getDistribution();
        for (double x : getCumulativeTestValues()) {
            assertEquals(1,dist.inverseCumulativeProbability(x), 0);
        }
    }

    //--- Additional test cases -----------------------------------------------

    @Test
    public void testMeanVariance() {
        ConstantRealDistribution dist;

        dist = new ConstantRealDistribution(-1);
        assertEquals(dist.getNumericalMean(), -1, 0d);
        assertEquals(dist.getNumericalVariance(), 0, 0d);
    }

}
