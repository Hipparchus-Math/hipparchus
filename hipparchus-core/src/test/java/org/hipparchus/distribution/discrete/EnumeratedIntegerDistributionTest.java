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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.Pair;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test class for {@link EnumeratedIntegerDistribution}.
 */
public class EnumeratedIntegerDistributionTest {

    /**
     * The distribution object used for testing.
     */
    private final EnumeratedIntegerDistribution testDistribution;

    /**
     * Creates the default distribution object used for testing.
     */
    public EnumeratedIntegerDistributionTest() {
        // Non-sorted singleton array with duplicates should be allowed.
        // Values with zero-probability do not extend the support.
        testDistribution = new EnumeratedIntegerDistribution(
                new int[]{3, -1, 3, 7, -2, 8},
                new double[]{0.2, 0.2, 0.3, 0.3, 0.0, 0.0});
    }

    /**
     * Tests if the EnumeratedIntegerDistribution constructor throws
     * exceptions for invalid data.
     */
    @Test
    public void testExceptions() {
        EnumeratedIntegerDistribution invalid = null;
        try {
            new EnumeratedIntegerDistribution(new int[]{1, 2}, new double[]{0.0});
            Assertions.fail("Expected MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
        }
        try {
            new EnumeratedIntegerDistribution(new int[]{1, 2}, new double[]{0.0, -1.0});
            Assertions.fail("Expected MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
        }
        try {
            new EnumeratedIntegerDistribution(new int[]{1, 2}, new double[]{0.0, 0.0});
            Assertions.fail("Expected MathRuntimeException");
        } catch (MathRuntimeException e) {
        }
        try {
          new EnumeratedIntegerDistribution(new int[]{1, 2}, new double[]{0.0, Double.NaN});
            Assertions.fail("Expected MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
        }
        try {
        new EnumeratedIntegerDistribution(new int[]{1, 2}, new double[]{0.0, Double.POSITIVE_INFINITY});
            Assertions.fail("Expected NotFiniteNumberException");
        } catch (MathIllegalArgumentException e) {
        }
        Assertions.assertNull(invalid, "Expected non-initialized DiscreteRealDistribution");
    }

    /**
     * Tests if the distribution returns proper probability values.
     */
    @Test
    public void testProbability() {
        int[] points = new int[]{-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8};
        double[] results = new double[]{0, 0.2, 0, 0, 0, 0.5, 0, 0, 0, 0.3, 0};
        for (int p = 0; p < points.length; p++) {
            double probability = testDistribution.probability(points[p]);
            Assertions.assertEquals(results[p], probability, 0.0);
        }
    }

    /**
     * Tests if the distribution returns proper cumulative probability values.
     */
    @Test
    public void testCumulativeProbability() {
        int[] points = new int[]{-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8};
        double[] results = new double[]{0, 0.2, 0.2, 0.2, 0.2, 0.7, 0.7, 0.7, 0.7, 1.0, 1.0};
        for (int p = 0; p < points.length; p++) {
            double probability = testDistribution.cumulativeProbability(points[p]);
            Assertions.assertEquals(results[p], probability, 1e-10);
        }
    }

    /**
     * Tests if the distribution returns proper mean value.
     */
    @Test
    public void testGetNumericalMean() {
        Assertions.assertEquals(3.4, testDistribution.getNumericalMean(), 1e-10);
    }

    /**
     * Tests if the distribution returns proper variance.
     */
    @Test
    public void testGetNumericalVariance() {
        Assertions.assertEquals(7.84, testDistribution.getNumericalVariance(), 1e-10);
    }

    /**
     * Tests if the distribution returns proper lower bound.
     */
    @Test
    public void testGetSupportLowerBound() {
        Assertions.assertEquals(-1, testDistribution.getSupportLowerBound());
    }

    /**
     * Tests if the distribution returns proper upper bound.
     */
    @Test
    public void testGetSupportUpperBound() {
        Assertions.assertEquals(7, testDistribution.getSupportUpperBound());
    }

    /**
     * Tests if the distribution returns properly that the support is connected.
     */
    @Test
    public void testIsSupportConnected() {
        Assertions.assertTrue(testDistribution.isSupportConnected());
    }

    @Test
    public void testCreateFromIntegers() {
        final int[] data = new int[] {0, 1, 1, 2, 2, 2};
        EnumeratedIntegerDistribution distribution = new EnumeratedIntegerDistribution(data);
        assertEquals(0.5, distribution.probability(2), 0);
        assertEquals(0.5, distribution.cumulativeProbability(1), 0);
    }

    @Test
    public void testGetPmf() {
        final int[] values = new int[] {0,1,2,3,4};
        final double[] masses = new double[] {0.2, 0.2, 0.4, 0.1, 0.1};
        final EnumeratedIntegerDistribution distribution = new EnumeratedIntegerDistribution(values, masses);
        final List<Pair<Integer, Double>> pmf = distribution.getPmf();
        assertEquals(5, pmf.size());
        final Map<Integer, Double> pmfMap = new HashMap<Integer, Double>();
        for (int i = 0; i < 5; i++) {
            pmfMap.put(i, masses[i]);
        }
        for (int i = 0; i < 5; i++) {
            assertEquals(pmf.get(i).getSecond(), pmfMap.get(pmf.get(i).getFirst()), 0);
        }
    }
}
