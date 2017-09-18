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
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.Pair;
import org.junit.Test;

/**
 * Test class for {@link EnumeratedRealDistribution}.
 */
public class EnumeratedRealDistributionTest {

    /**
     * The distribution object used for testing.
     */
    private final EnumeratedRealDistribution testDistribution;

    /**
     * Creates the default distribution object used for testing.
     */
    public EnumeratedRealDistributionTest() {
        // Non-sorted singleton array with duplicates should be allowed.
        // Values with zero-probability do not extend the support.
        testDistribution = new EnumeratedRealDistribution(
                new double[]{3.0, -1.0, 3.0, 7.0, -2.0, 8.0},
                new double[]{0.2, 0.2, 0.3, 0.3, 0.0, 0.0});
    }

    /**
     * Tests if the {@link EnumeratedRealDistribution} constructor throws
     * exceptions for invalid data.
     */
    @Test
    public void testExceptions() {
        EnumeratedRealDistribution invalid = null;
        try {
            invalid = new EnumeratedRealDistribution(new double[]{1.0, 2.0}, new double[]{0.0});
            fail("Expected MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
        }
        try{
        invalid = new EnumeratedRealDistribution(new double[]{1.0, 2.0}, new double[]{0.0, -1.0});
            fail("Expected MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
        }
        try {
            invalid = new EnumeratedRealDistribution(new double[]{1.0, 2.0}, new double[]{0.0, 0.0});
            fail("Expected MathRuntimeException");
        } catch (MathRuntimeException e) {
        }
        try {
            invalid = new EnumeratedRealDistribution(new double[]{1.0, 2.0}, new double[]{0.0, Double.NaN});
            fail("Expected MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
        }
        try {
            invalid = new EnumeratedRealDistribution(new double[]{1.0, 2.0}, new double[]{0.0, Double.POSITIVE_INFINITY});
            fail("Expected MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
        }
        assertNull("Expected non-initialized DiscreteRealDistribution", invalid);
    }

    /**
     * Tests if the distribution returns proper probability values.
     */
    @Test
    public void testProbability() {
        double[] points = new double[]{-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
        double[] results = new double[]{0, 0.2, 0, 0, 0, 0.5, 0, 0, 0, 0.3, 0};
        for (int p = 0; p < points.length; p++) {
            double density = testDistribution.probability(points[p]);
            assertEquals(results[p], density, 0.0);
        }
    }

    /**
     * Tests if the distribution returns proper density values.
     */
    @Test
    public void testDensity() {
        double[] points = new double[]{-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
        double[] results = new double[]{0, 0.2, 0, 0, 0, 0.5, 0, 0, 0, 0.3, 0};
        for (int p = 0; p < points.length; p++) {
            double density = testDistribution.density(points[p]);
            assertEquals(results[p], density, 0.0);
        }
    }

    /**
     * Tests if the distribution returns proper cumulative probability values.
     */
    @Test
    public void testCumulativeProbability() {
        double[] points = new double[]{-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
        double[] results = new double[]{0, 0.2, 0.2, 0.2, 0.2, 0.7, 0.7, 0.7, 0.7, 1.0, 1.0};
        for (int p = 0; p < points.length; p++) {
            double probability = testDistribution.cumulativeProbability(points[p]);
            assertEquals(results[p], probability, 1e-10);
        }
    }

    /**
     * Tests if the distribution returns proper mean value.
     */
    @Test
    public void testGetNumericalMean() {
        assertEquals(3.4, testDistribution.getNumericalMean(), 1e-10);
    }

    /**
     * Tests if the distribution returns proper variance.
     */
    @Test
    public void testGetNumericalVariance() {
        assertEquals(7.84, testDistribution.getNumericalVariance(), 1e-10);
    }

    /**
     * Tests if the distribution returns proper lower bound.
     */
    @Test
    public void testGetSupportLowerBound() {
        assertEquals(-1, testDistribution.getSupportLowerBound(), 0);
    }

    /**
     * Tests if the distribution returns proper upper bound.
     */
    @Test
    public void testGetSupportUpperBound() {
        assertEquals(7, testDistribution.getSupportUpperBound(), 0);
    }

    /**
     * Tests if the distribution returns properly that the support is connected.
     */
    @Test
    public void testIsSupportConnected() {
        assertTrue(testDistribution.isSupportConnected());
    }


    @Test
    public void testIssue1065() {
        // Test Distribution for inverseCumulativeProbability
        //
        //         ^
        //         |
        // 1.000   +--------------------------------o===============
        //         |                               3|
        //         |                                |
        //         |                             1o=
        // 0.750   +-------------------------> o==  .
        //         |                          3|  . .
        //         |                   0       |  . .
        // 0.5625  +---------------> o==o======   . .
        //         |                 |  .      .  . .
        //         |                 |  .      .  . .
        //         |                5|  .      .  . .
        //         |                 |  .      .  . .
        //         |             o===   .      .  . .
        //         |             |   .  .      .  . .
        //         |            4|   .  .      .  . .
        //         |             |   .  .      .  . .
        // 0.000   +=============----+--+------+--+-+--------------->
        //                      14  18 21     28 31 33
        //
        // sum  = 4+5+0+3+1+3 = 16

        EnumeratedRealDistribution distribution = new EnumeratedRealDistribution(
                new double[] { 14.0, 18.0, 21.0, 28.0, 31.0, 33.0 },
                new double[] { 4.0 / 16.0, 5.0 / 16.0, 0.0 / 16.0, 3.0 / 16.0, 1.0 / 16.0, 3.0 / 16.0 });

        assertEquals(14.0, distribution.inverseCumulativeProbability(0.0000), 0.0);
        assertEquals(14.0, distribution.inverseCumulativeProbability(0.2500), 0.0);
        assertEquals(33.0, distribution.inverseCumulativeProbability(1.0000), 0.0);

        assertEquals(18.0, distribution.inverseCumulativeProbability(0.5000), 0.0);
        assertEquals(18.0, distribution.inverseCumulativeProbability(0.5624), 0.0);
        assertEquals(28.0, distribution.inverseCumulativeProbability(0.5626), 0.0);
        assertEquals(31.0, distribution.inverseCumulativeProbability(0.7600), 0.0);
        assertEquals(18.0, distribution.inverseCumulativeProbability(0.5625), 0.0);
        assertEquals(28.0, distribution.inverseCumulativeProbability(0.7500), 0.0);
    }

    @Test
    public void testCreateFromDoubles() {
        final double[] data = new double[] {0, 1, 1, 2, 2, 2};
        EnumeratedRealDistribution distribution = new EnumeratedRealDistribution(data);
        assertEquals(0.5, distribution.probability(2), 0);
        assertEquals(0.5, distribution.cumulativeProbability(1), 0);
    }

    @Test
    public void testGetPmf() {
        final double[] data = new double[] {0,0,1,1,2,2,2,2,3,4};
        final EnumeratedRealDistribution distribution = new EnumeratedRealDistribution(data);
        final List<Pair<Double, Double>> pmf = distribution.getPmf();
        assertEquals(5, pmf.size());
        final Map<Double, Double> pmfMap = new HashMap<Double, Double>();
        pmfMap.put(0d, 0.2);
        pmfMap.put(1d, 0.2);
        pmfMap.put(2d, 0.4);
        pmfMap.put(3d, 0.1);
        pmfMap.put(4d, 0.1);
        for (int i = 0; i < 5; i++) {
            assertEquals(pmf.get(i).getSecond(), pmfMap.get(pmf.get(i).getFirst()), 0);
        }
    }
}
