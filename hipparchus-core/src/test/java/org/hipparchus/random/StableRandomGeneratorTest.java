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
package org.hipparchus.random;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.Assert;
import org.junit.Test;

/**
 * The class <code>StableRandomGeneratorTest</code> contains tests for the class
 * {@link StableRandomGenerator}
 *
 */
public class StableRandomGeneratorTest {

    private RandomGenerator rg = new Well19937c(100);
    private final static int sampleSize = 10000;

    /**
     * Run the double nextDouble() method test Due to leptokurtic property the
     * acceptance range is widened.
     *
     * TODO: verify that tolerance this wide is really OK
     */
    @Test
    public void testNextDouble() {
        StableRandomGenerator generator = new StableRandomGenerator(rg, 1.3,
                0.1);
        double[] sample = new double[2 * sampleSize];
        for (int i = 0; i < sample.length; ++i) {
            sample[i] = generator.nextNormalizedDouble();
        }
        Assert.assertEquals(0.0, UnitTestUtils.mean(sample), 0.3);
    }

    /**
     * If alpha = 2, than it must be Gaussian distribution
     */
    @Test
    public void testGaussianCase() {
        StableRandomGenerator generator = new StableRandomGenerator(rg, 2d, 0.0);

        double[] sample = new double[sampleSize];
        for (int i = 0; i < sample.length; ++i) {
            sample[i] = generator.nextNormalizedDouble();
        }
        Assert.assertEquals(0.0, UnitTestUtils.mean(sample), 0.02);
        Assert.assertEquals(1.0, UnitTestUtils.variance(sample), 0.02);
    }

    /**
     * If alpha = 1, than it must be Cauchy distribution
     */
    @Test
    public void testCauchyCase() {
        StableRandomGenerator generator = new StableRandomGenerator(rg, 1d, 0.0);

        final double[] values = new double[sampleSize];
        for (int i = 0; i < sampleSize; ++i) {
            values[i] = generator.nextNormalizedDouble();
        }

        // Standard Cauchy distribution should have zero median and mode
        double median = UnitTestUtils.median(values);
        Assert.assertEquals(0.0, median, 0.2);
    }

    /**
     * Input parameter range tests
     */
    @Test
    public void testAlphaRangeBelowZero() {
        try {
            new StableRandomGenerator(rg, -1.0, 0.0);
            Assert.fail("Expected MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            Assert.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_LEFT, e.getSpecifier());
            Assert.assertEquals(-1.0, ((Double) e.getParts()[0]).doubleValue(), 1.0e-10);
        }
    }

    @Test
    public void testAlphaRangeAboveTwo() {
        try {
            new StableRandomGenerator(rg, 3.0, 0.0);
            Assert.fail("Expected MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            Assert.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_LEFT, e.getSpecifier());
            Assert.assertEquals(3.0, ((Double) e.getParts()[0]).doubleValue(), 1.0e-10);
        }
    }

    @Test
    public void testBetaRangeBelowMinusOne() {
        try {
            new StableRandomGenerator(rg, 1.0, -2.0);
            Assert.fail("Expected MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            Assert.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, e.getSpecifier());
            Assert.assertEquals(-2.0, ((Double) e.getParts()[0]).doubleValue(), 1.0e-10);
        }
    }

    @Test
    public void testBetaRangeAboveOne() {
        try {
            new StableRandomGenerator(rg, 1.0, 2.0);
            Assert.fail("Expected MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            Assert.assertEquals(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, e.getSpecifier());
            Assert.assertEquals(2.0, ((Double) e.getParts()[0]).doubleValue(), 1.0e-10);
        }
    }
}
