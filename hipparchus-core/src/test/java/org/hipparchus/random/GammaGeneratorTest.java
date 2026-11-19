/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.random;

import org.hipparchus.distribution.RealDistribution;
import org.hipparchus.distribution.continuous.GammaDistribution;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.SplittableRandom;

public class GammaGeneratorTest extends AbstractNonUniformGeneratorTest
{

    /** {@inheritDoc} */
    @Override
    protected NonUniformGenerator buildGenerator(final RandomGenerator uniform,
                                                 final double... parameters)
    {
        return new GammaGenerator(uniform, parameters[0], parameters[1]);
    }

    /** {@inheritDoc} */
    @Override
    protected RealDistribution buildDistribution(double... parameters)
    {
        return new GammaDistribution(parameters[0], parameters[1]);
    }

    @Test
    public void testDistributionGreaterThan1()
    {
        Assertions.assertEquals(0.0, ksStatistic(0x0678dad954a74560L, 3.0, 2.0), 2.1e-3);
    }

    @Test
    public void testDistributionSmallerThan1()
    {
        Assertions.assertEquals(0.0, ksStatistic(0x35a7b28cf4ddc7c7L, 0.25, 3.0), 2.0e-3);
    }

    @Test
    public void testMean()
    {
        final SplittableRandom seedsGenerator = new SplittableRandom(0x79e52d38915419bbL);
        for (double alpha = 0.75; alpha < 5.0; alpha += 0.125)
        {
            for (double theta = 1.0; theta < 4.0; theta += 0.5)
            {
                final double mean = mean(seedsGenerator.nextLong(), alpha, theta);
                final double theoretical = alpha * theta;
                Assertions.assertEquals(theoretical, mean, 0.0078 * theoretical);
            }
        }
    }

    @Test
    public void testVariance()
    {
        final SplittableRandom seedsGenerator = new SplittableRandom(0x14fec21eec8186ccL);
        for (double alpha = 0.75; alpha < 5.0; alpha += 0.125)
        {
            for (double theta = 1.0; theta < 4.0; theta += 0.5)
            {
                final double variance = variance(seedsGenerator.nextLong(), alpha, theta);
                final double theoretical = alpha * theta * theta;
                Assertions.assertEquals(theoretical, variance, 0.021 * theoretical);
            }
        }
    }

    @Test
    public void testMaxRejection()
    {
        try {
            new GammaGenerator(new RandomAdaptorTest.ConstantGenerator(1.23), 2.3, 4.5).nextVariate();
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assertions.assertEquals(LocalizedCoreFormats.INTERNAL_ERROR, mise.getSpecifier());
        }
    }

}
