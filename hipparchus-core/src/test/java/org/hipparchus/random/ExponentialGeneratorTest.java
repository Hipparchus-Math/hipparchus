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
import org.hipparchus.distribution.continuous.ExponentialDistribution;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.SplittableRandom;

public class ExponentialGeneratorTest extends AbstractNonUniformGeneratorTest
{

    /** {@inheritDoc} */
    @Override
    protected NonUniformGenerator buildGenerator(final RandomGenerator uniform,
                                                 final double... parameters)
    {
        return new ExponentialGenerator(uniform, parameters[0]);
    }

    /** {@inheritDoc} */
    @Override
    protected RealDistribution buildDistribution(double... parameters)
    {
        return new ExponentialDistribution(1.0 / parameters[0]);
    }

    @Test
    public void testDistribution()
    {
        Assertions.assertEquals(0.0, ksStatistic(0x60c0c3a455e3a8e8L, 5.0), 8.0e-3);
    }

    @Test
    public void testMean()
    {
        final SplittableRandom seedsGenerator = new SplittableRandom(0x44365ff34bd58c36L);
        for (double rate = 5.0; rate < 20.0; rate += 0.125)
        {
            final double mean = mean(seedsGenerator.nextLong(), rate);
            Assertions.assertEquals(1.0 / rate, mean, 1.2e-3);
        }
    }

    @Test
    public void testVariance()
    {
        final SplittableRandom seedsGenerator = new SplittableRandom(0x4e10509f9526c087L);
        for (double rate = 5.0; rate < 20.0; rate += 0.125)
        {
            final double variance = variance(seedsGenerator.nextLong(), rate);
            Assertions.assertEquals(1.0 / (rate * rate), variance, 6.0e-4);
        }
    }

}
