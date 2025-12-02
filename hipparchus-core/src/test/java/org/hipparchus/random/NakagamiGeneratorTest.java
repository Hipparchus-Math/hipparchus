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
import org.hipparchus.distribution.continuous.NakagamiDistribution;
import org.hipparchus.special.Gamma;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.SplittableRandom;

public class NakagamiGeneratorTest extends AbstractNonUniformGeneratorTest
{

    /** {@inheritDoc} */
    @Override
    protected NonUniformGenerator buildGenerator(final RandomGenerator uniform, final double... parameters)
    {
        return new NakagamiGenerator(uniform, parameters[0], parameters[1]);
    }

    /** {@inheritDoc} */
    @Override
    protected RealDistribution buildDistribution(double... parameters)
    {
        return new NakagamiDistribution(parameters[0], parameters[1]);
    }

    @Test
    public void testDistribution()
    {
        Assertions.assertEquals(0.0, ksStatistic(0x099a9800e248b9f4L, 3.0, 2.0), 3.3e-3);
    }

    @Test
    public void testMean()
    {
        final SplittableRandom seedsGenerator = new SplittableRandom(0x7b2173aba272b5e7L);
        for (double m = 0.75; m < 5.0; m += 0.125)
        {
            for (double omega = 1.0; omega < 4.0; omega += 0.5)
            {
                final double mean = mean(seedsGenerator.nextLong(), m, omega);
                final double r = Gamma.gamma(m + 0.5) / Gamma.gamma(m);
                final double theoretical = r * FastMath.sqrt(omega / m);
                Assertions.assertEquals(theoretical, mean, 0.004 * theoretical);
            }
        }
    }

    @Test
    public void testVariance()
    {
        final SplittableRandom seedsGenerator = new SplittableRandom(0x748ffa6ce93daaebL);
        for (double m = 0.75; m < 5.0; m += 0.125)
        {
            for (double omega = 1.0; omega < 4.0; omega += 0.5)
            {
                final double variance = variance(seedsGenerator.nextLong(), m, omega);
                final double r = Gamma.gamma(m + 0.5) / Gamma.gamma(m);
                final double theoretical = omega * (1 - r * r / m);
                Assertions.assertEquals(theoretical, variance, 0.014 * theoretical);
            }
        }
    }

}
