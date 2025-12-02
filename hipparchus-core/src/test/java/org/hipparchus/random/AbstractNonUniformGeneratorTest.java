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

import java.util.SplittableRandom;

public abstract class AbstractNonUniformGeneratorTest
{

    /** Sample size for generators. */
    private static final int SAMPLE_SIZE = 100000;

    /**
     * Build a generator for non-uniform variates.
     * @param uniform underlying generator for uniform variates
     * @param parameters non-uniform distribution parameters
     * @return generator for non-uniform variates
     */
    protected abstract NonUniformGenerator buildGenerator(RandomGenerator uniform, double... parameters);

    /**
     * Build non-uniform distribution.
     * @param parameters non-uniform distribution parameters
     * @return non-uniform distribution
     */
    protected abstract RealDistribution buildDistribution(double... parameters);

    protected double ksStatistic(final long seed, double... parameters)
    {

        // build generator and distribution
        final SplittableRandom    seedsGenerator = new SplittableRandom(seed);
        final RandomGenerator     uniform        = new Well19937a(seedsGenerator.nextLong());
        final NonUniformGenerator generator      = buildGenerator(uniform, parameters);
        final RealDistribution    distribution   = buildDistribution(parameters);

        final KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest();

        // generate sample
        final double[] data = new double[SAMPLE_SIZE];
        for (int i = 0; i < SAMPLE_SIZE; ++i)
        {
            data[i] = generator.nextVariate();
        }

        return ksTest.kolmogorovSmirnovStatistic(distribution, data);

    }

    protected double mean(final long seed, double... parameters)
    {
        final NonUniformGenerator generator = buildGenerator(new Well19937a(seed), parameters);
        double sum = 0.0;
        for (int i = 0; i < SAMPLE_SIZE; ++i)
        {
            sum += generator.nextVariate();
        }
        return sum / SAMPLE_SIZE;
    }

    protected double variance(final long seed, double... parameters)
    {
        // compute variance using one-pass Welford method
        final NonUniformGenerator generator = buildGenerator(new Well19937a(seed), parameters);
        double m = 0;
        double s = 0;
        for (int i = 0; i < SAMPLE_SIZE; ++i)
        {
            final double x = generator.nextVariate();
            final double oldM = m;
            m += (x - m) / (i + 1);
            s += (x - m) * (x - oldM);
        }
        return s / (SAMPLE_SIZE - 1);
    }

}
