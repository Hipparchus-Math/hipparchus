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
package org.hipparchus.distribution.multivariate;

import org.hipparchus.distribution.MultivariateRealDistribution;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.random.RandomGenerator;

/**
 * Base class for multivariate probability distributions.
 */
public abstract class AbstractMultivariateRealDistribution
    implements MultivariateRealDistribution {
    /** RNG instance used to generate samples from the distribution. */
    protected final RandomGenerator random;
    /** The number of dimensions or columns in the multivariate distribution. */
    private final int dimension;

    /** Simple constructor.
     * @param rng Random number generator.
     * @param n Number of dimensions.
     */
    protected AbstractMultivariateRealDistribution(RandomGenerator rng,
                                                   int n) {
        random = rng;
        dimension = n;
    }

    /** {@inheritDoc} */
    @Override
    public void reseedRandomGenerator(long seed) {
        random.setSeed(seed);
    }

    /** {@inheritDoc} */
    @Override
    public int getDimension() {
        return dimension;
    }

    /** {@inheritDoc} */
    @Override
    public abstract double[] sample();

    /** {@inheritDoc} */
    @Override
    public double[][] sample(final int sampleSize) {
        if (sampleSize <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_OF_SAMPLES,
                                                   sampleSize);
        }
        final double[][] out = new double[sampleSize][dimension];
        for (int i = 0; i < sampleSize; i++) {
            out[i] = sample();
        }
        return out;
    }
}
