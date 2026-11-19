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

/**
 * Base class for generating non-uniform random variates from an underlying uniform generator.
 */
public abstract class AbstractNonUniformGenerator implements NonUniformGenerator
{

    /** Underlying random generator. */
    private final RandomGenerator random;

    /**
     * Simple constructor.
     * @param random underlying random generator
     */
    public AbstractNonUniformGenerator(final RandomGenerator random)
    {
        this.random = random;
    }

    /**
     * Generate next uniform variate.
     * @return next uniform variate between 0 and 1
     */
    protected double nextUniform()
    {
        return random.nextDouble();
    }

    /**
     * Generate next normal variate.
     * @return next normal variate witm mean 0 and standard deviation 1
     */
    protected double nextNormal()
    {
        return random.nextGaussian();
    }

}
