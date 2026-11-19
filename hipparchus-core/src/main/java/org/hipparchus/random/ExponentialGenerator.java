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

import org.hipparchus.util.FastMath;

/**
 * Procedure to generate variates for exponential distribution.
 */
public class ExponentialGenerator extends AbstractNonUniformGenerator
{

    /** Rate parameter */
    private final double lambda;

    /**
     * Simple constructor.
     * @param random underlying random generator
     * @param lambda rate parameter
     */
    public ExponentialGenerator(final RandomGenerator random, final double lambda)
    {
        super(random);
        this.lambda = lambda;
    }

    /** {@inheritDoc} */
    @Override
    public double nextVariate()
    {
        return -FastMath.log(nextUniform()) / lambda;
    }

}
