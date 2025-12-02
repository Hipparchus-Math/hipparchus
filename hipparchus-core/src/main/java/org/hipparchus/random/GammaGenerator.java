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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;

/**
 * Procedure to generate variates for gamma distribution.
 * <p>
 * The method used for generation is based on Marsaglia and Tsang
 * 2000 paper <a href="https://doi.org/10.1145/358407.358414">A
 * Simple Method for Generating Gamma Variables</a>
 * </p>
 */
public class GammaGenerator extends AbstractNonUniformGenerator
{

    /** Maximum number of rejections allowed before failing. */
    private static final int MAX_REJECTION = 100;

    /** Factor used in squeeze method. */
    private static final double SQUEEZE_FACTOR = 0.0331;

    /** Indicator for shape factor smaller than 1. */
    private final boolean smallShape;

    /** D factor from Marsaglia and Tsang paper. */
    private final double d;

    /** C factor from Marsaglia and Tsang paper. */
    private final double c;

    /** Shape parameter */
    private final double alpha;

    /** Scale parameter */
    private final double theta;

    /**
     * Simple constructor.
     * @param random underlying random generator
     * @param alpha shape parameter
     * @param theta scale parameter (this is the inverse of the rate parameter λ)
     */
    public GammaGenerator(final RandomGenerator random, final double alpha, final double theta)
    {
        super(random);

        // fix small shape factor to ensure we use at least 1
        final double fixedShape;
        if (alpha < 1.0)
        {
            fixedShape = alpha + 1.0;
            smallShape = true;
        }
        else
        {
            fixedShape = alpha;
            smallShape = false;
        }
        this.d     = (3.0 * fixedShape - 1.0) / 3.0;
        this.c     = 1.0 / FastMath.sqrt(9.0 * fixedShape - 3.0);

        this.alpha = alpha;
        this.theta = theta;

    }

    /** {@inheritDoc} */
    @Override
    public double nextVariate()
    {
        final double canonical = canonicalVariate();
        return theta * (smallShape ? canonical * FastMath.pow(nextUniform(), 1.0 / alpha): canonical);
    }

    /**
     * Generate canonical variate for shape at least 1 and unit scale
     * @return canonical variate
     */
    private double canonicalVariate()
    {
        for (int i = 0; i < MAX_REJECTION; ++i)
        {

            // generate the normal variate
            double x = 0;
            double n = -1;
            for (int j = 0; n <= 0 && j < MAX_REJECTION; ++j)
            {
                x = nextNormal();
                n = 1 + c * x;
            }

            // intermediate variate v
            final double v  = n * n * n;
            final double dv = d * v;

            // generate uniform variate
            final double u = nextUniform();
            final double x2 = x * x;

            // the first alternative (with SQUEEZE_FACTOR) is evaluated first;
            // as Java short-circuits the || boolean operator, we often can
            // avoid the costly computation of logarithms in the second alternative
            if (u < 1 - SQUEEZE_FACTOR * x2 * x2 ||
                FastMath.log(u) < 0.5 * x * x + d - dv + d * FastMath.log(v))
            {
                return dv;
            }

            // reject this attempt, we need to try again with new random variables

        }

        // this should never happen
        throw new MathIllegalStateException(LocalizedCoreFormats.INTERNAL_ERROR);

    }

}
