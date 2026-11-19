/*
 * Licensed to the Hipparchus project under one or more
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

package org.hipparchus.distribution.continuous;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.special.Erf;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/**
 * Implementation of the half normal (gaussian) distribution.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Half-normal_distribution">Half-Normal distribution (Wikipedia)</a>
 * @see <a href="https://mathworld.wolfram.com/Half-NormalDistribution.html">Half-Normal distribution (MathWorld)</a>
 */
public class HalfNormalDistribution extends AbstractRealDistribution {
    /** Serializable version identifier. */
    private static final long serialVersionUID = 20250915L;
    /** &radic;(2) */
    private static final double SQRT2 = FastMath.sqrt(2.0);
    /** Parameter of this distribution (known as scale parameter). */
    private final double scale;
    /** Mean of this distribution. */
    private final double mean;
    /** Standard deviation of this distribution. */
    private final double standardDeviation;
    /** The value of {@code 0.5*log(2/pi)} stored for faster computation. */
    private final double halfLog2OverPi;

    /**
     * Create a half normal distribution using the scale parameter of the distribution.
     *
     * @param scale Scale parameter of the distribution.
     * @throws MathIllegalArgumentException if {@code scale <= 0}.
     */
    public HalfNormalDistribution(double scale)
        throws MathIllegalArgumentException {
        if (scale <= 0.) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.SCALE, scale);
        }
        this.scale = scale;
        this.mean  = scale*SQRT2/FastMath.sqrt(FastMath.PI);
        this.standardDeviation = scale*FastMath.sqrt(1.0-2.0/FastMath.PI);
        this.halfLog2OverPi = 0.5 * FastMath.log(2.0/FastMath.PI);
    }

    /**
     * Access the mean.
     *
     * @return the mean for this distribution.
     */
    public double getMean() {
        return mean;
    }

    /**
     * Access the standard deviation.
     *
     * @return the standard deviation for this distribution.
     */
    public double getStandardDeviation() {
        return standardDeviation;
    }

    /** {@inheritDoc} */
    @Override
    public double density(double x) {
        if (x < 0) {
            return 0;
        }
        return FastMath.exp(logDensity(x));
    }

    /** {@inheritDoc} */
    @Override
    public double logDensity(double x) {
        if (x < 0) {
            return Double.NEGATIVE_INFINITY;
        }
        final double x0 = x / scale;
        return -0.5 * x0 * x0 - FastMath.log(scale) + halfLog2OverPi;
    }

    /**
     * {@inheritDoc}
     *
     */
    @Override
    public double cumulativeProbability(double x)  {
        if (x < 0) {
            return 0;
        }
        return Erf.erf(x / (scale * SQRT2));
    }

    /** {@inheritDoc} */
    @Override
    public double inverseCumulativeProbability(final double p) throws MathIllegalArgumentException {
        MathUtils.checkRangeInclusive(p, 0, 1);
        return scale * SQRT2 * Erf.erfInv(p);
    }

    /** {@inheritDoc} */
    @Override
    public double probability(double x0,
                              double x1)
        throws MathIllegalArgumentException {

        if (x0 < 0 || x1 < 0) {
            return 0;
        }
        if (x0 > x1) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT,
                                                   x0, x1, true);
        }
        final double denom = scale * SQRT2;
        final double v0 = x0 / denom;
        final double v1 = x1 / denom;
        return Erf.erf(v0, v1);
    }

    /**
     * {@inheritDoc}
     *
     * Get the mean of the distribution.
     */
    @Override
    public double getNumericalMean() {
        return getMean();
    }

    /**
     * {@inheritDoc}
     *
     * Get the variance of the distribution.
     */
    @Override
    public double getNumericalVariance() {
        final double s = getStandardDeviation();
        return s * s;
    }

    /**
     * {@inheritDoc}
     *
     * The lower bound of the support is always zero
     * no matter the parameters.
     *
     * @return lower bound of the support (always zero)
     *
     */
    @Override
    public double getSupportLowerBound() {
        return 0.0;
    }

    /**
     * {@inheritDoc}
     *
     * The upper bound of the support is always positive infinity
     * no matter the parameters.
     *
     * @return upper bound of the support (always
     * {@code Double.POSITIVE_INFINITY})
     */
    @Override
    public double getSupportUpperBound() {
        return Double.POSITIVE_INFINITY;
    }

    /**
     * {@inheritDoc}
     *
     * The support of this distribution is connected.
     *
     * @return {@code true}
     */
    @Override
    public boolean isSupportConnected() {
        return true;
    }
}

