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
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/**
 * Implementation of the Half Cauchy distribution.
 *
 * @see <a href="https://search.r-project.org/CRAN/refmans/LaplacesDemon/html/dist.HalfCauchy.html">Half-Cauchy distribution</a>
 */
public class HalfCauchyDistribution extends AbstractRealDistribution {
    /** Serializable version identifier */
    private static final long serialVersionUID = 20250915L;
    /** The scale of this distribution. */
    private final double scale;

    /**
     * Creates a Half Cauchy distribution.
     *
     * @param scale Scale parameter for this distribution
     * @throws MathIllegalArgumentException if {@code scale <= 0}
     */
    public HalfCauchyDistribution(double scale)
        throws MathIllegalArgumentException {
        if (scale <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.SCALE, scale);
        }
        this.scale = scale;
    }

    /** {@inheritDoc} */
    @Override
    public double cumulativeProbability(double x) {
        if (x < 0) {
            return 0;
        }
        return (2.0/FastMath.PI) * FastMath.atan(x/scale) ;
    }

    /**
     * Access the median.
     *
     * @return the median for this distribution.
     */
    public double getMedian() {
        return Double.NaN;
    }

    /**
     * Access the scale parameter.
     *
     * @return the scale parameter for this distribution.
     */
    public double getScale() {
        return scale;
    }

    /** {@inheritDoc} */
    @Override
    public double density(double x) {
        if (x < 0) {
            return 0;
        }
        return (2.0 / FastMath.PI) * (scale / (x * x + scale * scale));
    }

    /**
     * {@inheritDoc}
     *
     * Returns zero when {@code p == 0}
     * and {@code Double.POSITIVE_INFINITY} when {@code p == 1}.
     */
    @Override
    public double inverseCumulativeProbability(double p) throws MathIllegalArgumentException {
        MathUtils.checkRangeInclusive(p, 0, 1);

        double ret;
        if (p == 0) {
            ret = 0.0;
        } else  if (p == 1) {
            ret = Double.POSITIVE_INFINITY;
        } else {
            ret = scale * FastMath.tan(FastMath.PI * p / 2.0);
        }
        return ret;

    }

    /**
     * {@inheritDoc}
     *
     * The mean is always undefined no matter the parameters.
     *
     * @return mean (always Double.NaN)
     */
    @Override
    public double getNumericalMean() {
        return Double.NaN;
    }

    /**
     * {@inheritDoc}
     *
     * The variance is always undefined no matter the parameters.
     *
     * @return variance (always Double.NaN)
     */
    @Override
    public double getNumericalVariance() {
        return Double.NaN;
    }

    /**
     * {@inheritDoc}
     *
     * The lower bound of the support is always zero no matter
     * the parameters.
     *
     * @return lower bound of the support (always zero)
     */
    @Override
    public double getSupportLowerBound() {
        return 0.0;
    }

    /**
     * {@inheritDoc}
     *
     * The upper bound of the support is always positive infinity no matter
     * the parameters.
     *
     * @return upper bound of the support (always Double.POSITIVE_INFINITY)
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

