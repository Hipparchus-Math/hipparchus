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

package org.hipparchus.distribution.continuous;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.special.Erf;
import org.hipparchus.util.FastMath;

/**
 * Implementation of the log-normal (gaussian) distribution.
 * <p>
 * <strong>Parameters:</strong>
 * {@code X} is log-normally distributed if its natural logarithm {@code log(X)}
 * is normally distributed. The probability distribution function of {@code X}
 * is given by (for {@code x > 0})
 * <p>
 * {@code exp(-0.5 * ((ln(x) - m) / s)^2) / (s * sqrt(2 * pi) * x)}
 * <ul>
 * <li>{@code m} is the <em>location</em> parameter: this is the mean of the
 * normally distributed natural logarithm of this distribution,</li>
 * <li>{@code s} is the <em>shape</em> parameter: this is the standard
 * deviation of the normally distributed natural logarithm of this
 * distribution.
 * </ul>
 *
 * @see <a href="http://en.wikipedia.org/wiki/Log-normal_distribution">
 * Log-normal distribution (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/LogNormalDistribution.html">
 * Log Normal distribution (MathWorld)</a>
 */
public class LogNormalDistribution extends AbstractRealDistribution {

    /** Serializable version identifier. */
    private static final long serialVersionUID = 20120112;

    /** &radic;(2 &pi;) */
    private static final double SQRT2PI = FastMath.sqrt(2 * FastMath.PI);

    /** &radic;(2) */
    private static final double SQRT2 = FastMath.sqrt(2.0);

    /** The location parameter of this distribution (named m in MathWorld and µ in Wikipedia). */
    private final double location;

    /** The shape parameter of this distribution. */
    private final double shape;
    /** The value of {@code log(shape) + 0.5 * log(2*PI)} stored for faster computation. */
    private final double logShapePlusHalfLog2Pi;

    /**
     * Create a log-normal distribution, where the mean and standard deviation
     * of the {@link NormalDistribution normally distributed} natural
     * logarithm of the log-normal distribution are equal to zero and one
     * respectively. In other words, the location of the returned distribution is
     * {@code 0}, while its shape is {@code 1}.
     */
    public LogNormalDistribution() {
        this(0, 1);
    }

    /**
     * Create a log-normal distribution using the specified location and shape.
     *
     * @param location the location parameter of this distribution
     * @param shape the shape parameter of this distribution
     * @throws MathIllegalArgumentException if {@code shape <= 0}.
     */
    public LogNormalDistribution(double location, double shape)
        throws MathIllegalArgumentException {
        this(location, shape, DEFAULT_SOLVER_ABSOLUTE_ACCURACY);
    }


    /**
     * Creates a log-normal distribution.
     *
     * @param location Location parameter of this distribution.
     * @param shape Shape parameter of this distribution.
     * @param inverseCumAccuracy Inverse cumulative probability accuracy.
     * @throws MathIllegalArgumentException if {@code shape <= 0}.
     */
    public LogNormalDistribution(double location,
                                 double shape,
                                 double inverseCumAccuracy)
        throws MathIllegalArgumentException {
        super(inverseCumAccuracy);

        if (shape <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.SHAPE, shape);
        }

        this.location = location;
        this.shape = shape;
        this.logShapePlusHalfLog2Pi = FastMath.log(shape) + 0.5 * FastMath.log(2 * FastMath.PI);
    }

    /**
     * Returns the location parameter of this distribution.
     *
     * @return the location parameter
     * @since 1.4
     */
    public double getLocation() {
        return location;
    }

    /**
     * Returns the shape parameter of this distribution.
     *
     * @return the shape parameter
     */
    public double getShape() {
        return shape;
    }

    /**
     * {@inheritDoc}
     *
     * For location {@code m}, and shape {@code s} of this distribution, the PDF
     * is given by
     * <ul>
     * <li>{@code 0} if {@code x <= 0},</li>
     * <li>{@code exp(-0.5 * ((ln(x) - m) / s)^2) / (s * sqrt(2 * pi) * x)}
     * otherwise.</li>
     * </ul>
     */
    @Override
    public double density(double x) {
        if (x <= 0) {
            return 0;
        }
        final double x0 = FastMath.log(x) - location;
        final double x1 = x0 / shape;
        return FastMath.exp(-0.5 * x1 * x1) / (shape * SQRT2PI * x);
    }

    /** {@inheritDoc}
     *
     * See documentation of {@link #density(double)} for computation details.
     */
    @Override
    public double logDensity(double x) {
        if (x <= 0) {
            return Double.NEGATIVE_INFINITY;
        }
        final double logX = FastMath.log(x);
        final double x0 = logX - location;
        final double x1 = x0 / shape;
        return -0.5 * x1 * x1 - (logShapePlusHalfLog2Pi + logX);
    }

    /**
     * {@inheritDoc}
     *
     * For location {@code m}, and shape {@code s} of this distribution, the CDF
     * is given by
     * <ul>
     * <li>{@code 0} if {@code x <= 0},</li>
     * <li>{@code 0} if {@code ln(x) - m < 0} and {@code m - ln(x) > 40 * s}, as
     * in these cases the actual value is within {@code Double.MIN_VALUE} of 0,
     * <li>{@code 1} if {@code ln(x) - m >= 0} and {@code ln(x) - m > 40 * s},
     * as in these cases the actual value is within {@code Double.MIN_VALUE} of 1,</li>
     * <li>{@code 0.5 + 0.5 * erf((ln(x) - m) / (s * sqrt(2))} otherwise.</li>
     * </ul>
     */
    @Override
    public double cumulativeProbability(double x)  {
        if (x <= 0) {
            return 0;
        }
        final double dev = FastMath.log(x) - location;
        if (FastMath.abs(dev) > 40 * shape) {
            return dev < 0 ? 0.0d : 1.0d;
        }
        return 0.5 + 0.5 * Erf.erf(dev / (shape * SQRT2));
    }

    /** {@inheritDoc} */
    @Override
    public double probability(double x0,
                              double x1)
        throws MathIllegalArgumentException {
        if (x0 > x1) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT,
                                                x0, x1, true);
        }
        if (x0 <= 0 || x1 <= 0) {
            return super.probability(x0, x1);
        }
        final double denom = shape * SQRT2;
        final double v0 = (FastMath.log(x0) - location) / denom;
        final double v1 = (FastMath.log(x1) - location) / denom;
        return 0.5 * Erf.erf(v0, v1);
    }

    /**
     * {@inheritDoc}
     *
     * For location {@code m} and shape {@code s}, the mean is
     * {@code exp(m + s^2 / 2)}.
     */
    @Override
    public double getNumericalMean() {
        double s = shape;
        return FastMath.exp(location + (s * s / 2));
    }

    /**
     * {@inheritDoc}
     *
     * For location {@code m} and shape {@code s}, the variance is
     * {@code (exp(s^2) - 1) * exp(2 * m + s^2)}.
     */
    @Override
    public double getNumericalVariance() {
        final double s = shape;
        final double ss = s * s;
        return (FastMath.expm1(ss)) * FastMath.exp(2 * location + ss);
    }

    /**
     * {@inheritDoc}
     *
     * The lower bound of the support is always 0 no matter the parameters.
     *
     * @return lower bound of the support (always 0)
     */
    @Override
    public double getSupportLowerBound() {
        return 0;
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
