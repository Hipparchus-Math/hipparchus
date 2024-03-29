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
package org.hipparchus.distribution.discrete;

import org.hipparchus.distribution.continuous.NormalDistribution;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.special.Gamma;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/**
 * Implementation of the Poisson distribution.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Poisson_distribution">Poisson distribution (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/PoissonDistribution.html">Poisson distribution (MathWorld)</a>
 */
public class PoissonDistribution extends AbstractIntegerDistribution {
    /** Default maximum number of iterations for cumulative probability calculations. */
    public static final int DEFAULT_MAX_ITERATIONS = 10000000;
    /** Default convergence criterion. */
    public static final double DEFAULT_EPSILON = 1e-12;
    /** Serializable version identifier. */
    private static final long serialVersionUID = 20160320L;
    /** Distribution used to compute normal approximation. */
    private final NormalDistribution normal;
    /** Mean of the distribution. */
    private final double mean;

    /**
     * Maximum number of iterations for cumulative probability. Cumulative
     * probabilities are estimated using either Lanczos series approximation
     * of {@link Gamma#regularizedGammaP(double, double, double, int)}
     * or continued fraction approximation of
     * {@link Gamma#regularizedGammaQ(double, double, double, int)}.
     */
    private final int maxIterations;

    /** Convergence criterion for cumulative probability. */
    private final double epsilon;

    /**
     * Creates a new Poisson distribution with specified mean.
     *
     * @param p the Poisson mean
     * @throws MathIllegalArgumentException if {@code p <= 0}.
     */
    public PoissonDistribution(double p) throws MathIllegalArgumentException {
        this(p, DEFAULT_EPSILON, DEFAULT_MAX_ITERATIONS);
    }

    /**
     * Creates a new Poisson distribution with specified mean, convergence
     * criterion and maximum number of iterations.
     *
     * @param p Poisson mean.
     * @param epsilon Convergence criterion for cumulative probabilities.
     * @param maxIterations the maximum number of iterations for cumulative
     * probabilities.
     * @throws MathIllegalArgumentException if {@code p <= 0}.
     */
    public PoissonDistribution(double p, double epsilon, int maxIterations)
        throws MathIllegalArgumentException {
        if (p <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.MEAN, p);
        }
        mean = p;
        this.epsilon = epsilon;
        this.maxIterations = maxIterations;

        // Use the same RNG instance as the parent class.
        normal = new NormalDistribution(p, FastMath.sqrt(p));
    }

    /**
     * Creates a new Poisson distribution with the specified mean and
     * convergence criterion.
     *
     * @param p Poisson mean.
     * @param epsilon Convergence criterion for cumulative probabilities.
     * @throws MathIllegalArgumentException if {@code p <= 0}.
     */
    public PoissonDistribution(double p, double epsilon)
        throws MathIllegalArgumentException {
        this(p, epsilon, DEFAULT_MAX_ITERATIONS);
    }

    /**
     * Creates a new Poisson distribution with the specified mean and maximum
     * number of iterations.
     *
     * @param p Poisson mean.
     * @param maxIterations Maximum number of iterations for cumulative probabilities.
     */
    public PoissonDistribution(double p, int maxIterations) {
        this(p, DEFAULT_EPSILON, maxIterations);
    }

    /**
     * Get the mean for the distribution.
     *
     * @return the mean for the distribution.
     */
    public double getMean() {
        return mean;
    }

    /** {@inheritDoc} */
    @Override
    public double probability(int x) {
        final double logProbability = logProbability(x);
        return logProbability == Double.NEGATIVE_INFINITY ? 0 : FastMath.exp(logProbability);
    }

    /** {@inheritDoc} */
    @Override
    public double logProbability(int x) {
        double ret;
        if (x < 0 || x == Integer.MAX_VALUE) {
            ret = Double.NEGATIVE_INFINITY;
        } else if (x == 0) {
            ret = -mean;
        } else {
            ret = -SaddlePointExpansion.getStirlingError(x) -
                  SaddlePointExpansion.getDeviancePart(x, mean) -
                  0.5 * FastMath.log(MathUtils.TWO_PI) - 0.5 * FastMath.log(x);
        }
        return ret;
    }

    /** {@inheritDoc} */
    @Override
    public double cumulativeProbability(int x) {
        if (x < 0) {
            return 0;
        }
        if (x == Integer.MAX_VALUE) {
            return 1;
        }
        return Gamma.regularizedGammaQ((double) x + 1, mean, epsilon,
                                       maxIterations);
    }

    /**
     * Calculates the Poisson distribution function using a normal
     * approximation. The {@code N(mean, sqrt(mean))} distribution is used
     * to approximate the Poisson distribution. The computation uses
     * "half-correction" (evaluating the normal distribution function at
     * {@code x + 0.5}).
     *
     * @param x Upper bound, inclusive.
     * @return the distribution function value calculated using a normal
     * approximation.
     */
    public double normalApproximateProbability(int x)  {
        // calculate the probability using half-correction
        return normal.cumulativeProbability(x + 0.5);
    }

    /**
     * {@inheritDoc}
     *
     * For mean parameter {@code p}, the mean is {@code p}.
     */
    @Override
    public double getNumericalMean() {
        return getMean();
    }

    /**
     * {@inheritDoc}
     *
     * For mean parameter {@code p}, the variance is {@code p}.
     */
    @Override
    public double getNumericalVariance() {
        return getMean();
    }

    /**
     * {@inheritDoc}
     *
     * The lower bound of the support is always 0 no matter the mean parameter.
     *
     * @return lower bound of the support (always 0)
     */
    @Override
    public int getSupportLowerBound() {
        return 0;
    }

    /**
     * {@inheritDoc}
     *
     * The upper bound of the support is positive infinity,
     * regardless of the parameter values. There is no integer infinity,
     * so this method returns {@code Integer.MAX_VALUE}.
     *
     * @return upper bound of the support (always {@code Integer.MAX_VALUE} for
     * positive infinity)
     */
    @Override
    public int getSupportUpperBound() {
        return Integer.MAX_VALUE;
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
