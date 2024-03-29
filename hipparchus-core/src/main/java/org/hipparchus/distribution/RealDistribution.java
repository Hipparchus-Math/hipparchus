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
package org.hipparchus.distribution;

import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * Base interface for continuous distributions.
 */
public interface RealDistribution {

    /**
     * For a random variable {@code X} whose values are distributed according
     * to this distribution, this method returns {@code P(x0 < X <= x1)}.
     *
     * @param x0 the exclusive lower bound
     * @param x1 the inclusive upper bound
     * @return the probability that a random variable with this distribution
     * takes a value between {@code x0} and {@code x1},
     * excluding the lower and including the upper endpoint
     * @throws MathIllegalArgumentException if {@code x0 > x1}
     */
    double probability(double x0, double x1) throws MathIllegalArgumentException;

    /**
     * Returns the probability density function (PDF) of this distribution
     * evaluated at the specified point {@code x}. In general, the PDF is
     * the derivative of the {@link #cumulativeProbability(double) CDF}.
     * If the derivative does not exist at {@code x}, then an appropriate
     * replacement should be returned, e.g. {@code Double.POSITIVE_INFINITY},
     * {@code Double.NaN}, or  the limit inferior or limit superior of the
     * difference quotient.
     *
     * @param x the point at which the PDF is evaluated
     * @return the value of the probability density function at point {@code x}
     */
    double density(double x);

    /**
     * Returns the natural logarithm of the probability density function
     * (PDF) of this distribution evaluated at the specified point {@code x}.
     * In general, the PDF is the derivative of the {@link #cumulativeProbability(double) CDF}.
     * If the derivative does not exist at {@code x}, then an appropriate replacement
     * should be returned, e.g. {@code Double.POSITIVE_INFINITY}, {@code Double.NaN},
     * or the limit inferior or limit superior of the difference quotient. Note that
     * due to the floating point precision and under/overflow issues, this method will
     * for some distributions be more precise and faster than computing the logarithm of
     * {@link #density(double)}.
     *
     * @param x the point at which the PDF is evaluated
     * @return the logarithm of the value of the probability density function at point {@code x}
     */
    double logDensity(double x);

    /**
     * For a random variable {@code X} whose values are distributed according
     * to this distribution, this method returns {@code P(X <= x)}. In other
     * words, this method represents the (cumulative) distribution function
     * (CDF) for this distribution.
     *
     * @param x the point at which the CDF is evaluated
     * @return the probability that a random variable with this
     * distribution takes a value less than or equal to {@code x}
     */
    double cumulativeProbability(double x);

    /**
     * Computes the quantile function of this distribution. For a random
     * variable {@code X} distributed according to this distribution, the
     * returned value is
     * <ul>
     * <li><code>inf{x in R | P(X&lt;=x) &gt;= p}</code> for {@code 0 < p <= 1},</li>
     * <li><code>inf{x in R | P(X&lt;=x) &gt; 0}</code> for {@code p = 0}.</li>
     * </ul>
     *
     * @param p the cumulative probability
     * @return the smallest {@code p}-quantile of this distribution
     * (largest 0-quantile for {@code p = 0})
     * @throws MathIllegalArgumentException if {@code p < 0} or {@code p > 1}
     */
    double inverseCumulativeProbability(double p) throws MathIllegalArgumentException;

    /**
     * Use this method to get the numerical value of the mean of this
     * distribution.
     *
     * @return the mean or {@code Double.NaN} if it is not defined
     */
    double getNumericalMean();

    /**
     * Use this method to get the numerical value of the variance of this
     * distribution.
     *
     * @return the variance (possibly {@code Double.POSITIVE_INFINITY} as
     * for certain cases in {@link org.hipparchus.distribution.continuous.TDistribution})
     * or {@code Double.NaN} if it is not defined
     */
    double getNumericalVariance();

    /**
     * Access the lower bound of the support. This method must return the same
     * value as {@code inverseCumulativeProbability(0)}. In other words, this
     * method must return
     * <p><code>inf {x in R | P(X &lt;= x) &gt; 0}</code>.</p>
     *
     * @return lower bound of the support (might be
     * {@code Double.NEGATIVE_INFINITY})
     */
    double getSupportLowerBound();

    /**
     * Access the upper bound of the support. This method must return the same
     * value as {@code inverseCumulativeProbability(1)}. In other words, this
     * method must return
     * <p><code>inf {x in R | P(X &lt;= x) = 1}</code>.</p>
     *
     * @return upper bound of the support (might be
     * {@code Double.POSITIVE_INFINITY})
     */
    double getSupportUpperBound();

    /**
     * Use this method to get information about whether the support is connected,
     * i.e. whether all values between the lower and upper bound of the support
     * are included in the support.
     *
     * @return whether the support is connected or not
     */
    boolean isSupportConnected();

}
