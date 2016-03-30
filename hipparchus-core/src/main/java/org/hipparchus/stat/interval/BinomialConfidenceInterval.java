/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.stat.interval;

import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * Interface to generate confidence intervals for a binomial proportion.
 *
 * @see
 * <a href="http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval">
 * Binomial proportion confidence interval (Wikipedia)</a>
 */
public interface BinomialConfidenceInterval {

    /**
     * Returns the Agresti-Coull method for creating a binomial proportion confidence interval.
     * <p>
     * The returns instance is thread-safe.
     */
    static BinomialConfidenceInterval agrestiCoull() {
        return BinomialConfidenceIntervals.AgrestiCoullInterval.INSTANCE;
    }

    /**
     * Returns the normal approximation method for creating a binomial proportion confidence interval.
     * <p>
     * The returns instance is thread-safe.
     */
    static BinomialConfidenceInterval normalApproximation() {
        return BinomialConfidenceIntervals.NormalApproximationInterval.INSTANCE;
    }

    /**
     * Returns the Clopper-Pearson method for creating a binomial proportion confidence interval.
     * <p>
     * The returns instance is thread-safe.
     */
    static BinomialConfidenceInterval clopperPearson() {
        return BinomialConfidenceIntervals.ClopperPearsonInterval.INSTANCE;
    }

    /**
     * Returns the Wilson score method for creating a binomial proportion confidence interval.
     * <p>
     * The returns instance is thread-safe.
     */
    static BinomialConfidenceInterval wilsonScore() {
        return BinomialConfidenceIntervals.WilsonScoreInterval.INSTANCE;
    }

    /**
     * Create a confidence interval for the true probability of success
     * of an unknown binomial distribution with the given observed number
     * of trials, probability of success and confidence level.
     * <p>
     * Preconditions:
     * <ul>
     * <li>{@code numberOfTrials} must be positive</li>
     * <li>{@code probabilityOfSuccess} must be between 0 and 1 (inclusive)</li>
     * <li>{@code confidenceLevel} must be strictly between 0 and 1 (exclusive)</li>
     * </ul>
     *
     * @param numberOfTrials number of trials
     * @param probabilityOfSuccess observed probability of success
     * @param confidenceLevel desired probability that the true probability of
     * success falls within the returned interval
     * @return Confidence interval containing the probability of success with
     * probability {@code confidenceLevel}
     * @throws MathIllegalArgumentException if {@code numberOfTrials <= 0}.
     * @throws MathIllegalArgumentException if {@code probabilityOfSuccess is not in the interval [0, 1]}.
     * @throws MathIllegalArgumentException if {@code confidenceLevel} is not in the interval {@code (0, 1)}.
     */
    ConfidenceInterval createInterval(int numberOfTrials,
                                      double probabilityOfSuccess,
                                      double confidenceLevel) throws MathIllegalArgumentException;

}
