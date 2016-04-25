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
package org.hipparchus.stat.descriptive;

import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.stat.descriptive.moment.GeometricMean;
import org.hipparchus.stat.descriptive.moment.Mean;
import org.hipparchus.stat.descriptive.moment.SecondMoment;
import org.hipparchus.stat.descriptive.moment.Variance;
import org.hipparchus.stat.descriptive.rank.Max;
import org.hipparchus.stat.descriptive.rank.Min;
import org.hipparchus.stat.descriptive.summary.Sum;
import org.hipparchus.stat.descriptive.summary.SumOfLogs;
import org.hipparchus.stat.descriptive.summary.SumOfSquares;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/**
 * Computes summary statistics for a stream of data values added using the
 * {@link #addValue(double) addValue} method. The data values are not stored in
 * memory, so this class can be used to compute statistics for very large data
 * streams.
 * <p>
 * To create a SummaryStatistics instance with default statistic implementations,
 * use:
 * <pre>
 *    SummaryStatistics stats = SummaryStatistics.create();
 * </pre>
 * <p>
 * The {@link StorelessUnivariateStatistic} instances used to maintain summary
 * state and compute statistics are configurable via a builder. For example, the
 * default implementation for the variance can be overridden by calling
 * <pre>
 *    SummaryStatistics stats =
 *        SummaryStatistics.builder()
 *                         .withVariance(new MyVarianceImpl())
 *                         .build();
 * </pre>
 * <p>
 * Note: This class is not thread-safe. Use
 * {@link #synchronizedSummaryStatistics(SummaryStatistics)} if concurrent
 * access from multiple threads is required.
 */
public interface SummaryStatistics extends StatisticalSummary {

    /**
     * Returns a builder for a {@link SummaryStatistics}.
     *
     * @return a summary statistics builder.
     */
    static Builder builder() {
        return new Builder();
    }

    /**
     * Returns a new SummaryStatistics instance using default implementations
     * for all statistical values.
     * <p>
     * The returned instance is <b>not</b> thread-safe.
     */
    static SummaryStatistics create() {
        return builder().build();
    }

    /**
     * Decorates another SummaryStatistics to synchronize its behaviour for a
     * multi-threaded environment.
     * <p>
     * Methods are synchronized, then forwarded to the decorated SummaryStatistics.
     *
     * @param statistic the summary statistic to decorate
     * @return a synchronized wrapper for the summary statistic
     */
    static SummaryStatistics synchronizedSummaryStatistics(SummaryStatistics statistic) {
        if (statistic instanceof SynchronizedSummaryStatistics) {
            return statistic;
        }
        return new SynchronizedSummaryStatistics(statistic);
    }

    /**
     * Return a {@link StatisticalSummaryValues} instance reporting current
     * statistics.
     * @return Current values of statistics
     */
    default StatisticalSummary getSummary() {
        return new StatisticalSummaryValues(getMean(), getVariance(), getN(),
                                            getMax(), getMin(), getSum());
    }

    /**
     * Returns a copy of this SummaryStatistics instance with the same internal state.
     *
     * @return a copy of this
     */
    SummaryStatistics copy();

    /**
     * Add a value to the data
     * @param value the value to add
     */
    void addValue(double value);

    /**
     * Resets all statistics and storage.
     */
    void clear();

    /**
     * Returns the quadratic mean, a.k.a.
     * <a href="http://mathworld.wolfram.com/Root-Mean-Square.html">
     * root-mean-square</a> of the available values
     *
     * @return The quadratic mean or {@code Double.NaN} if no values
     * have been added.
     */
    default double getQuadraticMean() {
        final long size = getN();
        return size > 0 ? FastMath.sqrt(getSumOfSquares() / size) : Double.NaN;
    }

    /**
     * Returns the standard deviation of the values that have been added.
     * <p>
     * Double.NaN is returned if no values have been added.
     *
     * @return the standard deviation
     */
    @Override
    default double getStandardDeviation() {
        final long size = getN();
        if (size > 0) {
            return size > 1 ? FastMath.sqrt(getVariance()) : 0.0;
        } else {
            return Double.NaN;
        }
    }

    /**
     * Returns the <a href="http://en.wikibooks.org/wiki/Statistics/Summary/Variance">
     * population variance</a> of the values that have been added.
     * <p>
     * Double.NaN is returned if no values have been added.
     *
     * @return the population variance
     */
    double getPopulationVariance();

    /**
     * Returns the geometric mean of the values that have been added.
     * <p>
     * Double.NaN is returned if no values have been added.
     *
     * @return the geometric mean
     */
    double getGeometricMean();

    /**
     * Returns the sum of the squares of the values that have been added.
     * <p>
     * Double.NaN is returned if no values have been added.
     *
     * @return The sum of squares
     */
    double getSumOfSquares();

    /**
     * Returns the sum of the logs of the values that have been added.
     * <p>
     * Double.NaN is returned if no values have been added.
     *
     * @return the sum of logs
     */
    double getSumOfLogs();

    /**
     * Returns a statistic related to the Second Central Moment. Specifically,
     * what is returned is the sum of squared deviations from the sample mean
     * among the values that have been added.
     * <p>
     * Returns <code>Double.NaN</code> if no data values have been added and
     * returns <code>0</code> if there is just one value in the data set.
     *
     * @return second central moment statistic
     */
    double getSecondMoment();

    /**
     * A mutable builder for a SummaryStatistics.
     */
    public static class Builder {
        /** Second moment can not be set externally, needed for other statistics. */
        final SecondMoment secondMoment = new SecondMoment();

        StorelessUnivariateStatistic maxImpl;
        StorelessUnivariateStatistic minImpl;
        StorelessUnivariateStatistic sumImpl;
        StorelessUnivariateStatistic sumOfSquaresImpl;
        StorelessUnivariateStatistic sumOfLogsImpl;
        StorelessUnivariateStatistic meanImpl;
        StorelessUnivariateStatistic varianceImpl;
        StorelessUnivariateStatistic geoMeanImpl;

        boolean meanUsesExternalMoment;
        boolean varianceUsesExternalMoment;
        boolean geometricMeanUsesExternalSumOfLogs;

        boolean threadsafe;

        Builder() {}

        /**
         * Sets the max implementation to use.
         *
         * @param impl the max implementation
         * @return the builder
         * @throws NullArgumentException if impl is null
         */
        public Builder withMaxImpl(StorelessUnivariateStatistic impl) {
            MathUtils.checkNotNull(impl);
            this.maxImpl = impl;
            return this;
        }

        /**
         * Sets the min implementation to use.
         *
         * @param impl the min implementation
         * @return the builder
         * @throws NullArgumentException if impl is null
         */
        public Builder withMinImpl(StorelessUnivariateStatistic impl) {
            MathUtils.checkNotNull(impl);
            this.minImpl = impl;
            return this;
        }

        /**
         * Sets the mean implementation to use.
         *
         * @param impl the mean implementation
         * @return the builder
         * @throws NullArgumentException if impl is null
         */
        public Builder withMeanImpl(StorelessUnivariateStatistic impl) {
            MathUtils.checkNotNull(impl);
            this.meanImpl = impl;
            return this;
        }

        /**
         * Sets the geometric mean implementation to use.
         *
         * @param impl the geometric mean implementation
         * @return the builder
         * @throws NullArgumentException if impl is null
         */
        public Builder withGeometricMeanImpl(StorelessUnivariateStatistic impl) {
            MathUtils.checkNotNull(impl);
            this.geoMeanImpl = impl;
            return this;
        }

        /**
         * Sets the variance implementation to use.
         *
         * @param impl the variance implementation
         * @return the builder
         * @throws NullArgumentException if impl is null
         */
        public Builder withVarianceImpl(StorelessUnivariateStatistic impl) {
            MathUtils.checkNotNull(impl);
            this.varianceImpl = impl;
            return this;
        }

        /**
         * Sets the sum implementation to use.
         *
         * @param impl the sum implementation
         * @return the builder
         * @throws NullArgumentException if impl is null
         */
        public Builder withSumImpl(StorelessUnivariateStatistic impl) {
            MathUtils.checkNotNull(impl);
            this.sumImpl = impl;
            return this;
        }

        /**
         * Sets the sum of squares implementation to use.
         *
         * @param impl the sum of squares implementation
         * @return the builder
         * @throws NullArgumentException if impl is null
         */
        public Builder withSumOfSquaresImpl(StorelessUnivariateStatistic impl) {
            MathUtils.checkNotNull(impl);
            this.sumOfSquaresImpl = impl;
            return this;
        }

        /**
         * Sets the sum of logs implementation to use.
         *
         * @param impl the sum of logs implementation
         * @return the builder
         * @throws NullArgumentException if impl is null
         */
        public Builder withSumOfLogsImpl(StorelessUnivariateStatistic impl) {
            MathUtils.checkNotNull(impl);
            this.sumOfLogsImpl = impl;
            return this;
        }

        /**
         * Indicates if the returned instance shall be thread-safe,
         * i.e. all access is synchronized.
         *
         * @return the builder
         */
        public Builder threadsafe() {
            this.threadsafe = true;
            return this;
        }

        /**
         * Sets up default implementations to use.
         */
        protected void setupDefaultsImpls() {
            if (maxImpl == null) {
                maxImpl = new Max();
            }
            if (minImpl == null) {
                minImpl = new Min();
            }
            if (sumImpl == null) {
                sumImpl = new Sum();
            }
            if (sumOfSquaresImpl == null) {
                sumOfSquaresImpl = new SumOfSquares();
            }
            if (sumOfLogsImpl == null) {
                sumOfLogsImpl = new SumOfLogs();
            }
            if (meanImpl == null) {
                meanImpl = new Mean(secondMoment);
                meanUsesExternalMoment = true;
            }
            if (varianceImpl == null) {
                varianceImpl = new Variance(secondMoment);
                varianceUsesExternalMoment = true;
            }
            if (geoMeanImpl == null) {
                if (sumOfLogsImpl instanceof SumOfLogs) {
                    geoMeanImpl = new GeometricMean((SumOfLogs) sumOfLogsImpl);
                    geometricMeanUsesExternalSumOfLogs = true;
                } else {
                    geoMeanImpl = new GeometricMean();
                }
            }
        }

        /**
         * Constructs a new SummaryStatistics instance with the values
         * stored in this builder.
         *
         * @return a new SummaryStatistics instance.
         */
        public SummaryStatistics build() {
            // setup default implementations.
            setupDefaultsImpls();
            final SummaryStatistics instance = new SummaryStatisticsImpl(this);
            return threadsafe ? synchronizedSummaryStatistics(instance) : instance;
        }
    }

}
