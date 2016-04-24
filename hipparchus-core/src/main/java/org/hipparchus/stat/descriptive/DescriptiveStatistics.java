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

import java.util.Arrays;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.stat.descriptive.moment.GeometricMean;
import org.hipparchus.stat.descriptive.moment.Kurtosis;
import org.hipparchus.stat.descriptive.moment.Mean;
import org.hipparchus.stat.descriptive.moment.Skewness;
import org.hipparchus.stat.descriptive.moment.Variance;
import org.hipparchus.stat.descriptive.rank.Max;
import org.hipparchus.stat.descriptive.rank.Min;
import org.hipparchus.stat.descriptive.rank.Percentile;
import org.hipparchus.stat.descriptive.summary.Sum;
import org.hipparchus.stat.descriptive.summary.SumOfSquares;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;


/**
 * Maintains a dataset of values of a single variable and computes descriptive
 * statistics based on stored data.
 * <p>
 * The {@link #getWindowSize() windowSize} property sets a limit on the number
 * of values that can be stored in the dataset. The default value, INFINITE_WINDOW,
 * puts no limit on the size of the dataset. This value should be used with
 * caution, as the backing store will grow without bound in this case.
 * <p>
 * For very large datasets, {@link SummaryStatistics}, which does not store
 * the dataset, should be used instead of this class. If <code>windowSize</code>
 * is not INFINITE_WINDOW and more values are added than can be stored in the
 * dataset, new values are added in a "rolling" manner, with new values replacing
 * the "oldest" values in the dataset.
 * <p>
 * To create a DescriptiveStatistics instance with default statistic implementations,
 * use:
 * <pre>
 *    DescriptiveStatistics stats = DescriptiveStatistics.create();
 * </pre>
 * <p>
 * The {@link UnivariateStatistic} instances used to compute statistics
 * are configurable via a builder. For example, the default implementation for
 * the variance can be overridden by calling
 * <pre>
 *    DescriptiveStatistics stats =
 *        DescriptiveStatistics.builder()
 *                             .withVariance(new MyVarianceImpl())
 *                             .build();
 * </pre>
 * <p>
 * Note: this class is not threadsafe.
 * Use {@link #synchronizedDescriptiveStatistics(DescriptiveStatistics)}
 * if concurrent access from multiple threads is required.
 */
public interface DescriptiveStatistics extends StatisticalSummary {

    /**
     * Represents an infinite window size.  When the {@link #getWindowSize()}
     * returns this value, there is no limit to the number of data values
     * that can be stored in the dataset.
     */
    static final int INFINITE_WINDOW = -1;

    /**
     * Returns a builder for a {@link DescriptiveStatistics}.
     *
     * @return a descriptive statistics builder.
     */
    static Builder builder() {
        return new Builder();
    }

    /**
     * Construct a DescriptiveStatistics instance with an infinite window.
     */
    static DescriptiveStatistics create() {
        return builder().build();
    }

    /**
     * Construct a DescriptiveStatistics instance with the specified window.
     *
     * @param window the window size.
     * @throws MathIllegalArgumentException if window size is less than 1 but
     * not equal to {@link #INFINITE_WINDOW}
     */
    static DescriptiveStatistics of(int window) throws MathIllegalArgumentException {
        return builder().withWindowSize(window).build();
    }

    /**
     * Construct a DescriptiveStatistics instance with an infinite window
     * and the initial data values in double[] initialDoubleArray.
     *
     * @param initialDoubleArray the initial double[].
     * @throws NullArgumentException if the input array is null
     */
    static DescriptiveStatistics of(double[] initialDoubleArray) {
        return builder().withInitialValues(initialDoubleArray).build();
    }

    /**
     * Decorates another DescriptiveStatistics to synchronize its behaviour for a
     * multi-threaded environment.
     * <p>
     * Methods are synchronized, then forwarded to the decorated DescriptiveStatistics.
     *
     * @param statistic the summary statistic to decorate
     * @return a synchronized wrapper for the descriptive statistic
     */
    static DescriptiveStatistics synchronizedDescriptiveStatistics(DescriptiveStatistics statistic) {
        if (statistic instanceof SynchronizedDescriptiveStatistics) {
            return statistic;
        }
        return new SynchronizedDescriptiveStatistics(statistic);
    }

    /**
     * Adds the value to the dataset. If the dataset is at the maximum size
     * (i.e., the number of stored elements equals the currently configured
     * windowSize), the first (oldest) element in the dataset is discarded
     * to make room for the new value.
     *
     * @param v the value to be added
     */
    void addValue(double v);

    /**
     * Removes the most recent value from the dataset.
     *
     * @throws MathIllegalStateException if there are no elements stored
     */
    void removeMostRecentValue() throws MathIllegalStateException;

    /**
     * Replaces the most recently stored value with the given value.
     * There must be at least one element stored to call this method.
     *
     * @param v the value to replace the most recent stored value
     * @return replaced value
     * @throws MathIllegalStateException if there are no elements stored
     */
    double replaceMostRecentValue(double v) throws MathIllegalStateException;

    /**
     * Returns the geometric mean of the available values.
     * <p>
     * See {@link GeometricMean} for details on the computing algorithm.
     *
     * @see <a href="http://www.xycoon.com/geometric_mean.htm">
     * Geometric mean</a>
     *
     * @return The geometricMean, Double.NaN if no values have been added,
     * or if any negative values have been added.
     */
    double getGeometricMean();

    /**
     * Returns the population variance of the available values.
     *
     * @see <a href="http://en.wikibooks.org/wiki/Statistics/Summary/Variance">
     * Population variance</a>
     *
     * @return The population variance, Double.NaN if no values have been added,
     * or 0.0 for a single value set.
     */
    double getPopulationVariance();

    /**
     * Returns the standard deviation of the available values.
     * @return The standard deviation, Double.NaN if no values have been added
     * or 0.0 for a single value set.
     */
    @Override
    default double getStandardDeviation() {
        double stdDev = Double.NaN;
        if (getN() > 0) {
            if (getN() > 1) {
                stdDev = FastMath.sqrt(getVariance());
            } else {
                stdDev = 0.0;
            }
        }
        return stdDev;
    }

    /**
     * Returns the quadratic mean of the available values.
     *
     * @see <a href="http://mathworld.wolfram.com/Root-Mean-Square.html">
     * Root Mean Square</a>
     *
     * @return The quadratic mean or {@code Double.NaN} if no values
     * have been added.
     */
    default double getQuadraticMean() {
        final long n = getN();
        return n > 0 ? FastMath.sqrt(getSumOfSquares() / n) : Double.NaN;
    }

    /**
     * Returns the skewness of the available values. Skewness is a
     * measure of the asymmetry of a given distribution.
     *
     * @return The skewness, Double.NaN if less than 3 values have been added.
     */
    double getSkewness();

    /**
     * Returns the Kurtosis of the available values. Kurtosis is a
     * measure of the "peakedness" of a distribution.
     *
     * @return The kurtosis, Double.NaN if less than 4 values have been added.
     */
    double getKurtosis();

    /**
     * Returns the sum of the squares of the available values.
     * @return The sum of the squares or Double.NaN if no
     * values have been added.
     */
    double getSumOfSquares();

    /**
     * Returns an estimate for the pth percentile of the stored values.
     * <p>
     * The implementation provided here follows the first estimation procedure presented
     * <a href="http://www.itl.nist.gov/div898/handbook/prc/section2/prc252.htm">here.</a>
     * </p><p>
     * <strong>Preconditions</strong>:<ul>
     * <li><code>0 &lt; p &le; 100</code> (otherwise an
     * <code>MathIllegalArgumentException</code> is thrown)</li>
     * <li>at least one value must be stored (returns <code>Double.NaN
     *     </code> otherwise)</li>
     * </ul>
     *
     * @param p the requested percentile (scaled from 0 - 100)
     * @return An estimate for the pth percentile of the stored data
     * @throws MathIllegalArgumentException if p is not a valid quantile
     */
    double getPercentile(final double p)
        throws MathIllegalArgumentException;

    /**
     * Resets all statistics and storage.
     */
    void clear();

    /**
     * Returns the maximum number of values that can be stored in the
     * dataset, or INFINITE_WINDOW (-1) if there is no limit.
     *
     * @return The current window size or -1 if its Infinite.
     */
    int getWindowSize();

    /**
     * WindowSize controls the number of values that contribute to the
     * reported statistics.  For example, if windowSize is set to 3 and the
     * values {1,2,3,4,5} have been added <strong> in that order</strong> then
     * the <i>available values</i> are {3,4,5} and all reported statistics will
     * be based on these values. If {@code windowSize} is decreased as a result
     * of this call and there are more than the new value of elements in the
     * current dataset, values from the front of the array are discarded to
     * reduce the dataset to {@code windowSize} elements.
     *
     * @param windowSize sets the size of the window.
     * @throws MathIllegalArgumentException if window size is less than 1 but
     * not equal to {@link #INFINITE_WINDOW}
     */
    void setWindowSize(int windowSize)
        throws MathIllegalArgumentException;

    /**
     * Returns the current set of values in an array of double primitives.
     * The order of addition is preserved.  The returned array is a fresh
     * copy of the underlying data -- i.e., it is not a reference to the
     * stored data.
     *
     * @return the current set of numbers in the order in which they
     * were added to this set
     */
    double[] getValues();

    /**
     * Returns the current set of values in an array of double primitives,
     * sorted in ascending order.  The returned array is a fresh
     * copy of the underlying data -- i.e., it is not a reference to the
     * stored data.
     * @return returns the current set of
     * numbers sorted in ascending order
     */
    default double[] getSortedValues() {
        double[] sort = getValues();
        Arrays.sort(sort);
        return sort;
    }

    /**
     * Returns the element at the specified index
     * @param index The Index of the element
     * @return return the element at the specified index
     */
    double getElement(int index);

    /**
     * Apply the given statistic to the data associated with this set of statistics.
     * @param stat the statistic to apply
     * @return the computed value of the statistic.
     */
    double apply(UnivariateStatistic stat);

    /**
     * Returns a copy of this DescriptiveStatistics instance with the same internal state.
     *
     * @return a copy of this
     */
    DescriptiveStatistics copy();

    /**
     * A mutable builder for a DescriptiveStatistics.
     */
    public static class Builder {
        UnivariateStatistic          maxImpl           = new Max();
        UnivariateStatistic          minImpl           = new Min();
        UnivariateStatistic          sumImpl           = new Sum();
        UnivariateStatistic          sumsqImpl         = new SumOfSquares();
        UnivariateStatistic          meanImpl          = new Mean();
        UnivariateStatistic          varianceImpl      = new Variance();
        UnivariateStatistic          geometricMeanImpl = new GeometricMean();
        UnivariateStatistic          kurtosisImpl      = new Kurtosis();
        UnivariateStatistic          skewnessImpl      = new Skewness();
        QuantiledUnivariateStatistic percentileImpl    = new Percentile();

        int      windowSize = INFINITE_WINDOW;
        double[] initialValues;
        boolean  threadsafe;

        protected Builder() {}

        /**
         * Sets the max implementation to use.
         *
         * @param maxImpl the max implementation
         * @return the builder
         * @throws NullArgumentException if maxImpl is null
         */
        public Builder withMaxImpl(UnivariateStatistic maxImpl) {
            MathUtils.checkNotNull(maxImpl);
            this.maxImpl = maxImpl;
            return this;
        }

        /**
         * Sets the min implementation to use.
         *
         * @param minImpl the min implementation
         * @return the builder
         * @throws NullArgumentException if minImpl is null
         */
        public Builder withMinImpl(UnivariateStatistic minImpl) {
            MathUtils.checkNotNull(minImpl);
            this.minImpl = minImpl;
            return this;
        }

        /**
         * Sets the mean implementation to use.
         *
         * @param meanImpl the mean implementation
         * @return the builder
         * @throws NullArgumentException if meanImpl is null
         */
        public Builder withMeanImpl(UnivariateStatistic meanImpl) {
            MathUtils.checkNotNull(meanImpl);
            this.meanImpl = meanImpl;
            return this;
        }

        /**
         * Sets the geometric mean implementation to use.
         *
         * @param geometricMeanImpl the geometric mean implementation
         * @return the builder
         * @throws NullArgumentException if geometricMeanImpl is null
         */
        public Builder withGeometricMeanImpl(UnivariateStatistic geometricMeanImpl) {
            MathUtils.checkNotNull(geometricMeanImpl);
            this.geometricMeanImpl = geometricMeanImpl;
            return this;
        }

        /**
         * Sets the kurtosis implementation to use.
         *
         * @param kurtosisImpl the kurtosis implementation
         * @return the builder
         * @throws NullArgumentException if kurtosisImpl is null
         */
        public Builder withKurtosisImpl(UnivariateStatistic kurtosisImpl) {
            MathUtils.checkNotNull(kurtosisImpl);
            this.kurtosisImpl = kurtosisImpl;
            return this;
        }

        /**
         * Sets the percentile implementation to use.
         *
         * @param percentileImpl the percentile implementation
         * @return the builder
         * @throws NullArgumentException if percentileImpl is null
         */
        public Builder withPercentileImpl(QuantiledUnivariateStatistic percentileImpl) {
            MathUtils.checkNotNull(percentileImpl);
            this.percentileImpl = percentileImpl;
            return this;
        }

        /**
         * Sets the skewness implementation to use.
         *
         * @param skewnessImpl the skewness implementation
         * @return the builder
         * @throws NullArgumentException if skewnessImpl is null
         */
        public Builder withSkewnessImpl(UnivariateStatistic skewnessImpl) {
            MathUtils.checkNotNull(skewnessImpl);
            this.skewnessImpl = skewnessImpl;
            return this;
        }

        /**
         * Sets the variance implementation to use.
         *
         * @param varianceImpl the variance implementation
         * @return the builder
         * @throws NullArgumentException if varianceImpl is null
         */
        public Builder withVarianceImpl(UnivariateStatistic varianceImpl) {
            MathUtils.checkNotNull(varianceImpl);
            this.varianceImpl = varianceImpl;
            return this;
        }

        /**
         * Sets the sum of squares implementation to use.
         *
         * @param sumSqImpl the sum of squares implementation
         * @return the builder
         * @throws NullArgumentException if sumSqImpl is null
         */
        public Builder withSumOfSquaresImpl(UnivariateStatistic sumSqImpl) {
            MathUtils.checkNotNull(sumSqImpl);
            this.sumsqImpl = sumSqImpl;
            return this;
        }

        /**
         * Sets the sum implementation to use.
         *
         * @param sumImpl the sum implementation
         * @return the builder
         * @throws NullArgumentException if sumImpl is null
         */
        public Builder withSumImpl(UnivariateStatistic sumImpl) {
            MathUtils.checkNotNull(sumImpl);
            this.sumImpl = sumImpl;
            return this;
        }

        /**
         * Sets the window of values that contribute to the
         * reported statistics.
         *
         * @see DescriptiveStatistics#setWindowSize(int)
         *
         * @param windowSize sets the size of the window.
         * @throws MathIllegalArgumentException if window size is less than 1 but
         * not equal to {@link #INFINITE_WINDOW}
         * @return the builder
         */
        public Builder withWindowSize(int windowSize) {
            if (windowSize < 1 && windowSize != INFINITE_WINDOW) {
                throw new MathIllegalArgumentException(
                        LocalizedCoreFormats.NOT_POSITIVE_WINDOW_SIZE, windowSize);
            }

            this.windowSize = windowSize;
            return this;
        }

        /**
         * Sets the initial values of the reported statistics.
         *
         * @param initialValues the initial values to use.
         * @throws NullArgumentException if initialValues is null
         * @return the builder
         */
        public Builder withInitialValues(double[] initialValues) {
            MathUtils.checkNotNull(initialValues, LocalizedCoreFormats.INPUT_ARRAY);
            this.initialValues = initialValues;
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
         * Constructs a new DescriptiveStatistics instance with the values
         * stored in this builder.
         *
         * @return a new DescriptiveStatistics instance.
         */
        public DescriptiveStatistics build() {
            DescriptiveStatistics instance = new DescriptiveStatisticsImpl(this);
            return threadsafe ? synchronizedDescriptiveStatistics(instance) : instance;
        }
    }
}
