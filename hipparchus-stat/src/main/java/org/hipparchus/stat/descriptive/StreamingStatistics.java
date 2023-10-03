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
package org.hipparchus.stat.descriptive;

import java.io.Serializable;
import java.util.function.DoubleConsumer;

import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.stat.descriptive.moment.GeometricMean;
import org.hipparchus.stat.descriptive.moment.Mean;
import org.hipparchus.stat.descriptive.moment.SecondMoment;
import org.hipparchus.stat.descriptive.moment.Variance;
import org.hipparchus.stat.descriptive.rank.Max;
import org.hipparchus.stat.descriptive.rank.Min;
import org.hipparchus.stat.descriptive.rank.RandomPercentile;
import org.hipparchus.stat.descriptive.summary.Sum;
import org.hipparchus.stat.descriptive.summary.SumOfLogs;
import org.hipparchus.stat.descriptive.summary.SumOfSquares;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;

/**
 * Computes summary statistics for a stream of data values added using the
 * {@link #addValue(double) addValue} method. The data values are not stored in
 * memory, so this class can be used to compute statistics for very large data
 * streams.
 * <p>
 * By default, all statistics other than percentiles are maintained.  Percentile
 * calculations use an embedded {@link RandomPercentile} which carries more memory
 * and compute overhead than the other statistics, so it is disabled by default.
 * To enable percentiles, either pass {@code true} to the constructor or use a
 * {@link StreamingStatisticsBuilder} to configure an instance with percentiles turned
 * on. Other stats can also be selectively disabled using
 * {@code StreamingStatisticsBulder}.
 * <p>
 * Note: This class is not thread-safe.
 */
public class StreamingStatistics
    implements StatisticalSummary, AggregatableStatistic<StreamingStatistics>,
               DoubleConsumer, Serializable {

    /** Serialization UID */
    private static final long serialVersionUID = 20160422L;

    /** count of values that have been added */
    private long n;

    /** SecondMoment is used to compute the mean and variance */
    private final SecondMoment secondMoment;
    /** min of values that have been added */
    private final Min minImpl;
    /** max of values that have been added */
    private final Max maxImpl;
    /** sum of values that have been added */
    private final Sum sumImpl;
    /** sum of the square of each value that has been added */
    private final SumOfSquares sumOfSquaresImpl;
    /** sumLog of values that have been added */
    private final SumOfLogs sumOfLogsImpl;
    /** mean of values that have been added */
    private final Mean meanImpl;
    /** variance of values that have been added */
    private final Variance varianceImpl;
    /** geoMean of values that have been added */
    private final GeometricMean geoMeanImpl;
    /** population variance of values that have been added */
    private final Variance populationVariance;
    /** source of percentiles */
    private final RandomPercentile randomPercentile;

    /** whether or not moment stats (sum, mean, variance) are maintained */
    private final boolean computeMoments;
    /** whether or not sum of squares and quadratic mean are maintained */
    private final boolean computeSumOfSquares;
    /** whether or not sum of logs and geometric mean are maintained */
    private final boolean computeSumOfLogs;
    /** whether or not min and max are maintained */
    private final boolean computeExtrema;

    /**
     * Construct a new StreamingStatistics instance, maintaining all statistics
     * other than percentiles.
     */
    public StreamingStatistics() {
       this(Double.NaN, null);
    }

    /**
     * Construct a new StreamingStatistics instance, maintaining all statistics
     * other than percentiles and with/without percentiles per the arguments.
     *
     * @param epsilon bound on quantile estimation error (see {@link RandomGenerator})
     * @param randomGenerator PRNG used in sampling and merge operations (null if percentiles should not be computed)
     * @since 2.3
     */
    public StreamingStatistics(final double epsilon, final RandomGenerator randomGenerator) {
       this(true, true, true, true, epsilon, randomGenerator);
    }

    /**
     * Private constructor used by {@link StreamingStatisticsBuilder}.
     *
     * @param computeMoments whether or not moment stats (mean, sum, variance) are maintained
     * @param computeSumOfLogs whether or not sum of logs and geometric mean are maintained
     * @param computeSumOfSquares whether or not sum of squares and quadratic mean are maintained
     * @param computeExtrema whether or not min and max are maintained
     * @param epsilon bound on quantile estimation error (see {@link RandomGenerator})
     * @param randomGenerator PRNG used in sampling and merge operations (null if percentiles should not be computed)
     * @since 2.3
     */
    private StreamingStatistics(final boolean computeMoments,
                                final boolean computeSumOfLogs, final boolean computeSumOfSquares,
                                final boolean computeExtrema,
                                final double epsilon, final RandomGenerator randomGenerator) {
        this.computeMoments = computeMoments;
        this.computeSumOfLogs = computeSumOfLogs;
        this.computeSumOfSquares = computeSumOfSquares;
        this.computeExtrema = computeExtrema;

        this.secondMoment = computeMoments ? new SecondMoment() : null;
        this.maxImpl = computeExtrema ? new Max() : null;
        this.minImpl = computeExtrema ? new Min() : null;
        this.sumImpl = computeMoments ? new Sum() : null;
        this.sumOfSquaresImpl = computeSumOfSquares ? new SumOfSquares() : null;
        this.sumOfLogsImpl = computeSumOfLogs ? new SumOfLogs() : null;
        this.meanImpl = computeMoments ? new Mean(this.secondMoment) : null;
        this.varianceImpl = computeMoments ?  new Variance(this.secondMoment) : null;
        this.geoMeanImpl = computeSumOfLogs ? new GeometricMean(this.sumOfLogsImpl) : null;
        this.populationVariance = computeMoments ? new Variance(false, this.secondMoment) : null;
        this.randomPercentile = randomGenerator == null ? null : new RandomPercentile(epsilon, randomGenerator);
    }

    /**
     * A copy constructor. Creates a deep-copy of the {@code original}.
     *
     * @param original the {@code StreamingStatistics} instance to copy
     * @throws NullArgumentException if original is null
     */
    StreamingStatistics(StreamingStatistics original) throws NullArgumentException {
        MathUtils.checkNotNull(original);

        this.n                = original.n;
        this.secondMoment     = original.computeMoments ? original.secondMoment.copy() : null;
        this.maxImpl          = original.computeExtrema ? original.maxImpl.copy() : null;
        this.minImpl          = original.computeExtrema ? original.minImpl.copy() : null;
        this.sumImpl          = original.computeMoments ? original.sumImpl.copy() : null;
        this.sumOfLogsImpl    = original.computeSumOfLogs ? original.sumOfLogsImpl.copy() : null;
        this.sumOfSquaresImpl = original.computeSumOfSquares ? original.sumOfSquaresImpl.copy() : null;

        // Keep statistics with embedded moments in synch
        this.meanImpl     = original.computeMoments ? new Mean(this.secondMoment) : null;
        this.varianceImpl = original.computeMoments ? new Variance(this.secondMoment) : null;
        this.geoMeanImpl  = original.computeSumOfLogs ? new GeometricMean(this.sumOfLogsImpl) : null;
        this.populationVariance = original.computeMoments ? new Variance(false, this.secondMoment) : null;
        this.randomPercentile = original.randomPercentile != null ? original.randomPercentile.copy() : null;

        this.computeMoments = original.computeMoments;
        this.computeSumOfLogs = original.computeSumOfLogs;
        this.computeSumOfSquares = original.computeSumOfSquares;
        this.computeExtrema = original.computeExtrema;
    }

    /**
     * Returns a copy of this StreamingStatistics instance with the same internal state.
     *
     * @return a copy of this
     */
    public StreamingStatistics copy() {
        return new StreamingStatistics(this);
    }

    /**
     * Return a {@link StatisticalSummaryValues} instance reporting current
     * statistics.
     * @return Current values of statistics
     */
    public StatisticalSummary getSummary() {
        return new StatisticalSummaryValues(getMean(), getVariance(), getN(),
                                            getMax(), getMin(), getSum());
    }

    /**
     * Add a value to the data
     * @param value the value to add
     */
    public void addValue(double value) {
        if (computeMoments) {
            secondMoment.increment(value);
            sumImpl.increment(value);
        }
        if (computeExtrema) {
            minImpl.increment(value);
            maxImpl.increment(value);
        }
        if (computeSumOfSquares) {
            sumOfSquaresImpl.increment(value);
        }
        if (computeSumOfLogs) {
            sumOfLogsImpl.increment(value);
        }
        if (randomPercentile != null) {
            randomPercentile.increment(value);
        }
        n++;
    }

    /** {@inheritDoc} */
    @Override
    public void accept(double value) {
        addValue(value);
    }

    /**
     * Resets all statistics and storage.
     */
    public void clear() {
        this.n = 0;
        if (computeExtrema) {
            minImpl.clear();
            maxImpl.clear();
        }
        if (computeMoments) {
            sumImpl.clear();
            secondMoment.clear();
        }
        if (computeSumOfLogs) {
            sumOfLogsImpl.clear();
        }
        if (computeSumOfSquares) {
            sumOfSquaresImpl.clear();
        }
        if (randomPercentile != null) {
            randomPercentile.clear();
        }
    }

    /** {@inheritDoc} */
    @Override
    public long getN() {
        return n;
    }

    /** {@inheritDoc} */
    @Override
    public double getMax() {
        return computeExtrema ? maxImpl.getResult() : Double.NaN;
    }

    /** {@inheritDoc} */
    @Override
    public double getMin() {
        return computeExtrema ? minImpl.getResult() : Double.NaN;
    }

    /** {@inheritDoc} */
    @Override
    public double getSum() {
        return computeMoments ? sumImpl.getResult() : Double.NaN;
    }

    /**
     * Returns the sum of the squares of the values that have been added.
     * <p>
     * Double.NaN is returned if no values have been added.
     *
     * @return The sum of squares
     */
    public double getSumOfSquares() {
        return computeSumOfSquares ? sumOfSquaresImpl.getResult() : Double.NaN;
    }

    /** {@inheritDoc} */
    @Override
    public double getMean() {
        return computeMoments ? meanImpl.getResult() : Double.NaN;
    }

    /** {@inheritDoc} */
    @Override
    public double getVariance() {
        return computeMoments ? varianceImpl.getResult() : Double.NaN;
    }

    /**
     * Returns the <a href="http://en.wikibooks.org/wiki/Statistics/Summary/Variance">
     * population variance</a> of the values that have been added.
     * <p>
     * Double.NaN is returned if no values have been added.
     *
     * @return the population variance
     */
    public double getPopulationVariance() {
        return computeMoments ? populationVariance.getResult() : Double.NaN;
    }

    /**
     * Returns the geometric mean of the values that have been added.
     * <p>
     * Double.NaN is returned if no values have been added.
     *
     * @return the geometric mean
     */
    public double getGeometricMean() {
        return computeSumOfLogs ? geoMeanImpl.getResult() : Double.NaN;
    }

    /**
     * Returns the sum of the logs of the values that have been added.
     * <p>
     * Double.NaN is returned if no values have been added.
     *
     * @return the sum of logs
     */
    public double getSumOfLogs() {
        return computeSumOfLogs ? sumOfLogsImpl.getResult() : Double.NaN;
    }

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
    public double getSecondMoment() {
        return computeMoments ? secondMoment.getResult() : Double.NaN;
    }

    /**
     * Returns the quadratic mean, a.k.a.
     * <a href="http://mathworld.wolfram.com/Root-Mean-Square.html">
     * root-mean-square</a> of the available values
     *
     * @return The quadratic mean or {@code Double.NaN} if no values
     * have been added.
     */
    public double getQuadraticMean() {
        if (computeSumOfSquares) {
            long size = getN();
            return size > 0 ? FastMath.sqrt(getSumOfSquares() / size) : Double.NaN;
        } else {
            return Double.NaN;
        }
    }

    /**
     * Returns the standard deviation of the values that have been added.
     * <p>
     * Double.NaN is returned if no values have been added.
     *
     * @return the standard deviation
     */
    @Override
    public double getStandardDeviation() {
        long size = getN();
        if (computeMoments) {
            if (size > 0) {
                return size > 1 ? FastMath.sqrt(getVariance()) : 0.0;
            } else {
                return Double.NaN;
            }
        } else {
            return Double.NaN;
        }
    }

    /**
     * Returns an estimate of the median of the values that have been entered.
     * See {@link RandomPercentile} for a description of the algorithm used for large
     * data streams.
     *
     * @return the median
     */
    public double getMedian() {
        return randomPercentile != null ? randomPercentile.getResult(50d) : Double.NaN;
    }

    /**
     * Returns an estimate of the given percentile of the values that have been entered.
     * See {@link RandomPercentile} for a description of the algorithm used for large
     * data streams.
     *
     * @param percentile the desired percentile (must be between 0 and 100)
     * @return estimated percentile
     */
    public double getPercentile(double percentile) {
        return randomPercentile == null ? Double.NaN : randomPercentile.getResult(percentile);
    }

    /**
     * {@inheritDoc}
     * Statistics are aggregated only when both this and other are maintaining them.  For example,
     * if this.computeMoments is false, but other.computeMoments is true, the moment data in other
     * will be lost.
     */
    @Override
    public void aggregate(StreamingStatistics other) {
        MathUtils.checkNotNull(other);

        if (other.n > 0) {
            this.n += other.n;
            if (computeMoments && other.computeMoments) {
                this.secondMoment.aggregate(other.secondMoment);
                this.sumImpl.aggregate(other.sumImpl);
            }
            if (computeExtrema && other.computeExtrema) {
                this.minImpl.aggregate(other.minImpl);
                this.maxImpl.aggregate(other.maxImpl);
            }
            if (computeSumOfLogs && other.computeSumOfLogs) {
                this.sumOfLogsImpl.aggregate(other.sumOfLogsImpl);
            }
            if (computeSumOfSquares && other.computeSumOfSquares) {
                this.sumOfSquaresImpl.aggregate(other.sumOfSquaresImpl);
            }
            if (randomPercentile != null && other.randomPercentile != null) {
                this.randomPercentile.aggregate(other.randomPercentile);
            }
        }
    }

    /**
     * Generates a text report displaying summary statistics from values that
     * have been added.
     *
     * @return String with line feeds displaying statistics
     */
    @Override
    public String toString() {
        StringBuilder outBuffer = new StringBuilder(200); // the size is just a wild guess
        String endl = "\n";
        outBuffer.append("StreamingStatistics:").append(endl).
                  append("n: ").append(getN()).append(endl).
                  append("min: ").append(getMin()).append(endl).
                  append("max: ").append(getMax()).append(endl).
                  append("sum: ").append(getSum()).append(endl).
                  append("mean: ").append(getMean()).append(endl).
                  append("variance: ").append(getVariance()).append(endl).
                  append("population variance: ").append(getPopulationVariance()).append(endl).
                  append("standard deviation: ").append(getStandardDeviation()).append(endl).
                  append("geometric mean: ").append(getGeometricMean()).append(endl).
                  append("second moment: ").append(getSecondMoment()).append(endl).
                  append("sum of squares: ").append(getSumOfSquares()).append(endl).
                  append("sum of logs: ").append(getSumOfLogs()).append(endl);
        return outBuffer.toString();
    }

    /**
     * Returns true iff <code>object</code> is a <code>StreamingStatistics</code>
     * instance and all statistics have the same values as this.
     *
     * @param object the object to test equality against.
     * @return true if object equals this
     */
    @Override
    public boolean equals(Object object) {
        if (object == this) {
            return true;
        }
        if (!(object instanceof StreamingStatistics)) {
            return false;
        }
        StreamingStatistics other = (StreamingStatistics)object;
        return other.getN() == getN()                                                     &&
               Precision.equalsIncludingNaN(other.getMax(),           getMax())           &&
               Precision.equalsIncludingNaN(other.getMin(),           getMin())           &&
               Precision.equalsIncludingNaN(other.getSum(),           getSum())           &&
               Precision.equalsIncludingNaN(other.getGeometricMean(), getGeometricMean()) &&
               Precision.equalsIncludingNaN(other.getMean(),          getMean())          &&
               Precision.equalsIncludingNaN(other.getSumOfSquares(),  getSumOfSquares())  &&
               Precision.equalsIncludingNaN(other.getSumOfLogs(),     getSumOfLogs())     &&
               Precision.equalsIncludingNaN(other.getVariance(),      getVariance())      &&
               Precision.equalsIncludingNaN(other.getMedian(),        getMedian());
    }

    /**
     * Returns hash code based on values of statistics.
     * @return hash code
     */
    @Override
    public int hashCode() {
        int result = 31 + MathUtils.hash(getN());
        result = result * 31 + MathUtils.hash(getMax());
        result = result * 31 + MathUtils.hash(getMin());
        result = result * 31 + MathUtils.hash(getSum());
        result = result * 31 + MathUtils.hash(getGeometricMean());
        result = result * 31 + MathUtils.hash(getMean());
        result = result * 31 + MathUtils.hash(getSumOfSquares());
        result = result * 31 + MathUtils.hash(getSumOfLogs());
        result = result * 31 + MathUtils.hash(getVariance());
        result = result * 31 + MathUtils.hash(getMedian());
        return result;
    }

    /**
     * Returns a {@link StreamingStatisticsBuilder} to source configured
     * {@code StreamingStatistics} instances.
     *
     * @return a StreamingStatisticsBuilder instance
     */
    public static StreamingStatisticsBuilder builder() {
        return new StreamingStatisticsBuilder();
    }

    /**
     * Builder for StreamingStatistics instances.
     */
    public static class StreamingStatisticsBuilder {
        /** whether or not moment statistics are maintained by instances created by this factory */
        private boolean computeMoments;
        /** whether or not sum of squares and quadratic mean are maintained by instances created by this factory */
        private boolean computeSumOfSquares;
        /** whether or not sum of logs and geometric mean are maintained by instances created by this factory */
        private boolean computeSumOfLogs;
        /** whether or not min and max are maintained by instances created by this factory */
        private boolean computeExtrema;
        /** bound on quantile estimation error for percentiles.
         * @since 2.3
         */
        private double epsilon;
        /** PRNG used in sampling and merge operations.
         * @since 2.3
         */
        private RandomGenerator randomGenerator;

        /** Simple constructor.
         */
        public StreamingStatisticsBuilder() {
            computeMoments      = true;
            computeSumOfSquares = true;
            computeSumOfLogs    = true;
            computeExtrema      = true;
            percentiles(Double.NaN, null);
        }

        /**
         * Sets the computeMoments setting of the factory
         *
         * @param arg whether or not instances created using {@link #build()} will
         * maintain moment statistics
         * @return a factory with the given computeMoments property set
         */
        public StreamingStatisticsBuilder moments(boolean arg) {
            this.computeMoments = arg;
            return this;
        }

        /**
         * Sets the computeSumOfLogs setting of the factory
         *
         * @param arg whether or not instances created using {@link #build()} will
         * maintain log sums
         * @return a factory with the given computeSumOfLogs property set
         */
        public StreamingStatisticsBuilder sumOfLogs(boolean arg) {
            this.computeSumOfLogs = arg;
            return this;
        }

        /**
         * Sets the computeSumOfSquares setting of the factory.
         *
         * @param arg whether or not instances created using {@link #build()} will
         * maintain sums of squares
         * @return a factory with the given computeSumOfSquares property set
         */
        public StreamingStatisticsBuilder sumOfSquares(boolean arg) {
            this.computeSumOfSquares = arg;
            return this;
        }

        /**
         * Sets the computePercentiles setting of the factory.
         * @param epsilonBound bound on quantile estimation error (see {@link RandomGenerator})
         * @param generator PRNG used in sampling and merge operations
         * @return a factory with the given computePercentiles property set
         * @since 2.3
         */
        public StreamingStatisticsBuilder percentiles(final double epsilonBound, final RandomGenerator generator) {
            this.epsilon         = epsilonBound;
            this.randomGenerator = generator;
            return this;
        }

        /**
         * Sets the computeExtrema setting of the factory.
         *
         * @param arg whether or not instances created using {@link #build()} will
         * compute min and max
         * @return a factory with the given computeExtrema property set
         */
        public StreamingStatisticsBuilder extrema(boolean arg) {
            this.computeExtrema = arg;
            return this;
        }

        /**
         * Builds a StreamingStatistics instance with currently defined properties.
         *
         * @return newly configured StreamingStatistics instance
         */
        public StreamingStatistics build() {
            return new StreamingStatistics(computeMoments,
                                           computeSumOfLogs, computeSumOfSquares,
                                           computeExtrema,
                                           epsilon, randomGenerator);
        }
    }
}
