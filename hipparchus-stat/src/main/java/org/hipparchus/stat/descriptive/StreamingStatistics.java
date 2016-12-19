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

import java.io.Serializable;
import java.util.function.DoubleConsumer;

import org.hipparchus.exception.NullArgumentException;
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
 * Note: This class is not thread-safe.
 */
public class StreamingStatistics
    implements StatisticalSummary, AggregatableStatistic<StreamingStatistics>,
               DoubleConsumer, Serializable {

    /** Serialization UID */
    private static final long serialVersionUID = 20160422L;

    /** count of values that have been added */
    private long n = 0;

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

    /**
     * Construct a new StreamingStatistics instance.
     */
    public StreamingStatistics() {
        this.secondMoment     = new SecondMoment();
        this.maxImpl          = new Max();
        this.minImpl          = new Min();
        this.sumImpl          = new Sum();
        this.sumOfSquaresImpl = new SumOfSquares();
        this.sumOfLogsImpl    = new SumOfLogs();
        this.meanImpl         = new Mean(this.secondMoment);
        this.varianceImpl     = new Variance(this.secondMoment);
        this.geoMeanImpl      = new GeometricMean(this.sumOfLogsImpl);

        // the population variance can not be overridden
        // it will always use the second moment.
        this.populationVariance = new Variance(false, this.secondMoment);
        this.randomPercentile = new RandomPercentile();
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
        this.secondMoment     = original.secondMoment.copy();
        this.maxImpl          = original.maxImpl.copy();
        this.minImpl          = original.minImpl.copy();
        this.sumImpl          = original.sumImpl.copy();
        this.sumOfLogsImpl    = original.sumOfLogsImpl.copy();
        this.sumOfSquaresImpl = original.sumOfSquaresImpl.copy();

        // Keep default statistics with embedded moments in synch
        this.meanImpl     = new Mean(this.secondMoment);
        this.varianceImpl = new Variance(this.secondMoment);
        this.geoMeanImpl  = new GeometricMean(this.sumOfLogsImpl);

        this.populationVariance = new Variance(false, this.secondMoment);
        this.randomPercentile = original.randomPercentile.copy();
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
        secondMoment.increment(value);
        minImpl.increment(value);
        maxImpl.increment(value);
        sumImpl.increment(value);
        sumOfSquaresImpl.increment(value);
        sumOfLogsImpl.increment(value);
        randomPercentile.increment(value);

        // Do not update mean/variance/geoMean
        // as they use external moments.

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
        minImpl.clear();
        maxImpl.clear();
        sumImpl.clear();
        sumOfLogsImpl.clear();
        sumOfSquaresImpl.clear();
        secondMoment.clear();
        randomPercentile.clear();

        // No need to clear mean/variance/geoMean
        // as they use external moments.
    }

    /** {@inheritDoc} */
    @Override
    public long getN() {
        return n;
    }

    /** {@inheritDoc} */
    @Override
    public double getMax() {
        return maxImpl.getResult();
    }

    /** {@inheritDoc} */
    @Override
    public double getMin() {
        return minImpl.getResult();
    }

    /** {@inheritDoc} */
    @Override
    public double getSum() {
        return sumImpl.getResult();
    }

    /**
     * Returns the sum of the squares of the values that have been added.
     * <p>
     * Double.NaN is returned if no values have been added.
     *
     * @return The sum of squares
     */
    public double getSumOfSquares() {
        return sumOfSquaresImpl.getResult();
    }

    /** {@inheritDoc} */
    @Override
    public double getMean() {
        return meanImpl.getResult();
    }

    /** {@inheritDoc} */
    @Override
    public double getVariance() {
        return varianceImpl.getResult();
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
        return populationVariance.getResult();
    }

    /**
     * Returns the geometric mean of the values that have been added.
     * <p>
     * Double.NaN is returned if no values have been added.
     *
     * @return the geometric mean
     */
    public double getGeometricMean() {
        return geoMeanImpl.getResult();
    }

    /**
     * Returns the sum of the logs of the values that have been added.
     * <p>
     * Double.NaN is returned if no values have been added.
     *
     * @return the sum of logs
     */
    public double getSumOfLogs() {
        return sumOfLogsImpl.getResult();
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
        return secondMoment.getResult();
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
        long size = getN();
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
    public double getStandardDeviation() {
        long size = getN();
        if (size > 0) {
            return size > 1 ? FastMath.sqrt(getVariance()) : 0.0;
        } else {
            return Double.NaN;
        }
    }

    public double getMedian() {
        return randomPercentile != null ? randomPercentile.getResult(50d) : Double.NaN;
    }

    public double getPercentile(double percentile) {
        return randomPercentile != null ? randomPercentile.getResult(percentile) : Double.NaN;
    }

    /** {@inheritDoc} */
    @Override
    public void aggregate(StreamingStatistics other) {
        MathUtils.checkNotNull(other);

        if (other.n > 0) {
            this.n += other.n;
            this.secondMoment.aggregate(other.secondMoment);
            this.minImpl.aggregate(other.minImpl);
            this.maxImpl.aggregate(other.maxImpl);
            this.sumImpl.aggregate(other.sumImpl);
            this.sumOfLogsImpl.aggregate(other.sumOfLogsImpl);
            this.sumOfSquaresImpl.aggregate(other.sumOfSquaresImpl);
            this.randomPercentile.aggregate(other.randomPercentile);
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
        StringBuilder outBuffer = new StringBuilder();
        String endl = "\n";
        outBuffer.append("StreamingStatistics:").append(endl);
        outBuffer.append("n: ").append(getN()).append(endl);
        outBuffer.append("min: ").append(getMin()).append(endl);
        outBuffer.append("max: ").append(getMax()).append(endl);
        outBuffer.append("sum: ").append(getSum()).append(endl);
        outBuffer.append("mean: ").append(getMean()).append(endl);
        outBuffer.append("variance: ").append(getVariance()).append(endl);
        outBuffer.append("population variance: ").append(getPopulationVariance()).append(endl);
        outBuffer.append("standard deviation: ").append(getStandardDeviation()).append(endl);
        outBuffer.append("geometric mean: ").append(getGeometricMean()).append(endl);
        outBuffer.append("second moment: ").append(getSecondMoment()).append(endl);
        outBuffer.append("sum of squares: ").append(getSumOfSquares()).append(endl);
        outBuffer.append("sum of logs: ").append(getSumOfLogs()).append(endl);
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
        if (object instanceof StreamingStatistics == false) {
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

}
