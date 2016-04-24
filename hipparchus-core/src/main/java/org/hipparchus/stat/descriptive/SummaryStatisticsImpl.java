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

import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.stat.descriptive.moment.GeometricMean;
import org.hipparchus.stat.descriptive.moment.Mean;
import org.hipparchus.stat.descriptive.moment.SecondMoment;
import org.hipparchus.stat.descriptive.moment.Variance;
import org.hipparchus.stat.descriptive.summary.SumOfLogs;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;

/**
 * Computes summary statistics for a stream of data values added using the
 * {@link #addValue(double) addValue} method. The data values are not stored in
 * memory, so this class can be used to compute statistics for very large data
 * streams.
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
 * Note: This class is not thread-safe.
 */
class SummaryStatisticsImpl implements SummaryStatistics, Serializable {

    /** Serialization UID */
    private static final long serialVersionUID = 20160422L;

    /** count of values that have been added */
    private long n = 0;

    /** SecondMoment is used to compute the mean and variance */
    private final SecondMoment secondMoment;
    /** min of values that have been added */
    private final StorelessUnivariateStatistic minImpl;
    /** max of values that have been added */
    private final StorelessUnivariateStatistic maxImpl;
    /** sum of values that have been added */
    private final StorelessUnivariateStatistic sumImpl;
    /** sum of the square of each value that has been added */
    private final StorelessUnivariateStatistic sumOfSquaresImpl;
    /** sumLog of values that have been added */
    private final StorelessUnivariateStatistic sumOfLogsImpl;
    /** mean of values that have been added */
    private final StorelessUnivariateStatistic meanImpl;
    /** variance of values that have been added */
    private final StorelessUnivariateStatistic varianceImpl;
    /** geoMean of values that have been added */
    private final StorelessUnivariateStatistic geoMeanImpl;
    /** population variance of values that have been added */
    private final StorelessUnivariateStatistic populationVariance;

    /** Indicates if the mean impl uses an external moment */
    private final boolean meanUsesExternalMoment;
    /** Indicates if the variance impl uses an external moment */
    private final boolean varianceUsesExternalMoment;
    /** Indicates if the geo mean impl uses an external sum of logs */
    private final boolean geoMeanUsesExternalSumOfLogs;

    /**
     * A copy constructor. Creates a deep-copy of the {@code original}.
     *
     * @param original the {@code SummaryStatistics} instance to copy
     * @throws NullArgumentException if original is null
     */
    public SummaryStatisticsImpl(SummaryStatisticsImpl original) throws NullArgumentException {
        MathUtils.checkNotNull(original);

        this.n                = original.n;
        this.secondMoment     = original.secondMoment.copy();
        this.maxImpl          = original.maxImpl.copy();
        this.minImpl          = original.minImpl.copy();
        this.sumImpl          = original.sumImpl.copy();
        this.sumOfLogsImpl    = original.sumOfLogsImpl.copy();
        this.sumOfSquaresImpl = original.sumOfSquaresImpl.copy();

        // Keep default statistics with embedded moments in synch
        if (original.meanUsesExternalMoment) {
            this.meanImpl = new Mean(this.secondMoment);
        } else {
            this.meanImpl = original.meanImpl.copy();
        }
        this.meanUsesExternalMoment = original.meanUsesExternalMoment;

        if (original.varianceUsesExternalMoment) {
            this.varianceImpl = new Variance(this.secondMoment);
        } else {
            this.varianceImpl = original.varianceImpl.copy();
        }
        this.varianceUsesExternalMoment = original.varianceUsesExternalMoment;

        if (original.geoMeanUsesExternalSumOfLogs &&
            this.sumOfLogsImpl instanceof SumOfLogs) {
            this.geoMeanImpl = new GeometricMean((SumOfLogs) this.sumOfLogsImpl);
        } else {
            this.geoMeanImpl = original.geoMeanImpl.copy();
        }
        this.geoMeanUsesExternalSumOfLogs = original.geoMeanUsesExternalSumOfLogs;

        this.populationVariance = new Variance(false, this.secondMoment);
    }

    /**
     * Construct a new SummaryStatistics instance based
     * on the data provided by a builder.
     *
     * @param builder the builder to use.
     */
    protected SummaryStatisticsImpl(Builder builder) {
        this.secondMoment       = builder.secondMoment;
        this.maxImpl            = builder.maxImpl;
        this.minImpl            = builder.minImpl;
        this.meanImpl           = builder.meanImpl;
        this.sumImpl            = builder.sumImpl;
        this.sumOfSquaresImpl   = builder.sumOfSquaresImpl;
        this.sumOfLogsImpl      = builder.sumOfLogsImpl;
        this.varianceImpl       = builder.varianceImpl;
        this.geoMeanImpl        = builder.geoMeanImpl;

        // the population variance can not be overridden
        // it will always use the second moment.
        this.populationVariance = new Variance(false, this.secondMoment);

        this.meanUsesExternalMoment       = builder.meanUsesExternalMoment;
        this.varianceUsesExternalMoment   = builder.varianceUsesExternalMoment;
        this.geoMeanUsesExternalSumOfLogs = builder.geometricMeanUsesExternalSumOfLogs;
    }

    /** {@inheritDoc} */
    @Override
    public SummaryStatistics copy() {
        return new SummaryStatisticsImpl(this);
    }

    /** {@inheritDoc} */
    @Override
    public StatisticalSummary getSummary() {
        return new StatisticalSummaryValues(getMean(), getVariance(), getN(),
                                            getMax(), getMin(), getSum());
    }

    /** {@inheritDoc} */
    @Override
    public void addValue(double value) {
        secondMoment.increment(value);
        minImpl.increment(value);
        maxImpl.increment(value);
        sumImpl.increment(value);
        sumOfSquaresImpl.increment(value);
        sumOfLogsImpl.increment(value);

        // update mean/variance/geoMean if they
        // do not use external moments / sumOfLogs.

        if (!meanUsesExternalMoment) {
            meanImpl.increment(value);
        }
        if (!varianceUsesExternalMoment) {
            varianceImpl.increment(value);
        }
        if (!geoMeanUsesExternalSumOfLogs) {
            geoMeanImpl.increment(value);
        }

        n++;
    }

    /** {@inheritDoc} */
    @Override
    public void clear() {
        this.n = 0;
        minImpl.clear();
        maxImpl.clear();
        sumImpl.clear();
        sumOfLogsImpl.clear();
        sumOfSquaresImpl.clear();
        secondMoment.clear();
        if (!meanUsesExternalMoment) {
            meanImpl.clear();
        }
        if (!varianceUsesExternalMoment) {
            varianceImpl.clear();
        }
        if (!geoMeanUsesExternalSumOfLogs) {
            geoMeanImpl.clear();
        }
    }

    /** {@inheritDoc} */
    @Override
    public long getN() {
        return n;
    }

    /** {@inheritDoc} */
    @Override
    public double getSum() {
        return sumImpl.getResult();
    }

    /** {@inheritDoc} */
    @Override
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

    /** {@inheritDoc} */
    @Override
    public double getPopulationVariance() {
        return populationVariance.getResult();
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
    public double getGeometricMean() {
        return geoMeanImpl.getResult();
    }

    /** {@inheritDoc} */
    @Override
    public double getSumOfLogs() {
        return sumOfLogsImpl.getResult();
    }

    /** {@inheritDoc} */
    @Override
    public double getSecondMoment() {
        return secondMoment.getResult();
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
        outBuffer.append("SummaryStatistics:").append(endl);
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
     * Returns true iff <code>object</code> is a <code>SummaryStatistics</code>
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
        if (object instanceof SummaryStatistics == false) {
            return false;
        }
        SummaryStatistics other = (SummaryStatistics)object;
        return Precision.equalsIncludingNaN(other.getN(),             getN())             &&
               Precision.equalsIncludingNaN(other.getMax(),           getMax())           &&
               Precision.equalsIncludingNaN(other.getMin(),           getMin())           &&
               Precision.equalsIncludingNaN(other.getSum(),           getSum())           &&
               Precision.equalsIncludingNaN(other.getGeometricMean(), getGeometricMean()) &&
               Precision.equalsIncludingNaN(other.getMean(),          getMean())          &&
               Precision.equalsIncludingNaN(other.getSumOfSquares(),  getSumOfSquares())  &&
               Precision.equalsIncludingNaN(other.getSumOfLogs(),     getSumOfLogs())     &&
               Precision.equalsIncludingNaN(other.getVariance(),      getVariance());
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
        return result;
    }

}
