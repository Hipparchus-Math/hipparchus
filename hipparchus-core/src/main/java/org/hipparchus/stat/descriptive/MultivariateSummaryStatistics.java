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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.util.FastMath;

/**
 * Computes summary statistics for a stream of n-tuples added using the
 * {@link #addValue(double[]) addValue} method. The data values are not stored
 * in memory, so this class can be used to compute statistics for very large
 * n-tuple streams.
 * <p>
 * To compute statistics for a stream of n-tuples, construct a
 * {@link MultivariateSummaryStatistics} instance with dimension n and then use
 * {@link #addValue(double[])} to add n-tuples. The <code>getXxx</code>
 * methods where Xxx is a statistic return an array of <code>double</code>
 * values, where for <code>i = 0,...,n-1</code> the i<sup>th</sup> array element
 * is the value of the given statistic for data range consisting of the i<sup>th</sup>
 * element of each of the input n-tuples.  For example, if <code>addValue</code> is
 * called with actual parameters {0, 1, 2}, then {3, 4, 5} and finally {6, 7, 8},
 * <code>getSum</code> will return a three-element array with values {0+3+6, 1+4+7, 2+5+8}
 * <p>
 * Note: This class is not thread-safe.
 * Use {@link #synchronizedSummaryStatistics(MultivariateSummaryStatistics)}
 * if concurrent access from multiple threads is required.
 */
public interface MultivariateSummaryStatistics
    extends StatisticalMultivariateSummary {

    /**
     * Returns a builder for a {@link MultivariateSummaryStatistics}.
     *
     * @param k the dimension of the data
     * @return a summary statistics builder.
     */
    static Builder builder(int k) {
        return new Builder(k);
    }

    /**
     * Construct a MultivariateSummaryStatistics instance for the given
     * dimension. The returned instance will compute the unbiased sample
     * covariance.
     * <p>
     * The returned instance is <b>not</b> thread-safe.
     *
     * @param k dimension of the data
     */
    static MultivariateSummaryStatistics of(int k) {
        return builder(k).withBiasCorrectedCovariance(true).build();
    }

    /**
     * Decorates another MultivariateSummaryStatistics to synchronize its behaviour
     * for a multi-threaded environment.
     * <p>
     * Methods are synchronized, then forwarded to the decorated MultivariateSummaryStatistics.
     *
     * @param statistic the summary statistic to decorate
     * @return a synchronized wrapper for the summary statistic
     */
    static MultivariateSummaryStatistics synchronizedSummaryStatistics(MultivariateSummaryStatistics statistic) {
        if (statistic instanceof SynchronizedMultivariateSummaryStatistics) {
            return statistic;
        }
        return new SynchronizedMultivariateSummaryStatistics(statistic);
    }

    /**
     * Add an n-tuple to the data
     *
     * @param value  the n-tuple to add
     * @throws MathIllegalArgumentException if the array is null or the length
     * of the array does not match the one used at construction
     */
    void addValue(double[] value) throws MathIllegalArgumentException;

    /**
     * Resets all statistics and storage.
     */
    void clear();

    /**
     * Returns an array whose i<sup>th</sup> entry is the standard deviation of the
     * i<sup>th</sup> entries of the arrays that have been added using
     * {@link #addValue(double[])}
     *
     * @return the array of component standard deviations
     */
    @Override
    default double[] getStandardDeviation() {
        final int k = getDimension();
        double[] stdDev = new double[k];
        if (getN() < 1) {
            Arrays.fill(stdDev, Double.NaN);
        } else if (getN() < 2) {
            Arrays.fill(stdDev, 0.0);
        } else {
            RealMatrix matrix = getCovariance();
            for (int i = 0; i < k; ++i) {
                stdDev[i] = FastMath.sqrt(matrix.getEntry(i, i));
            }
        }
        return stdDev;
    }

    /**
     * A mutable builder for a MultivariateSummaryStatistics.
     */
    public static class Builder {
        final int k;
        boolean   covarianceBiasCorrected = true;
        boolean   threadsafe;

        Builder(int k) {
            this.k = k;
        }

        /**
         * Indicates if the returned instance shall compute the unbiased
         * sample or biased population covariance.
         *
         * @param biasCorrected if true, the unbiased sample
         * covariance is computed, otherwise the biased population covariance
         * is computed
         * @return the builder
         */
        public Builder withBiasCorrectedCovariance(boolean biasCorrected) {
            this.covarianceBiasCorrected = biasCorrected;
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
         * Constructs a new MultivariateSummaryStatistics instance with the values
         * stored in this builder.
         *
         * @return a new MultivariateSummaryStatistics instance.
         */
        public MultivariateSummaryStatistics build() {
            final MultivariateSummaryStatistics instance = new MultivariateSummaryStatisticsImpl(this);
            return threadsafe ? synchronizedSummaryStatistics(instance) : instance;
        }
    }

}
