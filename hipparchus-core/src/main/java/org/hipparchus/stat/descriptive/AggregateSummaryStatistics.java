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

import java.util.Collection;

import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.MathUtils;

/**
 * An aggregator for {@code SummaryStatistics} from several data sets or
 * data set partitions.  In its simplest usage mode, the client creates an
 * instance via the zero-argument constructor, then uses
 * {@link #createContributingStatistics()} to obtain a {@code SummaryStatistics}
 * for each individual data set / partition.  The per-set statistics objects
 * are used as normal, and at any time the aggregate statistics for all the
 * contributors can be obtained from this object.
 * <p>
 * Clients with specialized requirements can use alternative constructors to
 * control the statistics implementations and initial values used by the
 * contributing and the internal aggregate {@code SummaryStatistics} objects.
 * <p>
 * A static {@link #aggregate(Collection)} method is also included that computes
 * aggregate statistics directly from a Collection of SummaryStatistics instances.
 * <p>
 * When {@link #createContributingStatistics()} is used to create SummaryStatistics
 * instances to be aggregated concurrently, the created instances'
 * {@link SummaryStatistics#addValue(double)} methods must synchronize on the
 * aggregating instance maintained by this class.  In multi-threaded environments,
 * if the functionality provided by {@link #aggregate(Collection)} is adequate,
 * that method should be used to avoid unnecessary computation and synchronization
 * delays.
 * <p>
 * <b>Note:</b> {@link AggregateSummaryStatistics} and its contributing statistics
 * are <b>not</b> serializable.
 */
public class AggregateSummaryStatistics
    extends ForwardingSummaryStatistics {

    /**
     * A SummaryStatistics serving as a prototype for creating SummaryStatistics
     * contributing to this aggregate.
     */
    private final SummaryStatistics statisticsPrototype;

    /**
     * Initializes a new AggregateSummaryStatistics with default statistics
     * implementations.
     */
    public AggregateSummaryStatistics() {
        this(SummaryStatistics.create());
    }

    /**
     * Initializes a new AggregateSummaryStatistics with the specified statistics
     * object as a prototype for contributing statistics and for the internal
     * aggregate statistics.
     * <p>
     * This provides for customized statistics implementations to be used by
     * contributing and aggregate statistics.
     *
     * @param prototypeStatistics a {@code SummaryStatistics} serving as a
     * prototype both for the internal aggregate statistics and for
     * contributing statistics obtained via the {@code createContributingStatistics()}
     * method.  Being a prototype means that other objects are initialized by copying
     * this object's state. Any statistic values in the prototype are propagated to
     * contributing statistics objects and (once) into these aggregate statistics.
     * @throws NullArgumentException if prototypeStatistics is null
     * @see #createContributingStatistics()
     */
    public AggregateSummaryStatistics(SummaryStatistics prototypeStatistics)
        throws NullArgumentException {
        this(prototypeStatistics,
             prototypeStatistics == null ? null : prototypeStatistics.copy());
    }

    /**
     * Initializes a new AggregateSummaryStatistics with the specified statistics
     * object as a prototype for contributing statistics and for the internal
     * aggregate statistics.
     * <p>
     * This provides for different statistics implementations to be used by
     * contributing and aggregate statistics and for an initial state to be
     * supplied for the aggregate statistics.
     *
     * @param prototypeStatistics a {@code SummaryStatistics} serving as a
     * prototype for contributing statistics obtained via the
     * {@code createContributingStatistics()} method.  Being a prototype means
     * that other objects are initialized by copying this object's state.
     * Any statistic values in the prototype are propagated to contributing
     * statistics objects and (once) into these aggregate statistics.
     * @param initialStatistics a {@code SummaryStatistics} to serve as the
     * internal aggregate statistics object.
     * @throws NullArgumentException if prototypeStatistics or initialStatistics is null
     * @see #createContributingStatistics()
     */
    public AggregateSummaryStatistics(SummaryStatistics prototypeStatistics,
                                      SummaryStatistics initialStatistics) {
        super(SummaryStatistics.synchronizedSummaryStatistics(initialStatistics));
        MathUtils.checkNotNull(prototypeStatistics);
        this.statisticsPrototype = prototypeStatistics;
    }

    /** {@inheritDoc} */
    @Override
    public AggregateSummaryStatistics copy() {
        return new AggregateSummaryStatistics(statisticsPrototype.copy(),
                                              delegate().copy());
    }

    /**
     * Creates and returns a {@code SummaryStatistics} whose data will be
     * aggregated with those of this {@code AggregateSummaryStatistics}.
     *
     * @return a {@code SummaryStatistics} whose data will be aggregated with
     * those of this {@code AggregateSummaryStatistics}.  The initial state
     * is a copy of the configured prototype statistics.
     */
    public SummaryStatistics createContributingStatistics() {
        SummaryStatistics contributingStatistics
            = new AggregatingSummaryStatistics(statisticsPrototype.copy(),
                                               delegate());
        return contributingStatistics;
    }

    /**
     * Computes aggregate summary statistics.
     * <p>
     * This method can be used to combine statistics computed over partitions or
     * subsamples - i.e., the StatisticalSummaryValues returned should contain
     * the same values that would have been obtained by computing a single
     * StatisticalSummary over the combined dataset.
     *
     * @param statistics collection of SummaryStatistics to aggregate
     * @return summary statistics for the combined dataset
     * @throws NullArgumentException if the input is null
     */
    public static StatisticalSummary aggregate(Collection<? extends StatisticalSummary> statistics) {
        MathUtils.checkNotNull(statistics);

        long n = 0;
        double min = Double.NaN;
        double max = Double.NaN;
        double sum = Double.NaN;
        double mean = Double.NaN;
        double m2 = Double.NaN;

        for (StatisticalSummary current : statistics) {
            if (current.getN() == 0) {
                continue;
            }

            if (n == 0) {
                n = current.getN();
                min = current.getMin();
                sum = current.getSum();
                max = current.getMax();
                m2 = current.getVariance() * (n - 1);
                mean = current.getMean();
            } else {
                if (current.getMin() < min) {
                    min = current.getMin();
                }
                if (current.getMax() > max) {
                    max = current.getMax();
                }

                sum += current.getSum();
                final double oldN = n;
                final double curN = current.getN();
                n += curN;
                final double meanDiff = current.getMean() - mean;
                mean = sum / n;
                final double curM2 = current.getVariance() * (curN - 1d);
                m2 = m2 + curM2 + meanDiff * meanDiff * oldN * curN / n;
            }
        }

        final double variance = n == 0 ? Double.NaN :
                                n == 1 ? 0d         :
                                         m2 / (n - 1);

        return new StatisticalSummaryValues(mean, variance, n, max, min, sum);
    }

    @Override
    public String toString() {
        return "Aggregate" + delegate().toString();
    }

    /**
     * A SummaryStatistics that also forwards all values added to it to a second
     * {@code SummaryStatistics} for aggregation.
     */
    private static class AggregatingSummaryStatistics extends ForwardingSummaryStatistics {

        /**
         * An additional SummaryStatistics into which values added to these
         * statistics (and possibly others) are aggregated
         */
        private final SummaryStatistics aggregateStatistics;

        /**
         * Initializes a new AggregatingSummaryStatistics with the specified
         * aggregate statistics object
         *
         * @param prototypeStatistics the prototype statistics
         * @param aggregateStatistics a {@code SummaryStatistics} into which
         *      values added to this statistics object should be aggregated
         */
        AggregatingSummaryStatistics(SummaryStatistics prototypeStatistics,
                                     SummaryStatistics aggregateStatistics) {
            super(prototypeStatistics);
            this.aggregateStatistics = aggregateStatistics;
        }

        @Override
        public SummaryStatistics copy() {
            return new AggregateSummaryStatistics(delegate().copy(),
                                                  aggregateStatistics);
        }

        /**
         * {@inheritDoc}.
         * <p>
         * This version adds the provided value to the configured
         * aggregate after adding it to these statistics.
         *
         * @see SummaryStatistics#addValue(double)
         */
        @Override
        public void addValue(double value) {
            super.addValue(value);
            // no need to additionally sync the access
            // as the aggregate statistics is already synchronized.
            aggregateStatistics.addValue(value);
        }

    }
}
