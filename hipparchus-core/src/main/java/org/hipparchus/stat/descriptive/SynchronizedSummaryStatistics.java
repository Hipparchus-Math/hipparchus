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

/**
 * Implementation of
 * {@link org.hipparchus.stat.descriptive.SummaryStatistics} that
 * is safe to use in a multithreaded environment.  Multiple threads can safely
 * operate on a single instance without causing runtime exceptions due to race
 * conditions.  In effect, this implementation makes modification and access
 * methods atomic operations for a single instance.  That is to say, as one
 * thread is computing a statistic from the instance, no other thread can modify
 * the instance nor compute another statistic.
 */
public class SynchronizedSummaryStatistics extends SummaryStatistics {

    /** Serialization UID */
    private static final long serialVersionUID = 20160422L;

    /**
     * Returns a builder for a {@link SynchronizedSummaryStatistics}.
     *
     * @return a descriptive statistics builder.
     */
    public static Builder builder() {
        return new Builder();
    }

    /**
     * Construct a SynchronizedSummaryStatistics instance
     */
    public SynchronizedSummaryStatistics() {
        super();
    }

    /**
     * A copy constructor. Creates a deep-copy of the {@code original}.
     *
     * @param original the {@code SynchronizedSummaryStatistics} instance to copy
     * @throws NullArgumentException if original is null
     */
    private SynchronizedSummaryStatistics(SynchronizedSummaryStatistics original)
        throws NullArgumentException {
        super(original);
    }

    /**
     * Construct a new SynchronizedSummaryStatistics instance based
     * on the data provided by a builder.
     *
     * @param builder the builder to use.
     */
    protected SynchronizedSummaryStatistics(Builder builder) {
        super(builder);
    }

    /** {@inheritDoc} */
    @Override
    public synchronized StatisticalSummary getSummary() {
        return super.getSummary();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized void addValue(double value) {
        super.addValue(value);
    }

    /** {@inheritDoc} */
    @Override
    public synchronized long getN() {
        return super.getN();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getSum() {
        return super.getSum();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getSumOfSquares() {
        return super.getSumOfSquares();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getMean() {
        return super.getMean();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getStandardDeviation() {
        return super.getStandardDeviation();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getQuadraticMean() {
        return super.getQuadraticMean();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getVariance() {
        return super.getVariance();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getPopulationVariance() {
        return super.getPopulationVariance();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getMax() {
        return super.getMax();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getMin() {
        return super.getMin();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getGeometricMean() {
        return super.getGeometricMean();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized String toString() {
        return super.toString();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized void clear() {
        super.clear();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized boolean equals(Object object) {
        return super.equals(object);
    }

    /** {@inheritDoc} */
    @Override
    public synchronized int hashCode() {
        return super.hashCode();
    }

    /**
     * Returns a copy of this SynchronizedSummaryStatistics instance with the
     * same internal state.
     *
     * @return a copy of this
     */
    @Override
    public synchronized SynchronizedSummaryStatistics copy() {
        return new SynchronizedSummaryStatistics(this);
    }

    /**
     * A mutable builder for a SynchronizedSummaryStatistics.
     */
    public static class Builder extends SummaryStatistics.Builder {

        @Override
        public SynchronizedSummaryStatistics build() {
            setupDefaultsImpls();
            return new SynchronizedSummaryStatistics(this);
        }

    }

}
