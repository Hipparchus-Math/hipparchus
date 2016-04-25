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
import org.hipparchus.util.MathUtils;

/**
 * Implementation of {@link SummaryStatistics} that is safe to use in a
 * multi-threaded environment.
 * <p>
 * Multiple threads can safely operate on a single instance without causing
 * runtime exceptions due to race conditions.  In effect, this implementation
 * makes modification and access methods atomic operations for a single instance.
 * That is to say, as one thread is computing a statistic from the instance,
 * no other thread can modify the instance nor compute another statistic.
 */
final class SynchronizedSummaryStatistics
    implements SummaryStatistics, Serializable {

    /** Serialization UID */
    private static final long serialVersionUID = 20160422L;

    /** The SummaryStatistics instance to delegate to. */
    private final SummaryStatistics delegate;

    /**
     * Constructor that wraps the given instance.
     *
     * @param original the {@code SummaryStatistics} instance to wrap
     * @throws NullArgumentException if original is null
     */
    SynchronizedSummaryStatistics(SummaryStatistics original)
        throws NullArgumentException {
        MathUtils.checkNotNull(original);
        this.delegate = original;
    }

    /** {@inheritDoc} */
    @Override
    public synchronized StatisticalSummary getSummary() {
        return delegate.getSummary();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized void addValue(double value) {
        delegate.addValue(value);
    }

    /** {@inheritDoc} */
    @Override
    public synchronized void clear() {
        delegate.clear();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized long getN() {
        return delegate.getN();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getSum() {
        return delegate.getSum();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getSumOfSquares() {
        return delegate.getSumOfSquares();
    }

    /** {@inheritDoc} */
    @Override
    public double getSumOfLogs() {
        return delegate.getSumOfLogs();
    }

    /** {@inheritDoc} */
    @Override
    public double getSecondMoment() {
        return delegate.getSecondMoment();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getMean() {
        return delegate.getMean();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getStandardDeviation() {
        return delegate.getStandardDeviation();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getQuadraticMean() {
        return delegate.getQuadraticMean();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getVariance() {
        return delegate.getVariance();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getPopulationVariance() {
        return delegate.getPopulationVariance();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getMax() {
        return delegate.getMax();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getMin() {
        return delegate.getMin();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getGeometricMean() {
        return delegate.getGeometricMean();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized String toString() {
        return delegate.toString();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized boolean equals(Object object) {
        return delegate.equals(object);
    }

    /** {@inheritDoc} */
    @Override
    public synchronized int hashCode() {
        return delegate.hashCode();
    }

    /**
     * Returns a copy of this SynchronizedSummaryStatistics instance with the
     * same internal state.
     *
     * @return a copy of this
     */
    @Override
    public synchronized SummaryStatistics copy() {
        return new SynchronizedSummaryStatistics(this.delegate.copy());
    }

}
