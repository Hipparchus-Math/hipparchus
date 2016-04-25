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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.util.MathUtils;

/**
 * Implementation of {@link MultivariateSummaryStatistics} that
 * is safe to use in a multithreaded environment.
 * <p>
 * Multiple threads can safely operate on a single instance without causing
 * runtime exceptions due to race conditions.  In effect, this implementation
 * makes modification and access methods atomic operations for a single instance.
 * <p>
 * That is to say, as one thread is computing a statistic from the instance,
 * no other thread can modify the instance nor compute another statistic.
 */
final class SynchronizedMultivariateSummaryStatistics
    implements MultivariateSummaryStatistics, Serializable {

    /** Serialization UID */
    private static final long serialVersionUID = 20160413L;

    /** The MultivariateSummaryStatistics instance to delegate to. */
    private final MultivariateSummaryStatistics delegate;

    /**
     * Constructor that wraps the given instance.
     *
     * @param original the {@code MultivariateSummaryStatistics} instance to wrap
     * @throws NullArgumentException if original is null
     */
    SynchronizedMultivariateSummaryStatistics(MultivariateSummaryStatistics original)
        throws NullArgumentException {
        MathUtils.checkNotNull(original);
        this.delegate = original;
    }

    /** {@inheritDoc} */
    @Override
    public synchronized void addValue(double[] value) throws MathIllegalArgumentException {
        delegate.addValue(value);
    }

    /** {@inheritDoc} */
    @Override
    public synchronized void clear() {
        delegate.clear();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized int getDimension() {
        return delegate.getDimension();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized long getN() {
        return delegate.getN();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double[] getSum() {
        return delegate.getSum();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double[] getSumSq() {
        return delegate.getSumSq();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double[] getSumLog() {
        return delegate.getSumLog();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double[] getMean() {
        return delegate.getMean();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double[] getStandardDeviation() {
        return delegate.getStandardDeviation();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized RealMatrix getCovariance() {
        return delegate.getCovariance();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double[] getMax() {
        return delegate.getMax();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double[] getMin() {
        return delegate.getMin();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double[] getGeometricMean() {
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

}
