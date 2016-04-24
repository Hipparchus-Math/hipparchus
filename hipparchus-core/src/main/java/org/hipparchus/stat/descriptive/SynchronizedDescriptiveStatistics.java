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
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.MathUtils;

/**
 * Implementation of
 * {@link org.hipparchus.stat.descriptive.DescriptiveStatistics} that
 * is safe to use in a multithreaded environment.
 * <p>
 * Multiple threads can safely operate on a single instance without causing
 * runtime exceptions due to race conditions.  In effect, this implementation
 * makes modification and access methods atomic operations for a single instance.
 * That is to say, as one thread is computing a statistic from the instance,
 * no other thread can modify the instance nor compute another statistic.
 */
final class SynchronizedDescriptiveStatistics
    implements DescriptiveStatistics, Serializable {

    /** Serialization UID */
    private static final long serialVersionUID = 20160424L;

    /** The DescriptiveStatistics instance to delegate to. */
    private final DescriptiveStatistics delegate;

    /**
     * Constructor that wraps the given instance.
     *
     * @param original the {@code DescriptiveStatistics} instance to wrap
     * @throws NullArgumentException if original is null
     */
    SynchronizedDescriptiveStatistics(DescriptiveStatistics original)
        throws NullArgumentException {
        MathUtils.checkNotNull(original);
        this.delegate = original;
    }

    /** {@inheritDoc} */
    @Override
    public synchronized SynchronizedDescriptiveStatistics copy() {
        return new SynchronizedDescriptiveStatistics(this.delegate.copy());
    }

    /** {@inheritDoc} */
    @Override
    public synchronized void addValue(double v) {
        delegate.addValue(v);
    }

    /** {@inheritDoc} */
    @Override
    public synchronized void removeMostRecentValue()
        throws MathIllegalStateException {
        delegate.removeMostRecentValue();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double replaceMostRecentValue(double v)
        throws MathIllegalStateException {
        return delegate.replaceMostRecentValue(v);
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double apply(UnivariateStatistic stat) {
        return delegate.apply(stat);
    }

    /** {@inheritDoc} */
    @Override
    public synchronized void clear() {
        delegate.clear();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getElement(int index) {
        return delegate.getElement(index);
    }

    /** {@inheritDoc} */
    @Override
    public synchronized long getN() {
        return delegate.getN();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getMean() {
        return delegate.getMean();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getVariance() {
        return delegate.getVariance();
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
    public synchronized double getSum() {
        return delegate.getSum();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getGeometricMean() {
        return delegate.getGeometricMean();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getPopulationVariance() {
        return delegate.getPopulationVariance();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getSkewness() {
        return delegate.getSkewness();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getKurtosis() {
        return delegate.getKurtosis();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double getSumOfSquares() {
        return delegate.getSumOfSquares();
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
    public synchronized double getPercentile(double p)
        throws MathIllegalArgumentException {
        return delegate.getPercentile(p);
    }

    /** {@inheritDoc} */
    @Override
    public synchronized double[] getValues() {
        return delegate.getValues();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized int getWindowSize() {
        return delegate.getWindowSize();
    }

    /** {@inheritDoc} */
    @Override
    public synchronized void setWindowSize(int windowSize) throws MathIllegalArgumentException {
        delegate.setWindowSize(windowSize);
    }

    /** {@inheritDoc} */
    @Override
    public synchronized String toString() {
        return delegate.toString();
    }

}
