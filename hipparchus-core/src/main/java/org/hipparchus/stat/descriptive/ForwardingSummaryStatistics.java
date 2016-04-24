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
import org.hipparchus.util.MathUtils;

/**
 * A {@link SummaryStatistics} which forwards all its method calls to another
 * {@link SummaryStatistics}.
 * <p>
 * Subclasses should override one or more methods to modify the behavior of the backing
 * summary statistics as desired per the decorator pattern.
 */
abstract class ForwardingSummaryStatistics
    implements SummaryStatistics {

    /** The SummaryStatistics instance to delegate to. */
    private final SummaryStatistics delegate;

    /**
     * Constructor that wraps the given instance.
     *
     * @param original the {@code SummaryStatistics} instance to wrap
     * @throws NullArgumentException if original is null
     */
    ForwardingSummaryStatistics(SummaryStatistics original)
        throws NullArgumentException {
        MathUtils.checkNotNull(original);
        this.delegate = original;
    }

    /**
     * Returns the delegate instance.
     * @return the delegate instance
     */
    protected SummaryStatistics delegate() {
        return delegate;
    }

    /** {@inheritDoc} */
    @Override
    public StatisticalSummary getSummary() {
        return delegate.getSummary();
    }

    /** {@inheritDoc} */
    @Override
    public void addValue(double value) {
        delegate.addValue(value);
    }

    /** {@inheritDoc} */
    @Override
    public void clear() {
        delegate.clear();
    }

    /** {@inheritDoc} */
    @Override
    public long getN() {
        return delegate.getN();
    }

    /** {@inheritDoc} */
    @Override
    public double getSum() {
        return delegate.getSum();
    }

    /** {@inheritDoc} */
    @Override
    public double getSumOfSquares() {
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
    public double getMean() {
        return delegate.getMean();
    }

    /** {@inheritDoc} */
    @Override
    public double getStandardDeviation() {
        return delegate.getStandardDeviation();
    }

    /** {@inheritDoc} */
    @Override
    public double getQuadraticMean() {
        return delegate.getQuadraticMean();
    }

    /** {@inheritDoc} */
    @Override
    public double getVariance() {
        return delegate.getVariance();
    }

    /** {@inheritDoc} */
    @Override
    public double getPopulationVariance() {
        return delegate.getPopulationVariance();
    }

    /** {@inheritDoc} */
    @Override
    public double getMax() {
        return delegate.getMax();
    }

    /** {@inheritDoc} */
    @Override
    public double getMin() {
        return delegate.getMin();
    }

    /** {@inheritDoc} */
    @Override
    public double getGeometricMean() {
        return delegate.getGeometricMean();
    }

    /** {@inheritDoc} */
    @Override
    public String toString() {
        return delegate.toString();
    }

    /** {@inheritDoc} */
    @Override
    public boolean equals(Object object) {
        return delegate.equals(object);
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return delegate.hashCode();
    }

}
