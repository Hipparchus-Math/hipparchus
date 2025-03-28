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
package org.hipparchus.stat.descriptive.moment;

import java.io.Serializable;

import org.hipparchus.exception.NullArgumentException;


/**
 * Computes a statistic related to the Third Central Moment.  Specifically,
 * what is computed is the sum of cubed deviations from the sample mean.
 * <p>
 * The following recursive updating formula is used:
 * <p>
 * Let <ul>
 * <li> dev = (current obs - previous mean) </li>
 * <li> m2 = previous value of {@link SecondMoment} </li>
 * <li> n = number of observations (including current obs) </li>
 * </ul>
 * Then
 * <p>
 * new value = old value - 3 * (dev/n) * m2 + (n-1) * (n -2) * (dev^3/n^2)
 * <p>
 * Returns <code>Double.NaN</code> if no data values have been added and
 * returns <code>0</code> if there is just one value in the data set.
 * Note that Double.NaN may also be returned if the input includes NaN
 * and / or infinite values.
 * <p>
 * <strong>Note that this implementation is not synchronized.</strong> If
 * multiple threads access an instance of this class concurrently, and at least
 * one of the threads invokes the <code>increment()</code> or
 * <code>clear()</code> method, it must be synchronized externally.
 */
class ThirdMoment extends SecondMoment implements Serializable {

    /** Serializable version identifier */
    private static final long serialVersionUID = 20150412L;

    /** third moment of values that have been added */
    protected double m3;

    /**
     * Square of deviation of most recently added value from previous first
     * moment, normalized by previous sample size.  Retained to prevent
     * repeated computation in higher order moments.  nDevSq = nDev * nDev.
     */
    protected double nDevSq;

    /**
     * Create a ThirdMoment instance.
     */
    ThirdMoment() {
        super();
        m3 = Double.NaN;
        nDevSq = Double.NaN;
    }

    /**
     * Copy constructor, creates a new {@code ThirdMoment} identical
     * to the {@code original}.
     *
     * @param original the {@code ThirdMoment} instance to copy
     * @throws NullArgumentException if original is null
     */
    ThirdMoment(ThirdMoment original) throws NullArgumentException {
        super(original);
        this.m3     = original.m3;
        this.nDevSq = original.nDevSq;
    }

    /** {@inheritDoc} */
    @Override
    public void increment(final double d) {
        if (n < 1) {
            m3 = m2 = m1 = 0.0;
        }

        double prevM2 = m2;
        super.increment(d);
        nDevSq = nDev * nDev;
        double n0 = n;
        m3 = m3 - 3.0 * nDev * prevM2 + (n0 - 1) * (n0 - 2) * nDevSq * dev;
    }

    /** {@inheritDoc} */
    @Override
    public double getResult() {
        return m3;
    }

    /** {@inheritDoc} */
    @Override
    public void clear() {
        super.clear();
        m3 = Double.NaN;
        nDevSq = Double.NaN;
    }

    /**
     * Throws {@link UnsupportedOperationException}.
     */
    @Override
    public void aggregate(SecondMoment other) {
        throw new UnsupportedOperationException();
    }

    /** {@inheritDoc} */
    @Override
    public ThirdMoment copy() {
        return new ThirdMoment(this);
    }

}
