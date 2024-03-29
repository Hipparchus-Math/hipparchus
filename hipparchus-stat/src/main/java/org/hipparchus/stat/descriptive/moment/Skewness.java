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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.stat.descriptive.AbstractStorelessUnivariateStatistic;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;

/**
 * Computes the skewness of the available values.
 * <p>
 * We use the following (unbiased) formula to define skewness:
 * <p>
 * skewness = [n / (n -1) (n - 2)] sum[(x_i - mean)^3] / std^3
 * <p>
 * where n is the number of values, mean is the {@link Mean} and std is the
 * {@link StandardDeviation}.
 * <p>
 * Note that this statistic is undefined for n &lt; 3.  <code>Double.Nan</code>
 * is returned when there is not sufficient data to compute the statistic.
 * Double.NaN may also be returned if the input includes NaN and / or
 * infinite values.
 * <p>
 * <strong>Note that this implementation is not synchronized.</strong> If
 * multiple threads access an instance of this class concurrently, and at least
 * one of the threads invokes the <code>increment()</code> or
 * <code>clear()</code> method, it must be synchronized externally.
 */
public class Skewness extends AbstractStorelessUnivariateStatistic implements Serializable {

    /** Serializable version identifier */
    private static final long serialVersionUID = 20150412L;

    /** Third moment on which this statistic is based */
    protected final ThirdMoment moment;

     /**
     * Determines whether or not this statistic can be incremented or cleared.
     * <p>
     * Statistics based on (constructed from) external moments cannot
     * be incremented or cleared.
    */
    protected final boolean incMoment;

    /**
     * Constructs a Skewness.
     */
    public Skewness() {
        moment = new ThirdMoment();
        incMoment = true;
    }

    /**
     * Constructs a Skewness with an external moment.
     * @param m3 external moment
     */
    public Skewness(final ThirdMoment m3) {
        this.moment = m3;
        incMoment = false;
    }

    /**
     * Copy constructor, creates a new {@code Skewness} identical
     * to the {@code original}.
     *
     * @param original the {@code Skewness} instance to copy
     * @throws NullArgumentException if original is null
     */
    public Skewness(Skewness original) throws NullArgumentException {
        MathUtils.checkNotNull(original);
        this.moment    = original.moment.copy();
        this.incMoment = original.incMoment;
    }

    /**
     * {@inheritDoc}
     * <p>Note that when {@link #Skewness(ThirdMoment)} is used to
     * create a Skewness, this method does nothing. In that case, the
     * ThirdMoment should be incremented directly.
     */
    @Override
    public void increment(final double d) {
        if (incMoment) {
            moment.increment(d);
        }
    }

    /**
     * Returns the value of the statistic based on the values that have been added.
     * <p>
     * See {@link Skewness} for the definition used in the computation.
     *
     * @return the skewness of the available values.
     */
    @Override
    public double getResult() {

        if (moment.n < 3) {
            return Double.NaN;
        }
        double variance = moment.m2 / (moment.n - 1);
        if (variance < 10E-20) {
            return 0.0d;
        } else {
            double n0 = moment.getN();
            return  (n0 * moment.m3) /
            ((n0 - 1) * (n0 -2) * FastMath.sqrt(variance) * variance);
        }
    }

    /** {@inheritDoc} */
    @Override
    public long getN() {
        return moment.getN();
    }

    /** {@inheritDoc} */
    @Override
    public void clear() {
        if (incMoment) {
            moment.clear();
        }
    }

    /**
     * Returns the Skewness of the entries in the specified portion of the
     * input array.
     * <p>
     * See {@link Skewness} for the definition used in the computation.
     * <p>
     * Throws <code>IllegalArgumentException</code> if the array is null.
     *
     * @param values the input array
     * @param begin the index of the first array element to include
     * @param length the number of elements to include
     * @return the skewness of the values or Double.NaN if length is less than 3
     * @throws MathIllegalArgumentException if the array is null or the array index
     *  parameters are not valid
     */
    @Override
    public double evaluate(final double[] values, final int begin, final int length)
        throws MathIllegalArgumentException {

        // Initialize the skewness
        double skew = Double.NaN;

        if (MathArrays.verifyValues(values, begin, length) && length > 2 ) {
            Mean mean = new Mean();
            // Get the mean and the standard deviation
            double m = mean.evaluate(values, begin, length);

            // Calc the std, this is implemented here instead
            // of using the standardDeviation method eliminate
            // a duplicate pass to get the mean
            double accum = 0.0;
            double accum2 = 0.0;
            for (int i = begin; i < begin + length; i++) {
                final double d = values[i] - m;
                accum  += d * d;
                accum2 += d;
            }
            final double variance = (accum - (accum2 * accum2 / length)) / (length - 1);

            double accum3 = 0.0;
            for (int i = begin; i < begin + length; i++) {
                final double d = values[i] - m;
                accum3 += d * d * d;
            }
            accum3 /= variance * FastMath.sqrt(variance);

            // Get N
            double n0 = length;

            // Calculate skewness
            skew = (n0 / ((n0 - 1) * (n0 - 2))) * accum3;
        }
        return skew;
    }

    /** {@inheritDoc} */
    @Override
    public Skewness copy() {
        return new Skewness(this);
    }

}
