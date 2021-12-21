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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.fraction;

import java.util.Iterator;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;

/** Iterator for generating convergents.
 */
class ConvergentsIterator<T> implements Iterator<T> {

    private static final long OVERFLOW = Integer.MAX_VALUE;

    /** Value towards which convergents should converge. */
    private final double value;

    /** Maximum number of convergents to generate. */
    private final int maxConvergents;

    /** Builder for the convergents. */
    private final Builder<T> builder;

    /** Numerator of previous convergent. */
    private long    p0;

    /** Denominator of previous convergent. */
    private long    q0;

    /** Numerator of current convergent. */
    private long    p1;

    /** Denominator of current convergent. */
    private long    q1;

    /** Remainder of current convergent. */
    private double  r1;

    /** Stop indicator. */
    private boolean stop;

    /** Count of already generated convergents. */
    private int     n;

    /** Simple constructor.
     * @param value value towards which convergents should converge
     * @param maxConvergents maximum number of convergents to generate
     * @param builder builder for the convergents
     */
    ConvergentsIterator(final double value, final int maxConvergents, final Builder<T> builder) {

        this.value          = value;
        this.maxConvergents = maxConvergents;
        this.builder        = builder;

        // initialize iterations
        p0   = 0;
        q0   = 1;
        p1   = 1;
        q1   = 0;
        r1   = value;
        stop = false;
        n    = 0;

    }

    /** {@inheritDoc} */
    @Override
    public boolean hasNext() {
        return n < maxConvergents && !stop;
    }

    /** {@inheritDoc} */
    @Override
    public T next() {
        ++n;

        final long a1 = (long) FastMath.floor(r1);
        long p2 = (a1 * p1) + p0;
        long q2 = (a1 * q1) + q0;

        final double convergent = (double) p2 / (double) q2;
        stop = Precision.equals(convergent, value, 1);
        if ((p2 > OVERFLOW || q2 > OVERFLOW) && !stop) {
            throw new MathIllegalStateException(LocalizedCoreFormats.FRACTION_CONVERSION_OVERFLOW, value, p2, q2);
        }
        p0 = p1;
        p1 = p2;
        q0 = q1;
        q1 = q2;
        r1 = 1.0 / (r1 - a1);
        return builder.build(p2, q2);

    }

    /** Interface for building convergents.
     * @param <T> type of the convergent
     */
    public interface Builder<T> {
        /** Build a convergent.
         * @param numerator numerator of the convergent
         * @param denominator denominator of the convergent
         * @return convergent
         */
        T build(long numerator, long denominator);
    }

}
