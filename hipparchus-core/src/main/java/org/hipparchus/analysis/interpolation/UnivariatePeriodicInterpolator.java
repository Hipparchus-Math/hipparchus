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
package org.hipparchus.analysis.interpolation;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;

/**
 * Adapter for classes implementing the {@link UnivariateInterpolator}
 * interface.
 * The data to be interpolated is assumed to be periodic. Thus values that are
 * outside of the range can be passed to the interpolation function: They will
 * be wrapped into the initial range before being passed to the class that
 * actually computes the interpolation.
 *
 */
public class UnivariatePeriodicInterpolator
    implements UnivariateInterpolator {
    /** Default number of extension points of the samples array. */
    public static final int DEFAULT_EXTEND = 5;
    /** Interpolator. */
    private final UnivariateInterpolator interpolator;
    /** Period. */
    private final double period;
    /** Number of extension points. */
    private final int extend;

    /**
     * Builds an interpolator.
     *
     * @param interpolator Interpolator.
     * @param period Period.
     * @param extend Number of points to be appended at the beginning and
     * end of the sample arrays in order to avoid interpolation failure at
     * the (periodic) boundaries of the orginal interval. The value is the
     * number of sample points which the original {@code interpolator} needs
     * on each side of the interpolated point.
     */
    public UnivariatePeriodicInterpolator(UnivariateInterpolator interpolator,
                                          double period,
                                          int extend) {
        this.interpolator = interpolator;
        this.period = period;
        this.extend = extend;
    }

    /**
     * Builds an interpolator.
     * Uses {@link #DEFAULT_EXTEND} as the number of extension points on each side
     * of the original abscissae range.
     *
     * @param interpolator Interpolator.
     * @param period Period.
     */
    public UnivariatePeriodicInterpolator(UnivariateInterpolator interpolator,
                                          double period) {
        this(interpolator, period, DEFAULT_EXTEND);
    }

    /**
     * {@inheritDoc}
     *
     * @throws MathIllegalArgumentException if the number of extension points
     * is larger than the size of {@code xval}.
     */
    @Override
    public UnivariateFunction interpolate(double[] xval,
                                          double[] yval)
        throws MathIllegalArgumentException {
        if (xval.length < extend) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL,
                                                   xval.length, extend);
        }

        MathArrays.checkOrder(xval);
        final double offset = xval[0];

        final int len = xval.length + extend * 2;
        final double[] x = new double[len];
        final double[] y = new double[len];
        for (int i = 0; i < xval.length; i++) {
            final int index = i + extend;
            x[index] = MathUtils.reduce(xval[i], period, offset);
            y[index] = yval[i];
        }

        // Wrap to enable interpolation at the boundaries.
        for (int i = 0; i < extend; i++) {
            int index = xval.length - extend + i;
            x[i] = MathUtils.reduce(xval[index], period, offset) - period;
            y[i] = yval[index];

            index = len - extend + i;
            x[index] = MathUtils.reduce(xval[i], period, offset) + period;
            y[index] = yval[i];
        }

        MathArrays.sortInPlace(x, y);

        final UnivariateFunction f = interpolator.interpolate(x, y);
        return new UnivariateFunction() {
            /** {@inheritDoc} */
            @Override
            public double value(final double x) throws MathIllegalArgumentException {
                return f.value(MathUtils.reduce(x, period, offset));
            }
        };
    }
}
