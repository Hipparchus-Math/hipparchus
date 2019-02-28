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
package org.hipparchus.analysis.interpolation;

import java.lang.reflect.Array;

import org.hipparchus.RealFieldElement;
import org.hipparchus.analysis.polynomials.FieldPolynomialFunction;
import org.hipparchus.analysis.polynomials.FieldPolynomialSplineFunction;
import org.hipparchus.analysis.polynomials.PolynomialFunction;
import org.hipparchus.analysis.polynomials.PolynomialSplineFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;

/**
 * Implements a linear function for interpolation of real univariate functions.
 *
 */
public class LinearInterpolator implements UnivariateInterpolator, FieldUnivariateInterpolator {

    /**
     * Computes a linear interpolating function for the data set.
     *
     * @param x the arguments for the interpolation points
     * @param y the values for the interpolation points
     * @return a function which interpolates the data set
     * @throws MathIllegalArgumentException if {@code x} and {@code y}
     * have different sizes.
     * @throws MathIllegalArgumentException if {@code x} is not sorted in
     * strict increasing order.
     * @throws MathIllegalArgumentException if the size of {@code x} is smaller
     * than 2.
     */
    @Override
    public PolynomialSplineFunction interpolate(double x[], double y[])
        throws MathIllegalArgumentException {
        MathUtils.checkNotNull(x);
        MathUtils.checkNotNull(y);
        MathArrays.checkEqualLength(x, y);

        if (x.length < 2) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_OF_POINTS,
                                                x.length, 2, true);
        }

        // Number of intervals.  The number of data points is n + 1.
        int n = x.length - 1;

        MathArrays.checkOrder(x);

        // Slope of the lines between the datapoints.
        final double m[] = new double[n];
        for (int i = 0; i < n; i++) {
            m[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        }

        final PolynomialFunction polynomials[] = new PolynomialFunction[n];
        final double coefficients[] = new double[2];
        for (int i = 0; i < n; i++) {
            coefficients[0] = y[i];
            coefficients[1] = m[i];
            polynomials[i] = new PolynomialFunction(coefficients);
        }

        return new PolynomialSplineFunction(x, polynomials);
    }

    /**
     * Computes a linear interpolating function for the data set.
     *
     * @param x the arguments for the interpolation points
     * @param y the values for the interpolation points
     * @param <T> the type of the field elements
     * @return a function which interpolates the data set
     * @throws MathIllegalArgumentException if {@code x} and {@code y}
     * have different sizes.
     * @throws MathIllegalArgumentException if {@code x} is not sorted in
     * strict increasing order.
     * @throws MathIllegalArgumentException if the size of {@code x} is smaller
     * than 2.
     * @since 1.5
     */
    @Override
    public <T extends RealFieldElement<T>> FieldPolynomialSplineFunction<T> interpolate(final T x[], final T y[])
        throws MathIllegalArgumentException {
        MathUtils.checkNotNull(x);
        MathUtils.checkNotNull(y);
        MathArrays.checkEqualLength(x, y);

        if (x.length < 2) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_OF_POINTS,
                                                x.length, 2, true);
        }

        // Number of intervals.  The number of data points is n + 1.
        int n = x.length - 1;

        MathArrays.checkOrder(x);

        // Slope of the lines between the datapoints.
        final T m[] = MathArrays.buildArray(x[0].getField(), n);
        for (int i = 0; i < n; i++) {
            m[i] = y[i + 1].subtract(y[i]).divide(x[i + 1].subtract(x[i]));
        }

        @SuppressWarnings("unchecked")
        final FieldPolynomialFunction<T> polynomials[] =
                        (FieldPolynomialFunction<T>[]) Array.newInstance(FieldPolynomialFunction.class, n);
        final T coefficients[] = MathArrays.buildArray(x[0].getField(), 2);
        for (int i = 0; i < n; i++) {
            coefficients[0] = y[i];
            coefficients[1] = m[i];
            polynomials[i] = new FieldPolynomialFunction<>(coefficients);
        }

        return new FieldPolynomialSplineFunction<>(x, polynomials);
    }

}
