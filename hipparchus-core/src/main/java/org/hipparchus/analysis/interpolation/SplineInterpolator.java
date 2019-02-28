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

import org.hipparchus.Field;
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
 * Computes a natural (also known as "free", "unclamped") cubic spline interpolation for the data set.
 * <p>
 * The {@link #interpolate(double[], double[])} method returns a {@link PolynomialSplineFunction}
 * consisting of n cubic polynomials, defined over the subintervals determined by the x values,
 * {@code x[0] < x[i] ... < x[n].}  The x values are referred to as "knot points."</p>
 * <p>
 * The value of the PolynomialSplineFunction at a point x that is greater than or equal to the smallest
 * knot point and strictly less than the largest knot point is computed by finding the subinterval to which
 * x belongs and computing the value of the corresponding polynomial at <code>x - x[i] </code> where
 * <code>i</code> is the index of the subinterval.  See {@link PolynomialSplineFunction} for more details.
 * </p>
 * <p>
 * The interpolating polynomials satisfy: <ol>
 * <li>The value of the PolynomialSplineFunction at each of the input x values equals the
 *  corresponding y value.</li>
 * <li>Adjacent polynomials are equal through two derivatives at the knot points (i.e., adjacent polynomials
 *  "match up" at the knot points, as do their first and second derivatives).</li>
 * </ol></p>
 * <p>
 * The cubic spline interpolation algorithm implemented is as described in R.L. Burden, J.D. Faires,
 * <u>Numerical Analysis</u>, 4th Ed., 1989, PWS-Kent, ISBN 0-53491-585-X, pp 126-131.
 * </p>
 *
 */
public class SplineInterpolator implements UnivariateInterpolator, FieldUnivariateInterpolator {

    /**
     * Computes an interpolating function for the data set.
     * @param x the arguments for the interpolation points
     * @param y the values for the interpolation points
     * @return a function which interpolates the data set
     * @throws MathIllegalArgumentException if {@code x} and {@code y}
     * have different sizes.
     * @throws MathIllegalArgumentException if {@code x} is not sorted in
     * strict increasing order.
     * @throws MathIllegalArgumentException if the size of {@code x} is smaller
     * than 3.
     */
    @Override
    public PolynomialSplineFunction interpolate(double x[], double y[])
        throws MathIllegalArgumentException {

        MathUtils.checkNotNull(x);
        MathUtils.checkNotNull(y);
        MathArrays.checkEqualLength(x, y);
        if (x.length < 3) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_OF_POINTS,
                                                   x.length, 3, true);
        }

        // Number of intervals.  The number of data points is n + 1.
        final int n = x.length - 1;

        MathArrays.checkOrder(x);

        // Differences between knot points
        final double h[] = new double[n];
        for (int i = 0; i < n; i++) {
            h[i] = x[i + 1] - x[i];
        }

        final double mu[] = new double[n];
        final double z[] = new double[n + 1];
        mu[0] = 0d;
        z[0] = 0d;
        double g = 0;
        for (int i = 1; i < n; i++) {
            g = 2d * (x[i+1]  - x[i - 1]) - h[i - 1] * mu[i -1];
            mu[i] = h[i] / g;
            z[i] = (3d * (y[i + 1] * h[i - 1] - y[i] * (x[i + 1] - x[i - 1])+ y[i - 1] * h[i]) /
                    (h[i - 1] * h[i]) - h[i - 1] * z[i - 1]) / g;
        }

        // cubic spline coefficients --  b is linear, c quadratic, d is cubic (original y's are constants)
        final double b[] = new double[n];
        final double c[] = new double[n + 1];
        final double d[] = new double[n];

        z[n] = 0d;
        c[n] = 0d;

        for (int j = n -1; j >=0; j--) {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2d * c[j]) / 3d;
            d[j] = (c[j + 1] - c[j]) / (3d * h[j]);
        }

        final PolynomialFunction polynomials[] = new PolynomialFunction[n];
        final double coefficients[] = new double[4];
        for (int i = 0; i < n; i++) {
            coefficients[0] = y[i];
            coefficients[1] = b[i];
            coefficients[2] = c[i];
            coefficients[3] = d[i];
            polynomials[i] = new PolynomialFunction(coefficients);
        }

        return new PolynomialSplineFunction(x, polynomials);
    }

    /**
     * Computes an interpolating function for the data set.
     * @param x the arguments for the interpolation points
     * @param y the values for the interpolation points
     * @param <T> the type of the field elements
     * @return a function which interpolates the data set
     * @throws MathIllegalArgumentException if {@code x} and {@code y}
     * have different sizes.
     * @throws MathIllegalArgumentException if {@code x} is not sorted in
     * strict increasing order.
     * @throws MathIllegalArgumentException if the size of {@code x} is smaller
     * than 3.
     * @since 1.5
     */
    @Override
    public <T extends RealFieldElement<T>> FieldPolynomialSplineFunction<T> interpolate(T x[], T y[])
        throws MathIllegalArgumentException {

        MathUtils.checkNotNull(x);
        MathUtils.checkNotNull(y);
        MathArrays.checkEqualLength(x, y);
        if (x.length < 3) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_OF_POINTS,
                                                   x.length, 3, true);
        }

        // Number of intervals.  The number of data points is n + 1.
        final int n = x.length - 1;

        MathArrays.checkOrder(x);

        // Differences between knot points
        final Field<T> field = x[0].getField();
        final T h[] = MathArrays.buildArray(field, n);
        for (int i = 0; i < n; i++) {
            h[i] = x[i + 1].subtract(x[i]);
        }

        final T mu[] = MathArrays.buildArray(field, n);
        final T z[]  = MathArrays.buildArray(field, n + 1);
        mu[0] = field.getZero();
        z[0]  = field.getZero();
        for (int i = 1; i < n; i++) {
            final T g = x[i+1].subtract(x[i - 1]).multiply(2).subtract(h[i - 1].multiply(mu[i -1]));
            mu[i] = h[i].divide(g);
            z[i] =          y[i + 1].multiply(h[i - 1]).
                   subtract(y[i].multiply(x[i + 1].subtract(x[i - 1]))).
                        add(y[i - 1].multiply(h[i])).
                   multiply(3).
                   divide(h[i - 1].multiply(h[i])).
                   subtract(h[i - 1].multiply(z[i - 1])).
                   divide(g);
        }

        // cubic spline coefficients --  b is linear, c quadratic, d is cubic (original y's are constants)
        final T b[] = MathArrays.buildArray(field, n);
        final T c[] = MathArrays.buildArray(field, n + 1);
        final T d[] = MathArrays.buildArray(field, n);

        z[n] = field.getZero();
        c[n] = field.getZero();

        for (int j = n -1; j >=0; j--) {
            c[j] = z[j].subtract(mu[j].multiply(c[j + 1]));
            b[j] = y[j + 1].subtract(y[j]).divide(h[j]).
                   subtract(h[j].multiply(c[j + 1].add(c[j]).add(c[j])).divide(3));
            d[j] = c[j + 1].subtract(c[j]).divide(h[j].multiply(3));
        }

        @SuppressWarnings("unchecked")
        final FieldPolynomialFunction<T> polynomials[] =
                        (FieldPolynomialFunction<T>[]) Array.newInstance(FieldPolynomialFunction.class, n);
        final T coefficients[] = MathArrays.buildArray(field, 4);
        for (int i = 0; i < n; i++) {
            coefficients[0] = y[i];
            coefficients[1] = b[i];
            coefficients[2] = c[i];
            coefficients[3] = d[i];
            polynomials[i] = new FieldPolynomialFunction<>(coefficients);
        }

        return new FieldPolynomialSplineFunction<>(x, polynomials);
    }

}
