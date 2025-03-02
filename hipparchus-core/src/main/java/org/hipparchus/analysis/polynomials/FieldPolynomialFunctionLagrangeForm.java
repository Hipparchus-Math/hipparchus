/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.analysis.polynomials;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;

/**
 * Implements the representation of a real polynomial function in
 * <a href="http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html">
 * Lagrange Form</a>. For reference, see <b>Introduction to Numerical
 * Analysis</b>, ISBN 038795452X, chapter 2.
 * <p>
 * The approximated function should be smooth enough for Lagrange polynomial
 * to work well. Otherwise, consider using splines instead.</p>
 * @see PolynomialFunctionLagrangeForm
 * @since 4.0
 * @param <T> type of the field elements
 */
public class FieldPolynomialFunctionLagrangeForm<T extends CalculusFieldElement<T>>
        implements CalculusFieldUnivariateFunction<T> {
    /**
     * The coefficients of the polynomial, ordered by degree -- i.e.
     * coefficients[0] is the constant term and coefficients[n] is the
     * coefficient of x^n where n is the degree of the polynomial.
     */
    private T[] coefficients;
    /**
     * Interpolating points (abscissas).
     */
    private final T[] x;
    /**
     * Function values at interpolating points.
     */
    private final T[] y;
    /**
     * Whether the polynomial coefficients are available.
     */
    private boolean coefficientsComputed;

    /**
     * Construct a Lagrange polynomial with the given abscissas and function
     * values. The order of interpolating points is important.
     * <p>
     * The constructor makes copy of the input arrays and assigns them.</p>
     *
     * @param x interpolating points
     * @param y function values at interpolating points
     * @throws MathIllegalArgumentException if the array lengths are different.
     * @throws MathIllegalArgumentException if the number of points is less than 2.
     * @throws MathIllegalArgumentException if two abscissae have the same value.
     * @throws MathIllegalArgumentException if the abscissae are not sorted.
     */
    public FieldPolynomialFunctionLagrangeForm(final T[] x, final T[] y)
        throws MathIllegalArgumentException {
        this.x = x.clone();
        this.y = y.clone();
        coefficientsComputed = false;

        MathArrays.checkEqualLength(x, y);
        if (x.length < 2) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.WRONG_NUMBER_OF_POINTS, 2, x.length, true);
        }
        MathArrays.checkOrder(x, MathArrays.OrderDirection.INCREASING, true, true);
    }

    /**
     * Calculate the function value at the given point.
     *
     * @param z Point at which the function value is to be computed.
     * @return the function value.
     * @throws MathIllegalArgumentException if {@code x} and {@code y} have
     * different lengths.
     * @throws MathIllegalArgumentException
     * if {@code x} is not sorted in strictly increasing order.
     * @throws MathIllegalArgumentException if the size of {@code x} is less
     * than 2.
     */
    @Override
    public T value(final T z) {
        int nearest = 0;
        final int n = x.length;
        final T[] c = y.clone();
        final T[] d = c.clone();
        double minDist = Double.POSITIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            // find out the abscissa closest to z
            final double dist = FastMath.abs(z.subtract(x[i])).getReal();
            if (dist < minDist) {
                nearest = i;
                minDist = dist;
            }
        }

        // initial approximation to the function value at z
        T value = y[nearest];

        for (int i = 1; i < n; i++) {
            for (int j = 0; j < n-i; j++) {
                final T tc = x[j].subtract(z);
                final T td = x[i+j].subtract(z);
                final T divider = x[j].subtract(x[i+j]);
                // update the difference arrays
                final T w = (c[j+1].subtract(d[j])).divide(divider);
                c[j] = tc.multiply(w);
                d[j] = td.multiply(w);
            }
            // sum up the difference terms to get the final value
            if (nearest < 0.5*(n-i+1)) {
                value = value.add(c[nearest]);    // fork down
            } else {
                nearest--;
                value = value.add(d[nearest]);    // fork up
            }
        }

        return value;
    }

    /**
     * Returns the degree of the polynomial.
     *
     * @return the degree of the polynomial
     */
    public int degree() {
        return x.length - 1;
    }

    /**
     * Returns a copy of the interpolating points array.
     * <p>
     * Changes made to the returned copy will not affect the polynomial.</p>
     *
     * @return a fresh copy of the interpolating points array
     */
    public T[] getInterpolatingPoints() {
        return x.clone();
    }

    /**
     * Returns a copy of the interpolating values array.
     * <p>
     * Changes made to the returned copy will not affect the polynomial.</p>
     *
     * @return a fresh copy of the interpolating values array
     */
    public T[] getInterpolatingValues() {
        return y.clone();
    }

    /**
     * Returns a copy of the coefficients array.
     * <p>
     * Changes made to the returned copy will not affect the polynomial.</p>
     * <p>
     * Note that coefficients computation can be ill-conditioned. Use with caution
     * and only when it is necessary.</p>
     *
     * @return a fresh copy of the coefficients array
     */
    public T[] getCoefficients() {
        if (!coefficientsComputed) {
            computeCoefficients();
        }
        return coefficients.clone();
    }

    /**
     * Calculate the coefficients of Lagrange polynomial from the
     * interpolation data. It takes O(n^2) time.
     * Note that this computation can be ill-conditioned: Use with caution
     * and only when it is necessary.
     */
    protected void computeCoefficients() {
        final int n = degree() + 1;
        final Field<T> field = x[0].getField();
        coefficients = MathArrays.buildArray(field, n);

        // c[] are the coefficients of P(x) = (x-x[0])(x-x[1])...(x-x[n-1])
        final T[] c = MathArrays.buildArray(field, n + 1);
        c[0] = field.getOne();
        for (int i = 0; i < n; i++) {
            for (int j = i; j > 0; j--) {
                c[j] = c[j-1].subtract(c[j].multiply(x[i]));
            }
            c[0] = c[0].multiply(x[i].negate());
            c[i+1] = field.getOne();
        }

        final T[] tc = MathArrays.buildArray(field, n);
        for (int i = 0; i < n; i++) {
            // d = (x[i]-x[0])...(x[i]-x[i-1])(x[i]-x[i+1])...(x[i]-x[n-1])
            T d = field.getOne();
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    d = d.multiply(x[i].subtract(x[j]));
                }
            }
            final T t = y[i].divide(d);
            // Lagrange polynomial is the sum of n terms, each of which is a
            // polynomial of degree n-1. tc[] are the coefficients of the i-th
            // numerator Pi(x) = (x-x[0])...(x-x[i-1])(x-x[i+1])...(x-x[n-1]).
            tc[n-1] = c[n];     // actually c[n] = 1
            coefficients[n-1] = coefficients[n-1].add(t.multiply(tc[n-1]));
            for (int j = n-2; j >= 0; j--) {
                tc[j] = c[j+1].add(tc[j+1].multiply(x[i]));
                coefficients[j] = coefficients[j].add(t.multiply(tc[j]));
            }
        }

        coefficientsComputed = true;
    }
}
