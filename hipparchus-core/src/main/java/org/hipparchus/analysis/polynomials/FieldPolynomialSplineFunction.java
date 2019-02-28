/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.analysis.polynomials;

import java.lang.reflect.Array;
import java.util.Arrays;

import org.hipparchus.Field;
import org.hipparchus.RealFieldElement;
import org.hipparchus.analysis.RealFieldUnivariateFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;

/**
 * Represents a polynomial spline function.
 * <p>
 * A <strong>polynomial spline function</strong> consists of a set of
 * <i>interpolating polynomials</i> and an ascending array of domain
 * <i>knot points</i>, determining the intervals over which the spline function
 * is defined by the constituent polynomials.  The polynomials are assumed to
 * have been computed to match the values of another function at the knot
 * points.  The value consistency constraints are not currently enforced by
 * <code>PolynomialSplineFunction</code> itself, but are assumed to hold among
 * the polynomials and knot points passed to the constructor.</p>
 * <p>
 * N.B.:  The polynomials in the <code>polynomials</code> property must be
 * centered on the knot points to compute the spline function values.
 * See below.</p>
 * <p>
 * The domain of the polynomial spline function is
 * <code>[smallest knot, largest knot]</code>.  Attempts to evaluate the
 * function at values outside of this range generate IllegalArgumentExceptions.
 * </p>
 * <p>
 * The value of the polynomial spline function for an argument <code>x</code>
 * is computed as follows:
 * <ol>
 * <li>The knot array is searched to find the segment to which <code>x</code>
 * belongs.  If <code>x</code> is less than the smallest knot point or greater
 * than the largest one, an <code>IllegalArgumentException</code>
 * is thrown.</li>
 * <li> Let <code>j</code> be the index of the largest knot point that is less
 * than or equal to <code>x</code>.  The value returned is
 * {@code polynomials[j](x - knot[j])}</li></ol>
 *
 * @param <T> the type of the field elements
 * @since 1.5
 */
public class FieldPolynomialSplineFunction<T extends RealFieldElement<T>> implements RealFieldUnivariateFunction<T> {

    /**
     * Spline segment interval delimiters (knots).
     * Size is n + 1 for n segments.
     */
    private final T knots[];

    /**
     * The polynomial functions that make up the spline.  The first element
     * determines the value of the spline over the first subinterval, the
     * second over the second, etc.   Spline function values are determined by
     * evaluating these functions at {@code (x - knot[i])} where i is the
     * knot segment to which x belongs.
     */
    private final FieldPolynomialFunction<T> polynomials[];

    /**
     * Number of spline segments. It is equal to the number of polynomials and
     * to the number of partition points - 1.
     */
    private final int n;


    /**
     * Construct a polynomial spline function with the given segment delimiters
     * and interpolating polynomials.
     * The constructor copies both arrays and assigns the copies to the knots
     * and polynomials properties, respectively.
     *
     * @param knots Spline segment interval delimiters.
     * @param polynomials Polynomial functions that make up the spline.
     * @throws NullArgumentException if either of the input arrays is {@code null}.
     * @throws MathIllegalArgumentException if knots has length less than 2.
     * @throws MathIllegalArgumentException if {@code polynomials.length != knots.length - 1}.
     * @throws MathIllegalArgumentException if the {@code knots} array is not strictly increasing.
     *
     */
    @SuppressWarnings("unchecked")
    public FieldPolynomialSplineFunction(final T knots[], final FieldPolynomialFunction<T> polynomials[])
        throws MathIllegalArgumentException, NullArgumentException {
        if (knots == null ||
            polynomials == null) {
            throw new NullArgumentException();
        }
        if (knots.length < 2) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NOT_ENOUGH_POINTS_IN_SPLINE_PARTITION,
                                                   2, knots.length, false);
        }
        MathUtils.checkDimension(polynomials.length, knots.length - 1);
        MathArrays.checkOrder(knots);

        this.n = knots.length -1;
        this.knots = knots.clone();
        this.polynomials = (FieldPolynomialFunction<T>[]) Array.newInstance(FieldPolynomialFunction.class, n);
        System.arraycopy(polynomials, 0, this.polynomials, 0, n);
    }

    /** Get the {@link Field} to which the instance belongs.
     * @return {@link Field} to which the instance belongs
     */
    public Field<T> getField() {
        return knots[0].getField();
    }

    /**
     * Compute the value for the function.
     * See {@link FieldPolynomialSplineFunction} for details on the algorithm for
     * computing the value of the function.
     *
     * @param v Point for which the function value should be computed.
     * @return the value.
     * @throws MathIllegalArgumentException if {@code v} is outside of the domain of the
     * spline function (smaller than the smallest knot point or larger than the
     * largest knot point).
     */
    public T value(final double v) {
        return value(getField().getZero().add(v));
    }

    /**
     * Compute the value for the function.
     * See {@link FieldPolynomialSplineFunction} for details on the algorithm for
     * computing the value of the function.
     *
     * @param v Point for which the function value should be computed.
     * @return the value.
     * @throws MathIllegalArgumentException if {@code v} is outside of the domain of the
     * spline function (smaller than the smallest knot point or larger than the
     * largest knot point).
     */
    @Override
    public T value(final T v) {
        MathUtils.checkRangeInclusive(v.getReal(), knots[0].getReal(), knots[n].getReal());
        int i = Arrays.binarySearch(knots, v);
        if (i < 0) {
            i = -i - 2;
        }
        // This will handle the case where v is the last knot value
        // There are only n-1 polynomials, so if v is the last knot
        // then we will use the last polynomial to calculate the value.
        if ( i >= polynomials.length ) {
            i--;
        }
        return polynomials[i].value(v.subtract(knots[i]));
    }

    /**
     * Get the number of spline segments.
     * It is also the number of polynomials and the number of knot points - 1.
     *
     * @return the number of spline segments.
     */
    public int getN() {
        return n;
    }

    /**
     * Get a copy of the interpolating polynomials array.
     * It returns a fresh copy of the array. Changes made to the copy will
     * not affect the polynomials property.
     *
     * @return the interpolating polynomials.
     */
    public FieldPolynomialFunction<T>[] getPolynomials() {
        return polynomials.clone();
    }

    /**
     * Get an array copy of the knot points.
     * It returns a fresh copy of the array. Changes made to the copy
     * will not affect the knots property.
     *
     * @return the knot points.
     */
    public T[] getKnots() {
        return knots.clone();
    }

    /**
     * Indicates whether a point is within the interpolation range.
     *
     * @param x Point.
     * @return {@code true} if {@code x} is a valid point.
     */
    public boolean isValidPoint(T x) {
        if (x.getReal() < knots[0].getReal() ||
            x.getReal() > knots[n].getReal()) {
            return false;
        } else {
            return true;
        }
    }
    /**
     * Get the derivative of the polynomial spline function.
     *
     * @return the derivative function.
     */
    @SuppressWarnings("unchecked")
    public FieldPolynomialSplineFunction<T> polynomialSplineDerivative() {
        FieldPolynomialFunction<T> derivativePolynomials[] =
                        (FieldPolynomialFunction<T>[]) Array.newInstance(FieldPolynomialFunction.class, n);
        for (int i = 0; i < n; i++) {
            derivativePolynomials[i] = polynomials[i].polynomialDerivative();
        }
        return new FieldPolynomialSplineFunction<>(knots, derivativePolynomials);
    }


}
