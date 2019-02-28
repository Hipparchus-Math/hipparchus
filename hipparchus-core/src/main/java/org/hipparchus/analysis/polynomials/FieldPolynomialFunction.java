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

import org.hipparchus.Field;
import org.hipparchus.RealFieldElement;
import org.hipparchus.analysis.RealFieldUnivariateFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;

/**
 * Immutable representation of a real polynomial function with real coefficients.
 * <p>
 * <a href="http://mathworld.wolfram.com/HornersMethod.html">Horner's Method</a>
 * is used to evaluate the function.</p>
 * @param <T> the type of the field elements
 * @since 1.5
 *
 */
public class FieldPolynomialFunction<T extends RealFieldElement<T>> implements RealFieldUnivariateFunction<T> {

    /**
     * The coefficients of the polynomial, ordered by degree -- i.e.,
     * coefficients[0] is the constant term and coefficients[n] is the
     * coefficient of x^n where n is the degree of the polynomial.
     */
    private final T coefficients[];

    /**
     * Construct a polynomial with the given coefficients.  The first element
     * of the coefficients array is the constant term.  Higher degree
     * coefficients follow in sequence.  The degree of the resulting polynomial
     * is the index of the last non-null element of the array, or 0 if all elements
     * are null.
     * <p>
     * The constructor makes a copy of the input array and assigns the copy to
     * the coefficients property.</p>
     *
     * @param c Polynomial coefficients.
     * @throws NullArgumentException if {@code c} is {@code null}.
     * @throws MathIllegalArgumentException if {@code c} is empty.
     */
    public FieldPolynomialFunction(final T c[])
        throws MathIllegalArgumentException, NullArgumentException {
        super();
        MathUtils.checkNotNull(c);
        int n = c.length;
        if (n == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
        }
        while ((n > 1) && (c[n - 1].getReal() == 0)) {
            --n;
        }
        this.coefficients = MathArrays.buildArray(c[0].getField(), n);
        System.arraycopy(c, 0, this.coefficients, 0, n);
    }

    /**
     * Compute the value of the function for the given argument.
     * <p>
     *  The value returned is </p><p>
     *  {@code coefficients[n] * x^n + ... + coefficients[1] * x  + coefficients[0]}
     * </p>
     *
     * @param x Argument for which the function value should be computed.
     * @return the value of the polynomial at the given point.
     *
     * @see org.hipparchus.analysis.UnivariateFunction#value(double)
     */
    public T value(double x) {
       return evaluate(coefficients, getField().getZero().add(x));
    }

    /**
     * Compute the value of the function for the given argument.
     * <p>
     *  The value returned is </p><p>
     *  {@code coefficients[n] * x^n + ... + coefficients[1] * x  + coefficients[0]}
     * </p>
     *
     * @param x Argument for which the function value should be computed.
     * @return the value of the polynomial at the given point.
     *
     * @see org.hipparchus.analysis.UnivariateFunction#value(double)
     */
    @Override
    public T value(T x) {
       return evaluate(coefficients, x);
    }

    /** Get the {@link Field} to which the instance belongs.
     * @return {@link Field} to which the instance belongs
     */
    public Field<T> getField() {
        return coefficients[0].getField();
    }

    /**
     * Returns the degree of the polynomial.
     *
     * @return the degree of the polynomial.
     */
    public int degree() {
        return coefficients.length - 1;
    }

    /**
     * Returns a copy of the coefficients array.
     * <p>
     * Changes made to the returned copy will not affect the coefficients of
     * the polynomial.</p>
     *
     * @return a fresh copy of the coefficients array.
     */
    public T[] getCoefficients() {
        return coefficients.clone();
    }

    /**
     * Uses Horner's Method to evaluate the polynomial with the given coefficients at
     * the argument.
     *
     * @param coefficients Coefficients of the polynomial to evaluate.
     * @param argument Input value.
     * @param <T> the type of the field elements
     * @return the value of the polynomial.
     * @throws MathIllegalArgumentException if {@code coefficients} is empty.
     * @throws NullArgumentException if {@code coefficients} is {@code null}.
     */
    protected static <T extends RealFieldElement<T>> T evaluate(T[] coefficients, T argument)
        throws MathIllegalArgumentException, NullArgumentException {
        MathUtils.checkNotNull(coefficients);
        int n = coefficients.length;
        if (n == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
        }
        T result = coefficients[n - 1];
        for (int j = n - 2; j >= 0; j--) {
            result = argument.multiply(result).add(coefficients[j]);
        }
        return result;
    }

    /**
     * Add a polynomial to the instance.
     *
     * @param p Polynomial to add.
     * @return a new polynomial which is the sum of the instance and {@code p}.
     */
    public FieldPolynomialFunction<T> add(final FieldPolynomialFunction<T> p) {
        // identify the lowest degree polynomial
        final int lowLength  = FastMath.min(coefficients.length, p.coefficients.length);
        final int highLength = FastMath.max(coefficients.length, p.coefficients.length);

        // build the coefficients array
        T[] newCoefficients = MathArrays.buildArray(getField(), highLength);
        for (int i = 0; i < lowLength; ++i) {
            newCoefficients[i] = coefficients[i].add(p.coefficients[i]);
        }
        System.arraycopy((coefficients.length < p.coefficients.length) ?
                         p.coefficients : coefficients,
                         lowLength,
                         newCoefficients, lowLength,
                         highLength - lowLength);

        return new FieldPolynomialFunction<>(newCoefficients);
    }

    /**
     * Subtract a polynomial from the instance.
     *
     * @param p Polynomial to subtract.
     * @return a new polynomial which is the instance minus {@code p}.
     */
    public FieldPolynomialFunction<T> subtract(final FieldPolynomialFunction<T> p) {
        // identify the lowest degree polynomial
        int lowLength  = FastMath.min(coefficients.length, p.coefficients.length);
        int highLength = FastMath.max(coefficients.length, p.coefficients.length);

        // build the coefficients array
        T[] newCoefficients = MathArrays.buildArray(getField(), highLength);
        for (int i = 0; i < lowLength; ++i) {
            newCoefficients[i] = coefficients[i].subtract(p.coefficients[i]);
        }
        if (coefficients.length < p.coefficients.length) {
            for (int i = lowLength; i < highLength; ++i) {
                newCoefficients[i] = p.coefficients[i].negate();
            }
        } else {
            System.arraycopy(coefficients, lowLength, newCoefficients, lowLength,
                             highLength - lowLength);
        }

        return new FieldPolynomialFunction<>(newCoefficients);
    }

    /**
     * Negate the instance.
     *
     * @return a new polynomial with all coefficients negated
     */
    public FieldPolynomialFunction<T> negate() {
        final T[] newCoefficients = MathArrays.buildArray(getField(), coefficients.length);
        for (int i = 0; i < coefficients.length; ++i) {
            newCoefficients[i] = coefficients[i].negate();
        }
        return new FieldPolynomialFunction<>(newCoefficients);
    }

    /**
     * Multiply the instance by a polynomial.
     *
     * @param p Polynomial to multiply by.
     * @return a new polynomial equal to this times {@code p}
     */
    public FieldPolynomialFunction<T> multiply(final FieldPolynomialFunction<T> p) {
        final Field<T> field = getField();
        final T[] newCoefficients = MathArrays.buildArray(field, coefficients.length + p.coefficients.length - 1);

        for (int i = 0; i < newCoefficients.length; ++i) {
            newCoefficients[i] = field.getZero();
            for (int j = FastMath.max(0, i + 1 - p.coefficients.length);
                 j < FastMath.min(coefficients.length, i + 1);
                 ++j) {
                newCoefficients[i] = newCoefficients[i].add(coefficients[j].multiply(p.coefficients[i-j]));
            }
        }

        return new FieldPolynomialFunction<>(newCoefficients);
    }

    /**
     * Returns the coefficients of the derivative of the polynomial with the given coefficients.
     *
     * @param coefficients Coefficients of the polynomial to differentiate.
     * @param <T> the type of the field elements
     * @return the coefficients of the derivative or {@code null} if coefficients has length 1.
     * @throws MathIllegalArgumentException if {@code coefficients} is empty.
     * @throws NullArgumentException if {@code coefficients} is {@code null}.
     */
    protected static <T extends RealFieldElement<T>> T[] differentiate(T[] coefficients)
        throws MathIllegalArgumentException, NullArgumentException {
        MathUtils.checkNotNull(coefficients);
        int n = coefficients.length;
        if (n == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
        }
        final Field<T> field = coefficients[0].getField();
        final T[] result = MathArrays.buildArray(field, FastMath.max(1, n - 1));
        if (n == 1) {
            result[0] = field.getZero();
        } else {
            for (int i = n - 1; i > 0; i--) {
                result[i - 1] = coefficients[i].multiply(i);
            }
        }
        return result;
    }

    /**
     * Returns an anti-derivative of this polynomial, with 0 constant term.
     *
     * @return a polynomial whose derivative has the same coefficients as this polynomial
     */
    public FieldPolynomialFunction<T> antiDerivative() {
        final Field<T> field = getField();
        final int d = degree();
        final T[] anti = MathArrays.buildArray(field, d + 2);
        anti[0] = field.getZero();
        for (int i = 1; i <= d + 1; i++) {
            anti[i] = coefficients[i - 1].multiply(1.0 / i);
        }
        return new FieldPolynomialFunction<>(anti);
    }

    /**
     * Returns the definite integral of this polymomial over the given interval.
     * <p>
     * [lower, upper] must describe a finite interval (neither can be infinite
     * and lower must be less than or equal to upper).
     *
     * @param lower lower bound for the integration
     * @param upper upper bound for the integration
     * @return the integral of this polymomial over the given interval
     * @throws MathIllegalArgumentException if the bounds do not describe a finite interval
     */
    public T integrate(final double lower, final double upper) {
        final T zero = getField().getZero();
        return integrate(zero.add(lower), zero.add(upper));
    }

    /**
     * Returns the definite integral of this polymomial over the given interval.
     * <p>
     * [lower, upper] must describe a finite interval (neither can be infinite
     * and lower must be less than or equal to upper).
     *
     * @param lower lower bound for the integration
     * @param upper upper bound for the integration
     * @return the integral of this polymomial over the given interval
     * @throws MathIllegalArgumentException if the bounds do not describe a finite interval
     */
    public T integrate(final T lower, final T upper) {
        if (Double.isInfinite(lower.getReal()) || Double.isInfinite(upper.getReal())) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INFINITE_BOUND);
        }
        if (lower.getReal() > upper.getReal()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.LOWER_BOUND_NOT_BELOW_UPPER_BOUND);
        }
        final FieldPolynomialFunction<T> anti = antiDerivative();
        return anti.value(upper).subtract(anti.value(lower));
    }

    /**
     * Returns the derivative as a {@link FieldPolynomialFunction}.
     *
     * @return the derivative polynomial.
     */
    public FieldPolynomialFunction<T> polynomialDerivative() {
        return new FieldPolynomialFunction<>(differentiate(coefficients));
    }

}
