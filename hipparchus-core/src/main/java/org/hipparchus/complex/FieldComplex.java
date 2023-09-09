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
package org.hipparchus.complex;

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.FieldSinhCosh;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;

/**
 * Representation of a Complex number, i.e. a number which has both a
 * real and imaginary part.
 * <p>
 * Implementations of arithmetic operations handle {@code NaN} and
 * infinite values according to the rules for {@link java.lang.Double}, i.e.
 * {@link #equals} is an equivalence relation for all instances that have
 * a {@code NaN} in either real or imaginary part, e.g. the following are
 * considered equal:
 * <ul>
 *  <li>{@code 1 + NaNi}</li>
 *  <li>{@code NaN + i}</li>
 *  <li>{@code NaN + NaNi}</li>
 * </ul>
 * <p>
 * Note that this contradicts the IEEE-754 standard for floating
 * point numbers (according to which the test {@code x == x} must fail if
 * {@code x} is {@code NaN}). The method
 * {@link org.hipparchus.util.Precision#equals(double,double,int)
 * equals for primitive double} in {@link org.hipparchus.util.Precision}
 * conforms with IEEE-754 while this class conforms with the standard behavior
 * for Java object types.
 * @param <T> the type of the field elements
 * @since 2.0
 */
public class FieldComplex<T extends CalculusFieldElement<T>> implements CalculusFieldElement<FieldComplex<T>>  {

    /** A real number representing log(10). */
    private static final double LOG10 = 2.302585092994045684;

    /** The imaginary part. */
    private final T imaginary;

    /** The real part. */
    private final T real;

    /** Record whether this complex number is equal to NaN. */
    private final transient boolean isNaN;

    /** Record whether this complex number is infinite. */
    private final transient boolean isInfinite;

    /**
     * Create a complex number given only the real part.
     *
     * @param real Real part.
     */
    public FieldComplex(T real) {
        this(real, real.getField().getZero());
    }

    /**
     * Create a complex number given the real and imaginary parts.
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     */
    public FieldComplex(T real, T imaginary) {
        this.real = real;
        this.imaginary = imaginary;

        isNaN = real.isNaN() || imaginary.isNaN();
        isInfinite = !isNaN &&
            (real.isInfinite() || imaginary.isInfinite());
    }

    /** Get the square root of -1.
     * @param field field the complex components belong to
     * @return number representing "0.0 + 1.0i"
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> getI(final Field<T> field) {
        return new FieldComplex<>(field.getZero(), field.getOne());
    }

    /** Get the square root of -1.
     * @param field field the complex components belong to
     * @return number representing "0.0 _ 1.0i"
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> getMinusI(final Field<T> field) {
        return new FieldComplex<>(field.getZero(), field.getOne().negate());
    }

    /** Get a complex number representing "NaN + NaNi".
     * @param field field the complex components belong to
     * @return complex number representing "NaN + NaNi"
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> getNaN(final Field<T> field) {
        return new FieldComplex<>(field.getZero().add(Double.NaN), field.getZero().add(Double.NaN));
    }

    /** Get a complex number representing "+INF + INFi".
     * @param field field the complex components belong to
     * @return complex number representing "+INF + INFi"
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> getInf(final Field<T> field) {
        return new FieldComplex<>(field.getZero().add(Double.POSITIVE_INFINITY), field.getZero().add(Double.POSITIVE_INFINITY));
    }

    /** Get a complex number representing "1.0 + 0.0i".
     * @param field field the complex components belong to
     * @return complex number representing "1.0 + 0.0i"
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> getOne(final Field<T> field) {
        return new FieldComplex<>(field.getOne(), field.getZero());
    }

    /** Get a complex number representing "-1.0 + 0.0i".
     * @param field field the complex components belong to
     * @return complex number representing "-1.0 + 0.0i"
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> getMinusOne(final Field<T> field) {
        return new FieldComplex<>(field.getOne().negate(), field.getZero());
    }

    /** Get a complex number representing "0.0 + 0.0i".
     * @param field field the complex components belong to
     * @return complex number representing "0.0 + 0.0i
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> getZero(final Field<T> field) {
        return new FieldComplex<>(field.getZero(), field.getZero());
    }

    /** Get a complex number representing "π + 0.0i".
     * @param field field the complex components belong to
     * @return complex number representing "π + 0.0i
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> getPi(final Field<T> field) {
        return new FieldComplex<>(field.getZero().getPi(), field.getZero());
    }

    /**
     * Return the absolute value of this complex number.
     * Returns {@code NaN} if either real or imaginary part is {@code NaN}
     * and {@code Double.POSITIVE_INFINITY} if neither part is {@code NaN},
     * but at least one part is infinite.
     *
     * @return the absolute value.
     */
    @Override
    public FieldComplex<T> abs() {
        // we check NaN here because FastMath.hypot checks it after infinity
        return isNaN ? getNaN(getPartsField()) : createComplex(FastMath.hypot(real, imaginary), getPartsField().getZero());
    }

    /**
     * Returns a {@code Complex} whose value is
     * {@code (this + addend)}.
     * Uses the definitional formula
     * <p>
     *   {@code (a + bi) + (c + di) = (a+c) + (b+d)i}
     * </p>
     * If either {@code this} or {@code addend} has a {@code NaN} value in
     * either part, {@link #getNaN(Field)} is returned; otherwise {@code Infinite}
     * and {@code NaN} values are returned in the parts of the result
     * according to the rules for {@link java.lang.Double} arithmetic.
     *
     * @param  addend Value to be added to this {@code Complex}.
     * @return {@code this + addend}.
     * @throws NullArgumentException if {@code addend} is {@code null}.
     */
    @Override
    public FieldComplex<T> add(FieldComplex<T> addend) throws NullArgumentException {
        MathUtils.checkNotNull(addend);
        if (isNaN || addend.isNaN) {
            return getNaN(getPartsField());
        }

        return createComplex(real.add(addend.getRealPart()),
                             imaginary.add(addend.getImaginaryPart()));
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this + addend)},
     * with {@code addend} interpreted as a real number.
     *
     * @param addend Value to be added to this {@code Complex}.
     * @return {@code this + addend}.
     * @see #add(FieldComplex)
     */
    public FieldComplex<T> add(T addend) {
        if (isNaN || addend.isNaN()) {
            return getNaN(getPartsField());
        }

        return createComplex(real.add(addend), imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this + addend)},
     * with {@code addend} interpreted as a real number.
     *
     * @param addend Value to be added to this {@code Complex}.
     * @return {@code this + addend}.
     * @see #add(FieldComplex)
     */
    @Override
    public FieldComplex<T> add(double addend) {
        if (isNaN || Double.isNaN(addend)) {
            return getNaN(getPartsField());
        }

        return createComplex(real.add(addend), imaginary);
    }

    /**
     * Returns the conjugate of this complex number.
     * The conjugate of {@code a + bi} is {@code a - bi}.
     * <p>
     * {@link #getNaN(Field)} is returned if either the real or imaginary
     * part of this Complex number equals {@code Double.NaN}.
     * </p><p>
     * If the imaginary part is infinite, and the real part is not
     * {@code NaN}, the returned value has infinite imaginary part
     * of the opposite sign, e.g. the conjugate of
     * {@code 1 + POSITIVE_INFINITY i} is {@code 1 - NEGATIVE_INFINITY i}.
     * </p>
     * @return the conjugate of this Complex object.
     */
    public FieldComplex<T> conjugate() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        return createComplex(real, imaginary.negate());
    }

    /**
     * Returns a {@code Complex} whose value is
     * {@code (this / divisor)}.
     * Implements the definitional formula
     * <pre>
     *  <code>
     *    a + bi          ac + bd + (bc - ad)i
     *    ----------- = -------------------------
     *    c + di         c<sup>2</sup> + d<sup>2</sup>
     *  </code>
     * </pre>
     * but uses
     * <a href="http://doi.acm.org/10.1145/1039813.1039814">
     * prescaling of operands</a> to limit the effects of overflows and
     * underflows in the computation.
     * <p>
     * {@code Infinite} and {@code NaN} values are handled according to the
     * following rules, applied in the order presented:
     * <ul>
     *  <li>If either {@code this} or {@code divisor} has a {@code NaN} value
     *   in either part, {@link #getNaN(Field)} is returned.
     *  </li>
     *  <li>If {@code divisor} equals {@link #getZero(Field)}, {@link #getNaN(Field)} is returned.
     *  </li>
     *  <li>If {@code this} and {@code divisor} are both infinite,
     *   {@link #getNaN(Field)} is returned.
     *  </li>
     *  <li>If {@code this} is finite (i.e., has no {@code Infinite} or
     *   {@code NaN} parts) and {@code divisor} is infinite (one or both parts
     *   infinite), {@link #getZero(Field)} is returned.
     *  </li>
     *  <li>If {@code this} is infinite and {@code divisor} is finite,
     *   {@code NaN} values are returned in the parts of the result if the
     *   {@link java.lang.Double} rules applied to the definitional formula
     *   force {@code NaN} results.
     *  </li>
     * </ul>
     *
     * @param divisor Value by which this {@code Complex} is to be divided.
     * @return {@code this / divisor}.
     * @throws NullArgumentException if {@code divisor} is {@code null}.
     */
    @Override
    public FieldComplex<T> divide(FieldComplex<T> divisor)
        throws NullArgumentException {
        MathUtils.checkNotNull(divisor);
        if (isNaN || divisor.isNaN) {
            return getNaN(getPartsField());
        }

        final T c = divisor.getRealPart();
        final T d = divisor.getImaginaryPart();
        if (c.isZero() && d.isZero()) {
            return getNaN(getPartsField());
        }

        if (divisor.isInfinite() && !isInfinite()) {
            return getZero(getPartsField());
        }

        if (FastMath.abs(c).getReal() < FastMath.abs(d).getReal()) {
            T q = c.divide(d);
            T invDen = c.multiply(q).add(d).reciprocal();
            return createComplex(real.multiply(q).add(imaginary).multiply(invDen),
                                 imaginary.multiply(q).subtract(real).multiply(invDen));
        } else {
            T q = d.divide(c);
            T invDen = d.multiply(q).add(c).reciprocal();
            return createComplex(imaginary.multiply(q).add(real).multiply(invDen),
                                 imaginary.subtract(real.multiply(q)).multiply(invDen));
        }
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this / divisor)},
     * with {@code divisor} interpreted as a real number.
     *
     * @param  divisor Value by which this {@code Complex} is to be divided.
     * @return {@code this / divisor}.
     * @see #divide(FieldComplex)
     */
    public FieldComplex<T> divide(T divisor) {
        if (isNaN || divisor.isNaN()) {
            return getNaN(getPartsField());
        }
        if (divisor.isZero()) {
            return getNaN(getPartsField());
        }
        if (divisor.isInfinite()) {
            return !isInfinite() ? getZero(getPartsField()) : getNaN(getPartsField());
        }
        return createComplex(real.divide(divisor), imaginary.divide(divisor));
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this / divisor)},
     * with {@code divisor} interpreted as a real number.
     *
     * @param  divisor Value by which this {@code Complex} is to be divided.
     * @return {@code this / divisor}.
     * @see #divide(FieldComplex)
     */
    @Override
    public FieldComplex<T> divide(double divisor) {
        if (isNaN || Double.isNaN(divisor)) {
            return getNaN(getPartsField());
        }
        if (divisor == 0.0) {
            return getNaN(getPartsField());
        }
        if (Double.isInfinite(divisor)) {
            return !isInfinite() ? getZero(getPartsField()) : getNaN(getPartsField());
        }
        return createComplex(real.divide(divisor), imaginary.divide(divisor));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> reciprocal() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        if (real.isZero() && imaginary.isZero()) {
            return getInf(getPartsField());
        }

        if (isInfinite) {
            return getZero(getPartsField());
        }

        if (FastMath.abs(real).getReal() < FastMath.abs(imaginary).getReal()) {
            T q = real.divide(imaginary);
            T scale = real.multiply(q).add(imaginary).reciprocal();
            return createComplex(scale.multiply(q), scale.negate());
        } else {
            T q = imaginary.divide(real);
            T scale = imaginary.multiply(q).add(real).reciprocal();
            return createComplex(scale, scale.negate().multiply(q));
        }
    }

    /**
     * Test for equality with another object.
     * If both the real and imaginary parts of two complex numbers
     * are exactly the same, and neither is {@code Double.NaN}, the two
     * Complex objects are considered to be equal.
     * The behavior is the same as for JDK's {@link Double#equals(Object)
     * Double}:
     * <ul>
     *  <li>All {@code NaN} values are considered to be equal,
     *   i.e, if either (or both) real and imaginary parts of the complex
     *   number are equal to {@code Double.NaN}, the complex number is equal
     *   to {@code NaN}.
     *  </li>
     *  <li>
     *   Instances constructed with different representations of zero (i.e.
     *   either "0" or "-0") are <em>not</em> considered to be equal.
     *  </li>
     * </ul>
     *
     * @param other Object to test for equality with this instance.
     * @return {@code true} if the objects are equal, {@code false} if object
     * is {@code null}, not an instance of {@code Complex}, or not equal to
     * this instance.
     */
    @Override
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (other instanceof FieldComplex){
            @SuppressWarnings("unchecked")
            FieldComplex<T> c = (FieldComplex<T>) other;
            if (c.isNaN) {
                return isNaN;
            } else {
                return real.equals(c.real) && imaginary.equals(c.imaginary);
            }
        }
        return false;
    }

    /**
     * Test for the floating-point equality between Complex objects.
     * It returns {@code true} if both arguments are equal or within the
     * range of allowed error (inclusive).
     *
     * @param x First value (cannot be {@code null}).
     * @param y Second value (cannot be {@code null}).
     * @param maxUlps {@code (maxUlps - 1)} is the number of floating point
     * values between the real (resp. imaginary) parts of {@code x} and
     * {@code y}.
     * @param <T> the type of the field elements
     * @return {@code true} if there are fewer than {@code maxUlps} floating
     * point values between the real (resp. imaginary) parts of {@code x}
     * and {@code y}.
     *
     * @see Precision#equals(double,double,int)
     */
    public static <T extends CalculusFieldElement<T>>boolean equals(FieldComplex<T> x, FieldComplex<T> y, int maxUlps) {
        return Precision.equals(x.real.getReal(), y.real.getReal(), maxUlps) &&
               Precision.equals(x.imaginary.getReal(), y.imaginary.getReal(), maxUlps);
    }

    /**
     * Returns {@code true} iff the values are equal as defined by
     * {@link #equals(FieldComplex,FieldComplex,int) equals(x, y, 1)}.
     *
     * @param x First value (cannot be {@code null}).
     * @param y Second value (cannot be {@code null}).
     * @param <T> the type of the field elements
     * @return {@code true} if the values are equal.
     */
    public static <T extends CalculusFieldElement<T>>boolean equals(FieldComplex<T> x, FieldComplex<T> y) {
        return equals(x, y, 1);
    }

    /**
     * Returns {@code true} if, both for the real part and for the imaginary
     * part, there is no T value strictly between the arguments or the
     * difference between them is within the range of allowed error
     * (inclusive).  Returns {@code false} if either of the arguments is NaN.
     *
     * @param x First value (cannot be {@code null}).
     * @param y Second value (cannot be {@code null}).
     * @param eps Amount of allowed absolute error.
     * @param <T> the type of the field elements
     * @return {@code true} if the values are two adjacent floating point
     * numbers or they are within range of each other.
     *
     * @see Precision#equals(double,double,double)
     */
    public static <T extends CalculusFieldElement<T>>boolean equals(FieldComplex<T> x, FieldComplex<T> y,
                                                                    double eps) {
        return Precision.equals(x.real.getReal(), y.real.getReal(), eps) &&
               Precision.equals(x.imaginary.getReal(), y.imaginary.getReal(), eps);
    }

    /**
     * Returns {@code true} if, both for the real part and for the imaginary
     * part, there is no T value strictly between the arguments or the
     * relative difference between them is smaller or equal to the given
     * tolerance. Returns {@code false} if either of the arguments is NaN.
     *
     * @param x First value (cannot be {@code null}).
     * @param y Second value (cannot be {@code null}).
     * @param eps Amount of allowed relative error.
     * @param <T> the type of the field elements
     * @return {@code true} if the values are two adjacent floating point
     * numbers or they are within range of each other.
     *
     * @see Precision#equalsWithRelativeTolerance(double,double,double)
     */
    public static <T extends CalculusFieldElement<T>>boolean equalsWithRelativeTolerance(FieldComplex<T> x,
                                                                                         FieldComplex<T> y,
                                                                                         double eps) {
        return Precision.equalsWithRelativeTolerance(x.real.getReal(), y.real.getReal(), eps) &&
               Precision.equalsWithRelativeTolerance(x.imaginary.getReal(), y.imaginary.getReal(), eps);
    }

    /**
     * Get a hashCode for the complex number.
     * Any {@code Double.NaN} value in real or imaginary part produces
     * the same hash code {@code 7}.
     *
     * @return a hash code value for this object.
     */
    @Override
    public int hashCode() {
        if (isNaN) {
            return 7;
        }
        return 37 * (17 * imaginary.hashCode() + real.hashCode());
    }

    /** {@inheritDoc}
     * <p>
     * This implementation considers +0.0 and -0.0 to be equal for both
     * real and imaginary components.
     * </p>
     */
    @Override
    public boolean isZero() {
        return real.isZero() && imaginary.isZero();
    }

    /**
     * Access the imaginary part.
     *
     * @return the imaginary part.
     */
    public T getImaginary() {
        return imaginary;
    }

    /**
     * Access the imaginary part.
     *
     * @return the imaginary part.
     */
    public T getImaginaryPart() {
        return imaginary;
    }

    /**
     * Access the real part.
     *
     * @return the real part.
     */
    @Override
    public double getReal() {
        return real.getReal();
    }

    /**
     * Access the real part.
     *
     * @return the real part.
     */
    public T getRealPart() {
        return real;
    }

    /**
     * Checks whether either or both parts of this complex number is
     * {@code NaN}.
     *
     * @return true if either or both parts of this complex number is
     * {@code NaN}; false otherwise.
     */
    @Override
    public boolean isNaN() {
        return isNaN;
    }

    /** Check whether the instance is real (i.e. imaginary part is zero).
     * @return true if imaginary part is zero
      */
    public boolean isReal() {
        return imaginary.isZero();
    }

    /** Check whether the instance is an integer (i.e. imaginary part is zero and real part has no fractional part).
     * @return true if imaginary part is zero and real part has no fractional part
     */
    public boolean isMathematicalInteger() {
        return isReal() && Precision.isMathematicalInteger(real.getReal());
    }

    /**
     * Checks whether either the real or imaginary part of this complex number
     * takes an infinite value (either {@code Double.POSITIVE_INFINITY} or
     * {@code Double.NEGATIVE_INFINITY}) and neither part
     * is {@code NaN}.
     *
     * @return true if one or both parts of this complex number are infinite
     * and neither part is {@code NaN}.
     */
    @Override
    public boolean isInfinite() {
        return isInfinite;
    }

    /**
     * Returns a {@code Complex} whose value is {@code this * factor}.
     * Implements preliminary checks for {@code NaN} and infinity followed by
     * the definitional formula:
     * <p>
     *   {@code (a + bi)(c + di) = (ac - bd) + (ad + bc)i}
     * </p>
     * Returns {@link #getNaN(Field)} if either {@code this} or {@code factor} has one or
     * more {@code NaN} parts.
     * <p>
     * Returns {@link #getInf(Field)} if neither {@code this} nor {@code factor} has one
     * or more {@code NaN} parts and if either {@code this} or {@code factor}
     * has one or more infinite parts (same result is returned regardless of
     * the sign of the components).
     * </p><p>
     * Returns finite values in components of the result per the definitional
     * formula in all remaining cases.</p>
     *
     * @param  factor value to be multiplied by this {@code Complex}.
     * @return {@code this * factor}.
     * @throws NullArgumentException if {@code factor} is {@code null}.
     */
    @Override
    public FieldComplex<T> multiply(FieldComplex<T> factor)
        throws NullArgumentException {
        MathUtils.checkNotNull(factor);
        if (isNaN || factor.isNaN) {
            return getNaN(getPartsField());
        }
        if (real.isInfinite() ||
            imaginary.isInfinite() ||
            factor.real.isInfinite() ||
            factor.imaginary.isInfinite()) {
            // we don't use isInfinite() to avoid testing for NaN again
            return getInf(getPartsField());
        }
        return createComplex(real.linearCombination(real, factor.real, imaginary.negate(), factor.imaginary),
                             real.linearCombination(real, factor.imaginary, imaginary, factor.real));
    }

    /**
     * Returns a {@code Complex} whose value is {@code this * factor}, with {@code factor}
     * interpreted as a integer number.
     *
     * @param  factor value to be multiplied by this {@code Complex}.
     * @return {@code this * factor}.
     * @see #multiply(FieldComplex)
     */
    @Override
    public FieldComplex<T> multiply(final int factor) {
        if (isNaN) {
            return getNaN(getPartsField());
        }
        if (real.isInfinite() || imaginary.isInfinite()) {
            return getInf(getPartsField());
        }
        return createComplex(real.multiply(factor), imaginary.multiply(factor));
    }

    /**
     * Returns a {@code Complex} whose value is {@code this * factor}, with {@code factor}
     * interpreted as a real number.
     *
     * @param  factor value to be multiplied by this {@code Complex}.
     * @return {@code this * factor}.
     * @see #multiply(FieldComplex)
     */
    @Override
    public FieldComplex<T> multiply(double factor) {
        if (isNaN || Double.isNaN(factor)) {
            return getNaN(getPartsField());
        }
        if (real.isInfinite() ||
            imaginary.isInfinite() ||
            Double.isInfinite(factor)) {
            // we don't use isInfinite() to avoid testing for NaN again
            return getInf(getPartsField());
        }
        return createComplex(real.multiply(factor), imaginary.multiply(factor));
    }

    /**
     * Returns a {@code Complex} whose value is {@code this * factor}, with {@code factor}
     * interpreted as a real number.
     *
     * @param  factor value to be multiplied by this {@code Complex}.
     * @return {@code this * factor}.
     * @see #multiply(FieldComplex)
     */
    public FieldComplex<T> multiply(T factor) {
        if (isNaN || factor.isNaN()) {
            return getNaN(getPartsField());
        }
        if (real.isInfinite() ||
            imaginary.isInfinite() ||
            factor.isInfinite()) {
            // we don't use isInfinite() to avoid testing for NaN again
            return getInf(getPartsField());
        }
        return createComplex(real.multiply(factor), imaginary.multiply(factor));
    }

    /** Compute this * i.
     * @return this * i
     * @since 2.0
     */
    public FieldComplex<T> multiplyPlusI() {
        return createComplex(imaginary.negate(), real);
    }

    /** Compute this *- -i.
     * @return this * i
     * @since 2.0
     */
    public FieldComplex<T> multiplyMinusI() {
        return createComplex(imaginary, real.negate());
    }

    /**
     * Returns a {@code Complex} whose value is {@code (-this)}.
     * Returns {@code NaN} if either real or imaginary
     * part of this Complex number is {@code Double.NaN}.
     *
     * @return {@code -this}.
     */
    @Override
    public FieldComplex<T> negate() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        return createComplex(real.negate(), imaginary.negate());
    }

    /**
     * Returns a {@code Complex} whose value is
     * {@code (this - subtrahend)}.
     * Uses the definitional formula
     * <p>
     *  {@code (a + bi) - (c + di) = (a-c) + (b-d)i}
     * </p>
     * If either {@code this} or {@code subtrahend} has a {@code NaN]} value in either part,
     * {@link #getNaN(Field)} is returned; otherwise infinite and {@code NaN} values are
     * returned in the parts of the result according to the rules for
     * {@link java.lang.Double} arithmetic.
     *
     * @param  subtrahend value to be subtracted from this {@code Complex}.
     * @return {@code this - subtrahend}.
     * @throws NullArgumentException if {@code subtrahend} is {@code null}.
     */
    @Override
    public FieldComplex<T> subtract(FieldComplex<T> subtrahend)
        throws NullArgumentException {
        MathUtils.checkNotNull(subtrahend);
        if (isNaN || subtrahend.isNaN) {
            return getNaN(getPartsField());
        }

        return createComplex(real.subtract(subtrahend.getRealPart()),
                             imaginary.subtract(subtrahend.getImaginaryPart()));
    }

    /**
     * Returns a {@code Complex} whose value is
     * {@code (this - subtrahend)}.
     *
     * @param  subtrahend value to be subtracted from this {@code Complex}.
     * @return {@code this - subtrahend}.
     * @see #subtract(FieldComplex)
     */
    @Override
    public FieldComplex<T> subtract(double subtrahend) {
        if (isNaN || Double.isNaN(subtrahend)) {
            return getNaN(getPartsField());
        }
        return createComplex(real.subtract(subtrahend), imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is
     * {@code (this - subtrahend)}.
     *
     * @param  subtrahend value to be subtracted from this {@code Complex}.
     * @return {@code this - subtrahend}.
     * @see #subtract(FieldComplex)
     */
    public FieldComplex<T> subtract(T subtrahend) {
        if (isNaN || subtrahend.isNaN()) {
            return getNaN(getPartsField());
        }
        return createComplex(real.subtract(subtrahend), imaginary);
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/InverseCosine.html" TARGET="_top">
     * inverse cosine</a> of this complex number.
     * Implements the formula:
     * <p>
     *  {@code acos(z) = -i (log(z + i (sqrt(1 - z<sup>2</sup>))))}
     * </p>
     * Returns {@link #getNaN(Field)} if either real or imaginary part of the
     * input argument is {@code NaN} or infinite.
     *
     * @return the inverse cosine of this complex number.
     */
    @Override
    public FieldComplex<T> acos() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        return this.add(this.sqrt1z().multiplyPlusI()).log().multiplyMinusI();
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/InverseSine.html" TARGET="_top">
     * inverse sine</a> of this complex number.
     * Implements the formula:
     * <p>
     *  {@code asin(z) = -i (log(sqrt(1 - z<sup>2</sup>) + iz))}
     * </p><p>
     * Returns {@link #getNaN(Field)} if either real or imaginary part of the
     * input argument is {@code NaN} or infinite.</p>
     *
     * @return the inverse sine of this complex number.
     */
    @Override
    public FieldComplex<T> asin() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        return sqrt1z().add(this.multiplyPlusI()).log().multiplyMinusI();
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/InverseTangent.html" TARGET="_top">
     * inverse tangent</a> of this complex number.
     * Implements the formula:
     * <p>
     * {@code atan(z) = (i/2) log((1 - iz)/(1 + iz))}
     * </p><p>
     * Returns {@link #getNaN(Field)} if either real or imaginary part of the
     * input argument is {@code NaN} or infinite.</p>
     *
     * @return the inverse tangent of this complex number
     */
    @Override
    public FieldComplex<T> atan() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        final T one = getPartsField().getOne();
        if (real.isZero()) {

            // singularity at ±i
            if (imaginary.multiply(imaginary).subtract(one).isZero()) {
                return getNaN(getPartsField());
            }

            // branch cut on imaginary axis
            final T zero = getPartsField().getZero();
            final FieldComplex<T> tmp = createComplex(one.add(imaginary).divide(one.subtract(imaginary)), zero).
                                        log().multiplyPlusI().multiply(0.5);
            return createComplex(FastMath.copySign(tmp.real, real), tmp.imaginary);

        } else {
            // regular formula
            final FieldComplex<T> n = createComplex(one.add(imaginary), real.negate());
            final FieldComplex<T> d = createComplex(one.subtract(imaginary),  real);
            return n.divide(d).log().multiplyPlusI().multiply(0.5);
        }

    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/Cosine.html" TARGET="_top">
     * cosine</a> of this complex number.
     * Implements the formula:
     * <p>
     *  {@code cos(a + bi) = cos(a)cosh(b) - sin(a)sinh(b)i}
     * </p><p>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos},
     * {@link FastMath#cosh} and {@link FastMath#sinh}.
     * </p><p>
     * Returns {@link #getNaN(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p><p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.</p>
     * <pre>
     *  Examples:
     *  <code>
     *   cos(1 &plusmn; INFINITY i) = 1 \u2213 INFINITY i
     *   cos(&plusmn;INFINITY + i) = NaN + NaN i
     *   cos(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return the cosine of this complex number.
     */
    @Override
    public FieldComplex<T> cos() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        final FieldSinCos<T>   scr  = FastMath.sinCos(real);
        final FieldSinhCosh<T> schi = FastMath.sinhCosh(imaginary);
        return createComplex(scr.cos().multiply(schi.cosh()), scr.sin().negate().multiply(schi.sinh()));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/HyperbolicCosine.html" TARGET="_top">
     * hyperbolic cosine</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   cosh(a + bi) = cosh(a)cos(b) + sinh(a)sin(b)i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos},
     * {@link FastMath#cosh} and {@link FastMath#sinh}.
     * <p>
     * Returns {@link #getNaN(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   cosh(1 &plusmn; INFINITY i) = NaN + NaN i
     *   cosh(&plusmn;INFINITY + i) = INFINITY &plusmn; INFINITY i
     *   cosh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return the hyperbolic cosine of this complex number.
     */
    @Override
    public FieldComplex<T> cosh() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        final FieldSinhCosh<T> schr = FastMath.sinhCosh(real);
        final FieldSinCos<T>   sci  = FastMath.sinCos(imaginary);
        return createComplex(schr.cosh().multiply(sci.cos()), schr.sinh().multiply(sci.sin()));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/ExponentialFunction.html" TARGET="_top">
     * exponential function</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   exp(a + bi) = exp(a)cos(b) + exp(a)sin(b)i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#exp}, {@link FastMath#cos}, and
     * {@link FastMath#sin}.
     * <p>
     * Returns {@link #getNaN(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   exp(1 &plusmn; INFINITY i) = NaN + NaN i
     *   exp(INFINITY + i) = INFINITY + INFINITY i
     *   exp(-INFINITY + i) = 0 + 0i
     *   exp(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return <code><i>e</i><sup>this</sup></code>.
     */
    @Override
    public FieldComplex<T> exp() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        final T              expReal = FastMath.exp(real);
        final FieldSinCos<T> sc      = FastMath.sinCos(imaginary);
        return createComplex(expReal.multiply(sc.cos()), expReal.multiply(sc.sin()));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> expm1() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        final T              expm1Real = FastMath.expm1(real);
        final FieldSinCos<T> sc        = FastMath.sinCos(imaginary);
        return createComplex(expm1Real.multiply(sc.cos()), expm1Real.multiply(sc.sin()));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/NaturalLogarithm.html" TARGET="_top">
     * natural logarithm</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   log(a + bi) = ln(|a + bi|) + arg(a + bi)i
     *  </code>
     * </pre>
     * where ln on the right hand side is {@link FastMath#log},
     * {@code |a + bi|} is the modulus, {@link #abs},  and
     * {@code arg(a + bi) = }{@link FastMath#atan2}(b, a).
     * <p>
     * Returns {@link #getNaN(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * Infinite (or critical) values in real or imaginary parts of the input may
     * result in infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   log(1 &plusmn; INFINITY i) = INFINITY &plusmn; (&pi;/2)i
     *   log(INFINITY + i) = INFINITY + 0i
     *   log(-INFINITY + i) = INFINITY + &pi;i
     *   log(INFINITY &plusmn; INFINITY i) = INFINITY &plusmn; (&pi;/4)i
     *   log(-INFINITY &plusmn; INFINITY i) = INFINITY &plusmn; (3&pi;/4)i
     *   log(0 + 0i) = -INFINITY + 0i
     *  </code>
     * </pre>
     *
     * @return the value <code>ln &nbsp; this</code>, the natural logarithm
     * of {@code this}.
     */
    @Override
    public FieldComplex<T> log() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        return createComplex(FastMath.log(FastMath.hypot(real, imaginary)),
                             FastMath.atan2(imaginary, real));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> log1p() {
        return add(1.0).log();
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> log10() {
        return log().divide(LOG10);
    }

    /**
     * Returns of value of this complex number raised to the power of {@code x}.
     * <p>
     * If {@code x} is a real number whose real part has an integer value, returns {@link #pow(int)},
     * if both {@code this} and {@code x} are real and {@link FastMath#pow(double, double)}
     * with the corresponding real arguments would return a finite number (neither NaN
     * nor infinite), then returns the same value converted to {@code Complex},
     * with the same special cases.
     * In all other cases real cases, implements y<sup>x</sup> = exp(x&middot;log(y)).
     * </p>
     *
     * @param  x exponent to which this {@code Complex} is to be raised.
     * @return <code> this<sup>x</sup></code>.
     * @throws NullArgumentException if x is {@code null}.
     */
    @Override
    public FieldComplex<T> pow(FieldComplex<T> x)
        throws NullArgumentException {

        MathUtils.checkNotNull(x);

        if (x.imaginary.isZero()) {
            final int nx = (int) FastMath.rint(x.real.getReal());
            if (x.real.getReal() == nx) {
                // integer power
                return pow(nx);
            } else if (this.imaginary.isZero()) {
                // check real implementation that handles a bunch of special cases
                final T realPow = FastMath.pow(this.real, x.real);
                if (realPow.isFinite()) {
                    return createComplex(realPow, getPartsField().getZero());
                }
            }
        }

        // generic implementation
        return this.log().multiply(x).exp();

    }


    /**
     * Returns of value of this complex number raised to the power of {@code x}.
     * <p>
     * If {@code x} has an integer value, returns {@link #pow(int)},
     * if {@code this} is real and {@link FastMath#pow(double, double)}
     * with the corresponding real arguments would return a finite number (neither NaN
     * nor infinite), then returns the same value converted to {@code Complex},
     * with the same special cases.
     * In all other cases real cases, implements y<sup>x</sup> = exp(x&middot;log(y)).
     * </p>
     *
     * @param  x exponent to which this {@code Complex} is to be raised.
     * @return <code> this<sup>x</sup></code>.
     */
    public FieldComplex<T> pow(T x) {

        final int nx = (int) FastMath.rint(x.getReal());
        if (x.getReal() == nx) {
            // integer power
            return pow(nx);
        } else if (this.imaginary.isZero()) {
            // check real implementation that handles a bunch of special cases
            final T realPow = FastMath.pow(this.real, x);
            if (realPow.isFinite()) {
                return createComplex(realPow, getPartsField().getZero());
            }
        }

        // generic implementation
        return this.log().multiply(x).exp();

    }

    /**
     * Returns of value of this complex number raised to the power of {@code x}.
     * <p>
     * If {@code x} has an integer value, returns {@link #pow(int)},
     * if {@code this} is real and {@link FastMath#pow(double, double)}
     * with the corresponding real arguments would return a finite number (neither NaN
     * nor infinite), then returns the same value converted to {@code Complex},
     * with the same special cases.
     * In all other cases real cases, implements y<sup>x</sup> = exp(x&middot;log(y)).
     * </p>
     *
     * @param  x exponent to which this {@code Complex} is to be raised.
     * @return <code> this<sup>x</sup></code>.
     */
    @Override
    public FieldComplex<T> pow(double x) {

        final int nx = (int) FastMath.rint(x);
        if (x == nx) {
            // integer power
            return pow(nx);
        } else if (this.imaginary.isZero()) {
            // check real implementation that handles a bunch of special cases
            final T realPow = FastMath.pow(this.real, x);
            if (realPow.isFinite()) {
                return createComplex(realPow, getPartsField().getZero());
            }
        }

        // generic implementation
        return this.log().multiply(x).exp();

    }

     /** {@inheritDoc} */
    @Override
    public FieldComplex<T> pow(final int n) {

        FieldComplex<T> result = getField().getOne();
        final boolean invert;
        int p = n;
        if (p < 0) {
            invert = true;
            p = -p;
        } else {
            invert = false;
        }

        // Exponentiate by successive squaring
        FieldComplex<T> square = this;
        while (p > 0) {
            if ((p & 0x1) > 0) {
                result = result.multiply(square);
            }
            square = square.multiply(square);
            p = p >> 1;
        }

        return invert ? result.reciprocal() : result;

    }

     /**
      * Compute the
     * <a href="http://mathworld.wolfram.com/Sine.html" TARGET="_top">
     * sine</a>
     * of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   sin(a + bi) = sin(a)cosh(b) + cos(a)sinh(b)i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos},
     * {@link FastMath#cosh} and {@link FastMath#sinh}.
     * <p>
     * Returns {@link #getNaN(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p><p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or {@code NaN} values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   sin(1 &plusmn; INFINITY i) = 1 &plusmn; INFINITY i
     *   sin(&plusmn;INFINITY + i) = NaN + NaN i
     *   sin(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return the sine of this complex number.
     */
    @Override
    public FieldComplex<T> sin() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        final FieldSinCos<T>   scr  = FastMath.sinCos(real);
        final FieldSinhCosh<T> schi = FastMath.sinhCosh(imaginary);
        return createComplex(scr.sin().multiply(schi.cosh()), scr.cos().multiply(schi.sinh()));

    }

    /** {@inheritDoc}
     */
    @Override
    public FieldSinCos<FieldComplex<T>> sinCos() {
        if (isNaN) {
            return new FieldSinCos<>(getNaN(getPartsField()), getNaN(getPartsField()));
        }

        final FieldSinCos<T>   scr = FastMath.sinCos(real);
        final FieldSinhCosh<T> schi = FastMath.sinhCosh(imaginary);
        return new FieldSinCos<>(createComplex(scr.sin().multiply(schi.cosh()), scr.cos().multiply(schi.sinh())),
                                 createComplex(scr.cos().multiply(schi.cosh()), scr.sin().negate().multiply(schi.sinh())));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> atan2(FieldComplex<T> x) {

        // compute r = sqrt(x^2+y^2)
        final FieldComplex<T> r = x.multiply(x).add(multiply(this)).sqrt();

        if (x.real.getReal() >= 0) {
            // compute atan2(y, x) = 2 atan(y / (r + x))
            return divide(r.add(x)).atan().multiply(2);
        } else {
            // compute atan2(y, x) = +/- pi - 2 atan(y / (r - x))
            return divide(r.subtract(x)).atan().multiply(-2).add(x.real.getPi());
        }
    }

    /** {@inheritDoc}
     * <p>
     * Branch cuts are on the real axis, below +1.
     * </p>
     */
    @Override
    public FieldComplex<T> acosh() {
        final FieldComplex<T> sqrtPlus  = add(1).sqrt();
        final FieldComplex<T> sqrtMinus = subtract(1).sqrt();
        return add(sqrtPlus.multiply(sqrtMinus)).log();
    }

    /** {@inheritDoc}
     * <p>
     * Branch cuts are on the imaginary axis, above +i and below -i.
     * </p>
     */
    @Override
    public FieldComplex<T> asinh() {
        return add(multiply(this).add(1.0).sqrt()).log();
    }

    /** {@inheritDoc}
     * <p>
     * Branch cuts are on the real axis, above +1 and below -1.
     * </p>
     */
    @Override
    public FieldComplex<T> atanh() {
        final FieldComplex<T> logPlus  = add(1).log();
        final FieldComplex<T> logMinus = createComplex(getPartsField().getOne().subtract(real), imaginary.negate()).log();
        return logPlus.subtract(logMinus).multiply(0.5);
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/HyperbolicSine.html" TARGET="_top">
     * hyperbolic sine</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   sinh(a + bi) = sinh(a)cos(b)) + cosh(a)sin(b)i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos},
     * {@link FastMath#cosh} and {@link FastMath#sinh}.
     * <p>
     * Returns {@link #getNaN(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p><p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   sinh(1 &plusmn; INFINITY i) = NaN + NaN i
     *   sinh(&plusmn;INFINITY + i) = &plusmn; INFINITY + INFINITY i
     *   sinh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return the hyperbolic sine of {@code this}.
     */
    @Override
    public FieldComplex<T> sinh() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        final FieldSinhCosh<T> schr = FastMath.sinhCosh(real);
        final FieldSinCos<T>   sci  = FastMath.sinCos(imaginary);
        return createComplex(schr.sinh().multiply(sci.cos()), schr.cosh().multiply(sci.sin()));
    }

    /** {@inheritDoc}
     */
    @Override
    public FieldSinhCosh<FieldComplex<T>> sinhCosh() {
        if (isNaN) {
            return new FieldSinhCosh<>(getNaN(getPartsField()), getNaN(getPartsField()));
        }

        final FieldSinhCosh<T> schr = FastMath.sinhCosh(real);
        final FieldSinCos<T>   sci  = FastMath.sinCos(imaginary);
        return new FieldSinhCosh<>(createComplex(schr.sinh().multiply(sci.cos()), schr.cosh().multiply(sci.sin())),
                                   createComplex(schr.cosh().multiply(sci.cos()), schr.sinh().multiply(sci.sin())));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/SquareRoot.html" TARGET="_top">
     * square root</a> of this complex number.
     * Implements the following algorithm to compute {@code sqrt(a + bi)}:
     * <ol><li>Let {@code t = sqrt((|a| + |a + bi|) / 2)}</li>
     * <li><pre>if {@code  a ≥ 0} return {@code t + (b/2t)i}
     *  else return {@code |b|/2t + sign(b)t i }</pre></li>
     * </ol>
     * where <ul>
     * <li>{@code |a| = }{@link FastMath#abs(CalculusFieldElement) abs(a)}</li>
     * <li>{@code |a + bi| = }{@link FastMath#hypot(CalculusFieldElement, CalculusFieldElement) hypot(a, b)}</li>
     * <li>{@code sign(b) = }{@link FastMath#copySign(CalculusFieldElement, CalculusFieldElement) copySign(1, b)}
     * </ul>
     * The real part is therefore always nonnegative.
     * <p>
     * Returns {@link #getNaN(Field) NaN} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * <p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * </p>
     * <pre>
     *  Examples:
     *  <code>
     *   sqrt(1 ± ∞ i) = ∞ + NaN i
     *   sqrt(∞ + i) = ∞ + 0i
     *   sqrt(-∞ + i) = 0 + ∞ i
     *   sqrt(∞ ± ∞ i) = ∞ + NaN i
     *   sqrt(-∞ ± ∞ i) = NaN ± ∞ i
     *  </code>
     * </pre>
     *
     * @return the square root of {@code this} with nonnegative real part.
     */
    @Override
    public FieldComplex<T> sqrt() {
        if (isNaN) {
            return getNaN(getPartsField());
        }

        if (isZero()) {
            return getZero(getPartsField());
        }

        T t = FastMath.sqrt((FastMath.abs(real).add(FastMath.hypot(real, imaginary))).multiply(0.5));
        if (real.getReal() >= 0.0) {
            return createComplex(t, imaginary.divide(t.multiply(2)));
        } else {
            return createComplex(FastMath.abs(imaginary).divide(t.multiply(2)),
                                 FastMath.copySign(t, imaginary));
        }
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/SquareRoot.html" TARGET="_top">
     * square root</a> of <code>1 - this<sup>2</sup></code> for this complex
     * number.
     * Computes the result directly as
     * {@code sqrt(ONE.subtract(z.multiply(z)))}.
     * <p>
     * Returns {@link #getNaN(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     *
     * @return the square root of <code>1 - this<sup>2</sup></code>.
     */
    public FieldComplex<T> sqrt1z() {
        final FieldComplex<T> t2 = this.multiply(this);
        return createComplex(getPartsField().getOne().subtract(t2.real), t2.imaginary.negate()).sqrt();
    }

    /** {@inheritDoc}
     * <p>
     * This implementation compute the principal cube root by using a branch cut along real negative axis.
     * </p>
     */
    @Override
    public FieldComplex<T> cbrt() {
        final T              magnitude = FastMath.cbrt(abs().getRealPart());
        final FieldSinCos<T> sc        = FastMath.sinCos(getArgument().divide(3));
        return createComplex(magnitude.multiply(sc.cos()), magnitude.multiply(sc.sin()));
    }

    /** {@inheritDoc}
     * <p>
     * This implementation compute the principal n<sup>th</sup> root by using a branch cut along real negative axis.
     * </p>
     */
    @Override
    public FieldComplex<T> rootN(int n) {
        final T              magnitude = FastMath.pow(abs().getRealPart(), 1.0 / n);
        final FieldSinCos<T> sc        = FastMath.sinCos(getArgument().divide(n));
        return createComplex(magnitude.multiply(sc.cos()), magnitude.multiply(sc.sin()));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/Tangent.html" TARGET="_top">
     * tangent</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   tan(a + bi) = sin(2a)/(cos(2a)+cosh(2b)) + [sinh(2b)/(cos(2a)+cosh(2b))]i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos}, {@link FastMath#cosh} and
     * {@link FastMath#sinh}.
     * <p>
     * Returns {@link #getNaN(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * Infinite (or critical) values in real or imaginary parts of the input may
     * result in infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   tan(a &plusmn; INFINITY i) = 0 &plusmn; i
     *   tan(&plusmn;INFINITY + bi) = NaN + NaN i
     *   tan(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *   tan(&plusmn;&pi;/2 + 0 i) = &plusmn;INFINITY + NaN i
     *  </code>
     * </pre>
     *
     * @return the tangent of {@code this}.
     */
    @Override
    public FieldComplex<T> tan() {
        if (isNaN || real.isInfinite()) {
            return getNaN(getPartsField());
        }
        if (imaginary.getReal() > 20.0) {
            return getI(getPartsField());
        }
        if (imaginary.getReal() < -20.0) {
            return getMinusI(getPartsField());
        }

        final FieldSinCos<T> sc2r = FastMath.sinCos(real.multiply(2));
        T imaginary2 = imaginary.multiply(2);
        T d = sc2r.cos().add(FastMath.cosh(imaginary2));

        return createComplex(sc2r.sin().divide(d), FastMath.sinh(imaginary2).divide(d));

    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/HyperbolicTangent.html" TARGET="_top">
     * hyperbolic tangent</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   tan(a + bi) = sinh(2a)/(cosh(2a)+cos(2b)) + [sin(2b)/(cosh(2a)+cos(2b))]i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos}, {@link FastMath#cosh} and
     * {@link FastMath#sinh}.
     * <p>
     * Returns {@link #getNaN(Field)} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   tanh(a &plusmn; INFINITY i) = NaN + NaN i
     *   tanh(&plusmn;INFINITY + bi) = &plusmn;1 + 0 i
     *   tanh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *   tanh(0 + (&pi;/2)i) = NaN + INFINITY i
     *  </code>
     * </pre>
     *
     * @return the hyperbolic tangent of {@code this}.
     */
    @Override
    public FieldComplex<T> tanh() {
        if (isNaN || imaginary.isInfinite()) {
            return getNaN(getPartsField());
        }
        if (real.getReal() > 20.0) {
            return getOne(getPartsField());
        }
        if (real.getReal() < -20.0) {
            return getMinusOne(getPartsField());
        }
        T real2 = real.multiply(2);
        final FieldSinCos<T> sc2i = FastMath.sinCos(imaginary.multiply(2));
        T d = FastMath.cosh(real2).add(sc2i.cos());

        return createComplex(FastMath.sinh(real2).divide(d), sc2i.sin().divide(d));
    }



    /**
     * Compute the argument of this complex number.
     * The argument is the angle phi between the positive real axis and
     * the point representing this number in the complex plane.
     * The value returned is between -PI (not inclusive)
     * and PI (inclusive), with negative values returned for numbers with
     * negative imaginary parts.
     * <p>
     * If either real or imaginary part (or both) is NaN, NaN is returned.
     * Infinite parts are handled as {@code Math.atan2} handles them,
     * essentially treating finite parts as zero in the presence of an
     * infinite coordinate and returning a multiple of pi/4 depending on
     * the signs of the infinite parts.
     * See the javadoc for {@code Math.atan2} for full details.
     *
     * @return the argument of {@code this}.
     */
    public T getArgument() {
        return FastMath.atan2(getImaginaryPart(), getRealPart());
    }

    /**
     * Computes the n-th roots of this complex number.
     * The nth roots are defined by the formula:
     * <pre>
     *  <code>
     *   z<sub>k</sub> = abs<sup>1/n</sup> (cos(phi + 2&pi;k/n) + i (sin(phi + 2&pi;k/n))
     *  </code>
     * </pre>
     * for <i>{@code k=0, 1, ..., n-1}</i>, where {@code abs} and {@code phi}
     * are respectively the {@link #abs() modulus} and
     * {@link #getArgument() argument} of this complex number.
     * <p>
     * If one or both parts of this complex number is NaN, a list with just
     * one element, {@link #getNaN(Field)} is returned.
     * if neither part is NaN, but at least one part is infinite, the result
     * is a one-element list containing {@link #getInf(Field)}.
     *
     * @param n Degree of root.
     * @return a List of all {@code n}-th roots of {@code this}.
     * @throws MathIllegalArgumentException if {@code n <= 0}.
     */
    public List<FieldComplex<T>> nthRoot(int n) throws MathIllegalArgumentException {

        if (n <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.CANNOT_COMPUTE_NTH_ROOT_FOR_NEGATIVE_N,
                                                   n);
        }

        final List<FieldComplex<T>> result = new ArrayList<>();

        if (isNaN) {
            result.add(getNaN(getPartsField()));
            return result;
        }
        if (isInfinite()) {
            result.add(getInf(getPartsField()));
            return result;
        }

        // nth root of abs -- faster / more accurate to use a solver here?
        final T nthRootOfAbs = FastMath.pow(FastMath.hypot(real, imaginary), 1.0 / n);

        // Compute nth roots of complex number with k = 0, 1, ... n-1
        final T nthPhi = getArgument().divide(n);
        final double slice = 2 * FastMath.PI / n;
        T innerPart = nthPhi;
        for (int k = 0; k < n ; k++) {
            // inner part
            final FieldSinCos<T> scInner = FastMath.sinCos(innerPart);
            final T realPart = nthRootOfAbs.multiply(scInner.cos());
            final T imaginaryPart = nthRootOfAbs.multiply(scInner.sin());
            result.add(createComplex(realPart, imaginaryPart));
            innerPart = innerPart.add(slice);
        }

        return result;
    }

    /**
     * Create a complex number given the real and imaginary parts.
     *
     * @param realPart Real part.
     * @param imaginaryPart Imaginary part.
     * @return a new complex number instance.
     *
     * @see #valueOf(CalculusFieldElement, CalculusFieldElement)
     */
    protected FieldComplex<T> createComplex(final T realPart, final T imaginaryPart) {
        return new FieldComplex<>(realPart, imaginaryPart);
    }

    /**
     * Create a complex number given the real and imaginary parts.
     *
     * @param realPart Real part.
     * @param imaginaryPart Imaginary part.
     * @param <T> the type of the field elements
     * @return a Complex instance.
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T>
        valueOf(T realPart, T imaginaryPart) {
        if (realPart.isNaN() || imaginaryPart.isNaN()) {
            return getNaN(realPart.getField());
        }
        return new FieldComplex<>(realPart, imaginaryPart);
    }

    /**
     * Create a complex number given only the real part.
     *
     * @param realPart Real part.
     * @param <T> the type of the field elements
     * @return a Complex instance.
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T>
        valueOf(T realPart) {
        if (realPart.isNaN()) {
            return getNaN(realPart.getField());
        }
        return new FieldComplex<>(realPart);
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> newInstance(double realPart) {
        return valueOf(getPartsField().getZero().newInstance(realPart));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplexField<T> getField() {
        return FieldComplexField.getField(getPartsField());
    }

    /** Get the {@link Field} the real and imaginary parts belong to.
     * @return {@link Field} the real and imaginary parts belong to
     */
    public Field<T> getPartsField() {
        return real.getField();
    }

    /** {@inheritDoc} */
    @Override
    public String toString() {
        return "(" + real + ", " + imaginary + ")";
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> scalb(int n) {
        return createComplex(FastMath.scalb(real, n), FastMath.scalb(imaginary, n));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> ulp() {
        return createComplex(FastMath.ulp(real), FastMath.ulp(imaginary));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> hypot(FieldComplex<T> y) {
        if (isInfinite() || y.isInfinite()) {
            return getInf(getPartsField());
        } else if (isNaN() || y.isNaN()) {
            return getNaN(getPartsField());
        } else {
            return multiply(this).add(y.multiply(y)).sqrt();
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> linearCombination(final FieldComplex<T>[] a, final FieldComplex<T>[] b)
        throws MathIllegalArgumentException {
        final int n = 2 * a.length;
        final T[] realA      = MathArrays.buildArray(getPartsField(), n);
        final T[] realB      = MathArrays.buildArray(getPartsField(), n);
        final T[] imaginaryA = MathArrays.buildArray(getPartsField(), n);
        final T[] imaginaryB = MathArrays.buildArray(getPartsField(), n);
        for (int i = 0; i < a.length; ++i)  {
            final FieldComplex<T> ai = a[i];
            final FieldComplex<T> bi = b[i];
            realA[2 * i    ]      = ai.real;
            realA[2 * i + 1]      = ai.imaginary.negate();
            realB[2 * i    ]      = bi.real;
            realB[2 * i + 1]      = bi.imaginary;
            imaginaryA[2 * i    ] = ai.real;
            imaginaryA[2 * i + 1] = ai.imaginary;
            imaginaryB[2 * i    ] = bi.imaginary;
            imaginaryB[2 * i + 1] = bi.real;
        }
        return createComplex(real.linearCombination(realA,  realB),
                             real.linearCombination(imaginaryA, imaginaryB));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> linearCombination(final double[] a, final FieldComplex<T>[] b)
        throws MathIllegalArgumentException {
        final int n = a.length;
        final T[] realB      = MathArrays.buildArray(getPartsField(), n);
        final T[] imaginaryB = MathArrays.buildArray(getPartsField(), n);
        for (int i = 0; i < a.length; ++i)  {
            final FieldComplex<T> bi = b[i];
            realB[i]      = bi.real;
            imaginaryB[i] = bi.imaginary;
        }
        return createComplex(real.linearCombination(a,  realB),
                             real.linearCombination(a, imaginaryB));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> linearCombination(final FieldComplex<T> a1, final FieldComplex<T> b1, final FieldComplex<T> a2, final FieldComplex<T> b2) {
        return createComplex(real.linearCombination(a1.real, b1.real,
                                                    a1.imaginary.negate(), b1.imaginary,
                                                    a2.real, b2.real,
                                                    a2.imaginary.negate(), b2.imaginary),
                             real.linearCombination(a1.real, b1.imaginary,
                                                    a1.imaginary, b1.real,
                                                    a2.real, b2.imaginary,
                                                    a2.imaginary, b2.real));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> linearCombination(final double a1, final FieldComplex<T> b1, final double a2, final FieldComplex<T> b2) {
        return createComplex(real.linearCombination(a1, b1.real,
                                                    a2, b2.real),
                             real.linearCombination(a1, b1.imaginary,
                                                    a2, b2.imaginary));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> linearCombination(final FieldComplex<T> a1, final FieldComplex<T> b1,
                                                final FieldComplex<T> a2, final FieldComplex<T> b2,
                                                final FieldComplex<T> a3, final FieldComplex<T> b3) {
        FieldComplex<T>[] a = MathArrays.buildArray(getField(), 3);
        a[0] = a1;
        a[1] = a2;
        a[2] = a3;
        FieldComplex<T>[] b = MathArrays.buildArray(getField(), 3);
        b[0] = b1;
        b[1] = b2;
        b[2] = b3;
        return linearCombination(a, b);
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> linearCombination(final double a1, final FieldComplex<T> b1,
                                                final double a2, final FieldComplex<T> b2,
                                                final double a3, final FieldComplex<T> b3) {
        FieldComplex<T>[] b = MathArrays.buildArray(getField(), 3);
        b[0] = b1;
        b[1] = b2;
        b[2] = b3;
        return linearCombination(new double[]  { a1, a2, a3 }, b);
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> linearCombination(final FieldComplex<T> a1, final FieldComplex<T> b1,
                                                final FieldComplex<T> a2, final FieldComplex<T> b2,
                                                final FieldComplex<T> a3, final FieldComplex<T> b3,
                                                final FieldComplex<T> a4, final FieldComplex<T> b4) {
        FieldComplex<T>[] a = MathArrays.buildArray(getField(), 4);
        a[0] = a1;
        a[1] = a2;
        a[2] = a3;
        a[3] = a4;
        FieldComplex<T>[] b = MathArrays.buildArray(getField(), 4);
        b[0] = b1;
        b[1] = b2;
        b[2] = b3;
        b[3] = b4;
        return linearCombination(a, b);
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> linearCombination(final double a1, final FieldComplex<T> b1,
                                                final double a2, final FieldComplex<T> b2,
                                                final double a3, final FieldComplex<T> b3,
                                                final double a4, final FieldComplex<T> b4) {
        FieldComplex<T>[] b = MathArrays.buildArray(getField(), 4);
        b[0] = b1;
        b[1] = b2;
        b[2] = b3;
        b[3] = b4;
        return linearCombination(new double[]  { a1, a2, a3, a4 }, b);
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> ceil() {
        return createComplex(FastMath.ceil(getRealPart()), FastMath.ceil(getImaginaryPart()));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> floor() {
        return createComplex(FastMath.floor(getRealPart()), FastMath.floor(getImaginaryPart()));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> rint() {
        return createComplex(FastMath.rint(getRealPart()), FastMath.rint(getImaginaryPart()));
    }

    /** {@inheritDoc}
     * <p>
     * for complex numbers, the integer n corresponding to {@code this.subtract(remainder(a)).divide(a)}
     * is a <a href="https://en.wikipedia.org/wiki/Gaussian_integer">Wikipedia - Gaussian integer</a>.
     * </p>
     */
    @Override
    public FieldComplex<T> remainder(final double a) {
        return createComplex(FastMath.IEEEremainder(getRealPart(), a), FastMath.IEEEremainder(getImaginaryPart(), a));
    }

    /** {@inheritDoc}
     * <p>
     * for complex numbers, the integer n corresponding to {@code this.subtract(remainder(a)).divide(a)}
     * is a <a href="https://en.wikipedia.org/wiki/Gaussian_integer">Wikipedia - Gaussian integer</a>.
     * </p>
     */
    @Override
    public FieldComplex<T> remainder(final FieldComplex<T> a) {
        final FieldComplex<T> complexQuotient = divide(a);
        final T  qRInt           = FastMath.rint(complexQuotient.real);
        final T  qIInt           = FastMath.rint(complexQuotient.imaginary);
        return createComplex(real.subtract(qRInt.multiply(a.real)).add(qIInt.multiply(a.imaginary)),
                             imaginary.subtract(qRInt.multiply(a.imaginary)).subtract(qIInt.multiply(a.real)));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> sign() {
        if (isNaN() || isZero()) {
            return this;
        } else {
            return this.divide(FastMath.hypot(real, imaginary));
        }
    }

    /** {@inheritDoc}
     * <p>
     * The signs of real and imaginary parts are copied independently.
     * </p>
     */
    @Override
    public FieldComplex<T> copySign(final FieldComplex<T> z) {
        return createComplex(FastMath.copySign(getRealPart(), z.getRealPart()),
                             FastMath.copySign(getImaginaryPart(), z.getImaginaryPart()));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> copySign(double r) {
        return createComplex(FastMath.copySign(getRealPart(), r), FastMath.copySign(getImaginaryPart(), r));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> toDegrees() {
        return createComplex(FastMath.toDegrees(getRealPart()), FastMath.toDegrees(getImaginaryPart()));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> toRadians() {
        return createComplex(FastMath.toRadians(getRealPart()), FastMath.toRadians(getImaginaryPart()));
    }

    /** {@inheritDoc} */
    @Override
    public FieldComplex<T> getPi() {
        return getPi(getPartsField());
    }

}
