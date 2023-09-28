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

package org.hipparchus.complex;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.FieldSinhCosh;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;
import org.hipparchus.util.SinCos;
import org.hipparchus.util.SinhCosh;

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
 */
public class Complex implements CalculusFieldElement<Complex>, Comparable<Complex>, Serializable  {
    /** The square root of -1. A number representing "0.0 + 1.0i". */
    public static final Complex I = new Complex(0.0, 1.0);
    /** The square root of -1. A number representing "0.0 - 1.0i".
     * @since 1.7
     */
    public static final Complex MINUS_I = new Complex(0.0, -1.0);
    // CHECKSTYLE: stop ConstantName
    /** A complex number representing "NaN + NaNi". */
    public static final Complex NaN = new Complex(Double.NaN, Double.NaN);
    // CHECKSTYLE: resume ConstantName
    /** A complex number representing "+INF + INFi" */
    public static final Complex INF = new Complex(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
    /** A complex number representing "1.0 + 0.0i". */
    public static final Complex ONE = new Complex(1.0, 0.0);
    /** A complex number representing "-1.0 + 0.0i".
     * @since 1.7
     */
    public static final Complex MINUS_ONE = new Complex(-1.0, 0.0);
    /** A complex number representing "0.0 + 0.0i". */
    public static final Complex ZERO = new Complex(0.0, 0.0);
    /** A complex number representing "π + 0.0i". */
    public static final Complex PI   = new Complex(FastMath.PI, 0.0);

    /** A real number representing log(10). */
    private static final double LOG10 = 2.302585092994045684;

    /** Serializable version identifier */
    private static final long serialVersionUID = 20160305L;

    /** The imaginary part. */
    private final double imaginary;
    /** The real part. */
    private final double real;
    /** Record whether this complex number is equal to NaN. */
    private final transient boolean isNaN;
    /** Record whether this complex number is infinite. */
    private final transient boolean isInfinite;

    /**
     * Create a complex number given only the real part.
     *
     * @param real Real part.
     */
    public Complex(double real) {
        this(real, 0.0);
    }

    /**
     * Create a complex number given the real and imaginary parts.
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     */
    public Complex(double real, double imaginary) {
        this.real = real;
        this.imaginary = imaginary;

        isNaN = Double.isNaN(real) || Double.isNaN(imaginary);
        isInfinite = !isNaN &&
            (Double.isInfinite(real) || Double.isInfinite(imaginary));
    }

    /**
     * Return the absolute value of this complex number.
     * Returns {@code NaN} if either real or imaginary part is {@code NaN}
     * and {@code Double.POSITIVE_INFINITY} if neither part is {@code NaN},
     * but at least one part is infinite.
     *
     * @return the norm.
     * @since 2.0
     */
    @Override
    public Complex abs() {
        // we check NaN here because FastMath.hypot checks it after infinity
        return isNaN ? NaN : createComplex(FastMath.hypot(real, imaginary), 0.0);
    }

    /** {@inheritDoc} */
    @Override
    public double norm() {
        // we check NaN here because FastMath.hypot checks it after infinity
        return isNaN ? Double.NaN : FastMath.hypot(real, imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is
     * {@code (this + addend)}.
     * Uses the definitional formula
     * <p>
     *   {@code (a + bi) + (c + di) = (a+c) + (b+d)i}
     * </p>
     * If either {@code this} or {@code addend} has a {@code NaN} value in
     * either part, {@link #NaN} is returned; otherwise {@code Infinite}
     * and {@code NaN} values are returned in the parts of the result
     * according to the rules for {@link java.lang.Double} arithmetic.
     *
     * @param  addend Value to be added to this {@code Complex}.
     * @return {@code this + addend}.
     * @throws NullArgumentException if {@code addend} is {@code null}.
     */
    @Override
    public Complex add(Complex addend) throws NullArgumentException {
        MathUtils.checkNotNull(addend);
        if (isNaN || addend.isNaN) {
            return NaN;
        }

        return createComplex(real + addend.getRealPart(),
                             imaginary + addend.getImaginaryPart());
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this + addend)},
     * with {@code addend} interpreted as a real number.
     *
     * @param addend Value to be added to this {@code Complex}.
     * @return {@code this + addend}.
     * @see #add(Complex)
     */
    @Override
    public Complex add(double addend) {
        if (isNaN || Double.isNaN(addend)) {
            return NaN;
        }

        return createComplex(real + addend, imaginary);
    }

     /**
     * Returns the conjugate of this complex number.
     * The conjugate of {@code a + bi} is {@code a - bi}.
     * <p>
     * {@link #NaN} is returned if either the real or imaginary
     * part of this Complex number equals {@code Double.NaN}.
     * </p><p>
     * If the imaginary part is infinite, and the real part is not
     * {@code NaN}, the returned value has infinite imaginary part
     * of the opposite sign, e.g. the conjugate of
     * {@code 1 + POSITIVE_INFINITY i} is {@code 1 - NEGATIVE_INFINITY i}.
     * </p>
     * @return the conjugate of this Complex object.
     */
    public Complex conjugate() {
        if (isNaN) {
            return NaN;
        }

        return createComplex(real, -imaginary);
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
     *   in either part, {@link #NaN} is returned.
     *  </li>
     *  <li>If {@code divisor} equals {@link #ZERO}, {@link #NaN} is returned.
     *  </li>
     *  <li>If {@code this} and {@code divisor} are both infinite,
     *   {@link #NaN} is returned.
     *  </li>
     *  <li>If {@code this} is finite (i.e., has no {@code Infinite} or
     *   {@code NaN} parts) and {@code divisor} is infinite (one or both parts
     *   infinite), {@link #ZERO} is returned.
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
    public Complex divide(Complex divisor)
        throws NullArgumentException {
        MathUtils.checkNotNull(divisor);
        if (isNaN || divisor.isNaN) {
            return NaN;
        }

        final double c = divisor.getRealPart();
        final double d = divisor.getImaginaryPart();
        if (c == 0.0 && d == 0.0) {
            return NaN;
        }

        if (divisor.isInfinite() && !isInfinite()) {
            return ZERO;
        }

        if (FastMath.abs(c) < FastMath.abs(d)) {
            double q = c / d;
            double denominator = c * q + d;
            return createComplex((real * q + imaginary) / denominator,
                                 (imaginary * q - real) / denominator);
        } else {
            double q = d / c;
            double denominator = d * q + c;
            return createComplex((imaginary * q + real) / denominator,
                                 (imaginary - real * q) / denominator);
        }
    }

    /**
     * Returns a {@code Complex} whose value is {@code (this / divisor)},
     * with {@code divisor} interpreted as a real number.
     *
     * @param  divisor Value by which this {@code Complex} is to be divided.
     * @return {@code this / divisor}.
     * @see #divide(Complex)
     */
    @Override
    public Complex divide(double divisor) {
        if (isNaN || Double.isNaN(divisor)) {
            return NaN;
        }
        if (divisor == 0d) {
            return NaN;
        }
        if (Double.isInfinite(divisor)) {
            return !isInfinite() ? ZERO : NaN;
        }
        return createComplex(real / divisor,
                             imaginary  / divisor);
    }

    /** {@inheritDoc} */
    @Override
    public Complex reciprocal() {
        if (isNaN) {
            return NaN;
        }

        if (real == 0.0 && imaginary == 0.0) {
            return INF;
        }

        if (isInfinite) {
            return ZERO;
        }

        if (FastMath.abs(real) < FastMath.abs(imaginary)) {
            double q = real / imaginary;
            double scale = 1. / (real * q + imaginary);
            return createComplex(scale * q, -scale);
        } else {
            double q = imaginary / real;
            double scale = 1. / (imaginary * q + real);
            return createComplex(scale, -scale * q);
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
        if (other instanceof Complex){
            Complex c = (Complex) other;
            if (c.isNaN) {
                return isNaN;
            } else {
                return MathUtils.equals(real, c.real) &&
                       MathUtils.equals(imaginary, c.imaginary);
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
     * @return {@code true} if there are fewer than {@code maxUlps} floating
     * point values between the real (resp. imaginary) parts of {@code x}
     * and {@code y}.
     *
     * @see Precision#equals(double,double,int)
     */
    public static boolean equals(Complex x, Complex y, int maxUlps) {
        return Precision.equals(x.real, y.real, maxUlps) &&
               Precision.equals(x.imaginary, y.imaginary, maxUlps);
    }

    /**
     * Returns {@code true} iff the values are equal as defined by
     * {@link #equals(Complex,Complex,int) equals(x, y, 1)}.
     *
     * @param x First value (cannot be {@code null}).
     * @param y Second value (cannot be {@code null}).
     * @return {@code true} if the values are equal.
     */
    public static boolean equals(Complex x, Complex y) {
        return equals(x, y, 1);
    }

    /**
     * Returns {@code true} if, both for the real part and for the imaginary
     * part, there is no double value strictly between the arguments or the
     * difference between them is within the range of allowed error
     * (inclusive).  Returns {@code false} if either of the arguments is NaN.
     *
     * @param x First value (cannot be {@code null}).
     * @param y Second value (cannot be {@code null}).
     * @param eps Amount of allowed absolute error.
     * @return {@code true} if the values are two adjacent floating point
     * numbers or they are within range of each other.
     *
     * @see Precision#equals(double,double,double)
     */
    public static boolean equals(Complex x, Complex y, double eps) {
        return Precision.equals(x.real, y.real, eps) &&
               Precision.equals(x.imaginary, y.imaginary, eps);
    }

    /**
     * Returns {@code true} if, both for the real part and for the imaginary
     * part, there is no double value strictly between the arguments or the
     * relative difference between them is smaller or equal to the given
     * tolerance. Returns {@code false} if either of the arguments is NaN.
     *
     * @param x First value (cannot be {@code null}).
     * @param y Second value (cannot be {@code null}).
     * @param eps Amount of allowed relative error.
     * @return {@code true} if the values are two adjacent floating point
     * numbers or they are within range of each other.
     *
     * @see Precision#equalsWithRelativeTolerance(double,double,double)
     */
    public static boolean equalsWithRelativeTolerance(Complex x,
                                                      Complex y,
                                                      double eps) {
        return Precision.equalsWithRelativeTolerance(x.real, y.real, eps) &&
               Precision.equalsWithRelativeTolerance(x.imaginary, y.imaginary, eps);
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
        return 37 * (17 * MathUtils.hash(imaginary) +
            MathUtils.hash(real));
    }

    /** {@inheritDoc}
     * <p>
     * This implementation considers +0.0 and -0.0 to be equal for both
     * real and imaginary components.
     * </p>
     * @since 1.8
     */
    @Override
    public boolean isZero() {
        return real == 0.0 && imaginary == 0.0;
    }

    /**
     * Access the imaginary part.
     *
     * @return the imaginary part.
     */
    public double getImaginary() {
        return imaginary;
    }

    /**
     * Access the imaginary part.
     *
     * @return the imaginary part.
     * @since 2.0
     */
    public double getImaginaryPart() {
        return imaginary;
    }

    /**
     * Access the real part.
     *
     * @return the real part.
     */
    @Override
    public double getReal() {
        return real;
    }

    /**
     * Access the real part.
     *
     * @return the real part.
     * @since 2.0
     */
    public double getRealPart() {
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
     * @since 1.7
     */
    public boolean isReal() {
        return imaginary == 0.0;
    }

    /** Check whether the instance is an integer (i.e. imaginary part is zero and real part has no fractional part).
     * @return true if imaginary part is zero and real part has no fractional part
     * @since 1.7
     */
    public boolean isMathematicalInteger() {
        return isReal() && Precision.isMathematicalInteger(real);
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
     * Returns {@link #NaN} if either {@code this} or {@code factor} has one or
     * more {@code NaN} parts.
     * <p>
     * Returns {@link #INF} if neither {@code this} nor {@code factor} has one
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
    public Complex multiply(Complex factor)
        throws NullArgumentException {
        MathUtils.checkNotNull(factor);
        if (isNaN || factor.isNaN) {
            return NaN;
        }
        if (Double.isInfinite(real) ||
            Double.isInfinite(imaginary) ||
            Double.isInfinite(factor.real) ||
            Double.isInfinite(factor.imaginary)) {
            // we don't use isInfinite() to avoid testing for NaN again
            return INF;
        }
        return createComplex(MathArrays.linearCombination(real, factor.real, -imaginary, factor.imaginary),
                             MathArrays.linearCombination(real, factor.imaginary, imaginary, factor.real));
    }

    /**
     * Returns a {@code Complex} whose value is {@code this * factor}, with {@code factor}
     * interpreted as a integer number.
     *
     * @param  factor value to be multiplied by this {@code Complex}.
     * @return {@code this * factor}.
     * @see #multiply(Complex)
     */
    @Override
    public Complex multiply(final int factor) {
        if (isNaN) {
            return NaN;
        }
        if (Double.isInfinite(real) ||
            Double.isInfinite(imaginary)) {
            return INF;
        }
        return createComplex(real * factor, imaginary * factor);
    }

    /**
     * Returns a {@code Complex} whose value is {@code this * factor}, with {@code factor}
     * interpreted as a real number.
     *
     * @param  factor value to be multiplied by this {@code Complex}.
     * @return {@code this * factor}.
     * @see #multiply(Complex)
     */
    @Override
    public Complex multiply(double factor) {
        if (isNaN || Double.isNaN(factor)) {
            return NaN;
        }
        if (Double.isInfinite(real) ||
            Double.isInfinite(imaginary) ||
            Double.isInfinite(factor)) {
            // we don't use isInfinite() to avoid testing for NaN again
            return INF;
        }
        return createComplex(real * factor, imaginary * factor);
    }

    /** Compute this * i.
     * @return this * i
     * @since 2.0
     */
    public Complex multiplyPlusI() {
        return createComplex(-imaginary, real);
    }

    /** Compute this *- -i.
     * @return this * i
     * @since 2.0
     */
    public Complex multiplyMinusI() {
        return createComplex(imaginary, -real);
    }

    /**
     * Returns a {@code Complex} whose value is {@code (-this)}.
     * Returns {@code NaN} if either real or imaginary
     * part of this Complex number is {@code Double.NaN}.
     *
     * @return {@code -this}.
     */
    @Override
    public Complex negate() {
        if (isNaN) {
            return NaN;
        }

        return createComplex(-real, -imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is
     * {@code (this - subtrahend)}.
     * Uses the definitional formula
     * <p>
     *  {@code (a + bi) - (c + di) = (a-c) + (b-d)i}
     * </p>
     * If either {@code this} or {@code subtrahend} has a {@code NaN]} value in either part,
     * {@link #NaN} is returned; otherwise infinite and {@code NaN} values are
     * returned in the parts of the result according to the rules for
     * {@link java.lang.Double} arithmetic.
     *
     * @param  subtrahend value to be subtracted from this {@code Complex}.
     * @return {@code this - subtrahend}.
     * @throws NullArgumentException if {@code subtrahend} is {@code null}.
     */
    @Override
    public Complex subtract(Complex subtrahend)
        throws NullArgumentException {
        MathUtils.checkNotNull(subtrahend);
        if (isNaN || subtrahend.isNaN) {
            return NaN;
        }

        return createComplex(real - subtrahend.getRealPart(),
                             imaginary - subtrahend.getImaginaryPart());
    }

    /**
     * Returns a {@code Complex} whose value is
     * {@code (this - subtrahend)}.
     *
     * @param  subtrahend value to be subtracted from this {@code Complex}.
     * @return {@code this - subtrahend}.
     * @see #subtract(Complex)
     */
    @Override
    public Complex subtract(double subtrahend) {
        if (isNaN || Double.isNaN(subtrahend)) {
            return NaN;
        }
        return createComplex(real - subtrahend, imaginary);
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/InverseCosine.html" TARGET="_top">
     * inverse cosine</a> of this complex number.
     * Implements the formula:
     * <p>
     *  {@code acos(z) = -i (log(z + i (sqrt(1 - z<sup>2</sup>))))}
     * </p>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN} or infinite.
     *
     * @return the inverse cosine of this complex number.
     */
    @Override
    public Complex acos() {
        if (isNaN) {
            return NaN;
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
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN} or infinite.</p>
     *
     * @return the inverse sine of this complex number.
     */
    @Override
    public Complex asin() {
        if (isNaN) {
            return NaN;
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
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN} or infinite.</p>
     *
     * @return the inverse tangent of this complex number
     */
    @Override
    public Complex atan() {
        if (isNaN) {
            return NaN;
        }

        if (real == 0.0) {

            // singularity at ±i
            if (imaginary * imaginary - 1.0 == 0.0) {
                return NaN;
            }

            // branch cut on imaginary axis
            final Complex tmp = createComplex((1 + imaginary) / (1 - imaginary), 0.0).log().multiplyPlusI().multiply(0.5);
            return createComplex(FastMath.copySign(tmp.real, real), tmp.imaginary);

        } else {
            // regular formula
            final Complex n = createComplex(1 + imaginary, -real);
            final Complex d = createComplex(1 - imaginary,  real);
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
     * Returns {@link Complex#NaN} if either real or imaginary part of the
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
    public Complex cos() {
        if (isNaN) {
            return NaN;
        }

        final SinCos   scr  = FastMath.sinCos(real);
        final SinhCosh schi = FastMath.sinhCosh(imaginary);
        return createComplex(scr.cos() * schi.cosh(), -scr.sin() * schi.sinh());
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
     * Returns {@link Complex#NaN} if either real or imaginary part of the
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
    public Complex cosh() {
        if (isNaN) {
            return NaN;
        }

        final SinhCosh schr = FastMath.sinhCosh(real);
        final SinCos   sci  = FastMath.sinCos(imaginary);
        return createComplex(schr.cosh() * sci.cos(), schr.sinh() * sci.sin());
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
     * Returns {@link Complex#NaN} if either real or imaginary part of the
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
    public Complex exp() {
        if (isNaN) {
            return NaN;
        }

        final double expReal = FastMath.exp(real);
        final SinCos sc      = FastMath.sinCos(imaginary);
        return createComplex(expReal * sc.cos(), expReal * sc.sin());
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex expm1() {
        if (isNaN) {
            return NaN;
        }

        final double expm1Real = FastMath.expm1(real);
        final SinCos sc        = FastMath.sinCos(imaginary);
        return createComplex(expm1Real * sc.cos(), expm1Real * sc.sin());
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
     * {@code |a + bi|} is the modulus, {@link Complex#abs},  and
     * {@code arg(a + bi) = }{@link FastMath#atan2}(b, a).
     * <p>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
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
    public Complex log() {
        if (isNaN) {
            return NaN;
        }

        return createComplex(FastMath.log(FastMath.hypot(real, imaginary)),
                             FastMath.atan2(imaginary, real));
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex log1p() {
        return add(1.0).log();
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex log10() {
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
    public Complex pow(Complex x)
        throws NullArgumentException {

        MathUtils.checkNotNull(x);

        if (x.imaginary == 0.0) {
            final int nx = (int) FastMath.rint(x.real);
            if (x.real == nx) {
                // integer power
                return pow(nx);
            } else if (this.imaginary == 0.0) {
                // check real implementation that handles a bunch of special cases
                final double realPow = FastMath.pow(this.real, x.real);
                if (Double.isFinite(realPow)) {
                    return createComplex(realPow, 0);
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
    @Override
    public Complex pow(double x) {

        final int nx = (int) FastMath.rint(x);
        if (x == nx) {
            // integer power
            return pow(nx);
        } else if (this.imaginary == 0.0) {
            // check real implementation that handles a bunch of special cases
            final double realPow = FastMath.pow(this.real, x);
            if (Double.isFinite(realPow)) {
                return createComplex(realPow, 0);
            }
        }

        // generic implementation
        return this.log().multiply(x).exp();

    }

     /** {@inheritDoc}
      * @since 1.7
      */
    @Override
    public Complex pow(final int n) {

        Complex result = ONE;
        final boolean invert;
        int p = n;
        if (p < 0) {
            invert = true;
            p = -p;
        } else {
            invert = false;
        }

        // Exponentiate by successive squaring
        Complex square = this;
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
     * Returns {@link Complex#NaN} if either real or imaginary part of the
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
    public Complex sin() {
        if (isNaN) {
            return NaN;
        }

        final SinCos   scr  = FastMath.sinCos(real);
        final SinhCosh schi = FastMath.sinhCosh(imaginary);
        return createComplex(scr.sin() * schi.cosh(), scr.cos() * schi.sinh());

    }

    /** {@inheritDoc}
     */
    @Override
    public FieldSinCos<Complex> sinCos() {
        if (isNaN) {
            return new FieldSinCos<>(NaN, NaN);
        }

        final SinCos scr = FastMath.sinCos(real);
        final SinhCosh schi = FastMath.sinhCosh(imaginary);
        return new FieldSinCos<>(createComplex(scr.sin() * schi.cosh(),  scr.cos() * schi.sinh()),
                                 createComplex(scr.cos() * schi.cosh(), -scr.sin() * schi.sinh()));
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex atan2(Complex x) {

        // compute r = sqrt(x^2+y^2)
        final Complex r = x.multiply(x).add(multiply(this)).sqrt();

        if (FastMath.copySign(1.0, x.real) >= 0) {
            // compute atan2(y, x) = 2 atan(y / (r + x))
            return divide(r.add(x)).atan().multiply(2);
        } else {
            // compute atan2(y, x) = +/- pi - 2 atan(y / (r - x))
            return divide(r.subtract(x)).atan().multiply(-2).add(FastMath.PI);
        }
    }

    /** {@inheritDoc}
     * <p>
     * Branch cuts are on the real axis, below +1.
     * </p>
     * @since 1.7
     */
    @Override
    public Complex acosh() {
        final Complex sqrtPlus  = add(1).sqrt();
        final Complex sqrtMinus = subtract(1).sqrt();
        return add(sqrtPlus.multiply(sqrtMinus)).log();
    }

    /** {@inheritDoc}
     * <p>
     * Branch cuts are on the imaginary axis, above +i and below -i.
     * </p>
     * @since 1.7
     */
    @Override
    public Complex asinh() {
        return add(multiply(this).add(1.0).sqrt()).log();
    }

    /** {@inheritDoc}
     * <p>
     * Branch cuts are on the real axis, above +1 and below -1.
     * </p>
     * @since 1.7
     */
    @Override
    public Complex atanh() {
        final Complex logPlus  = add(1).log();
        final Complex logMinus = createComplex(1 - real, -imaginary).log();
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
     * Returns {@link Complex#NaN} if either real or imaginary part of the
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
    public Complex sinh() {
        if (isNaN) {
            return NaN;
        }

        final SinhCosh schr = FastMath.sinhCosh(real);
        final SinCos   sci  = FastMath.sinCos(imaginary);
        return createComplex(schr.sinh() * sci.cos(), schr.cosh() * sci.sin());
    }

    /** {@inheritDoc}
     */
    @Override
    public FieldSinhCosh<Complex> sinhCosh() {
        if (isNaN) {
            return new FieldSinhCosh<>(NaN, NaN);
        }

        final SinhCosh schr = FastMath.sinhCosh(real);
        final SinCos   sci  = FastMath.sinCos(imaginary);
        return new FieldSinhCosh<>(createComplex(schr.sinh() * sci.cos(), schr.cosh() * sci.sin()),
                                   createComplex(schr.cosh() * sci.cos(), schr.sinh() * sci.sin()));
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
     * <li>{@code |a| = }{@link FastMath#abs(double) abs(a)}</li>
     * <li>{@code |a + bi| = }{@link FastMath#hypot(double, double) hypot(a, b)}</li>
     * <li>{@code sign(b) = }{@link FastMath#copySign(double, double) copySign(1, b)}
     * </ul>
     * The real part is therefore always nonnegative.
     * <p>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
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
    public Complex sqrt() {
        if (isNaN) {
            return NaN;
        }

        if (real == 0.0 && imaginary == 0.0) {
            return ZERO;
        }

        double t = FastMath.sqrt((FastMath.abs(real) + FastMath.hypot(real, imaginary)) * 0.5);
        if (FastMath.copySign(1, real) >= 0.0) {
            return createComplex(t, imaginary / (2.0 * t));
        } else {
            return createComplex(FastMath.abs(imaginary) / (2.0 * t),
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
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * </p>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     *
     * @return the square root of <code>1 - this<sup>2</sup></code>.
     */
    public Complex sqrt1z() {
        final Complex t2 = this.multiply(this);
        return createComplex(1 - t2.real, -t2.imaginary).sqrt();
    }

    /** {@inheritDoc}
     * <p>
     * This implementation compute the principal cube root by using a branch cut along real negative axis.
     * </p>
     * @since 1.7
     */
    @Override
    public Complex cbrt() {
        final double magnitude = FastMath.cbrt(norm());
        final SinCos sc        = FastMath.sinCos(getArgument() / 3);
        return createComplex(magnitude * sc.cos(), magnitude * sc.sin());
    }

    /** {@inheritDoc}
     * <p>
     * This implementation compute the principal n<sup>th</sup> root by using a branch cut along real negative axis.
     * </p>
     * @since 1.7
     */
    @Override
    public Complex rootN(int n) {
        final double magnitude = FastMath.pow(norm(), 1.0 / n);
        final SinCos sc        = FastMath.sinCos(getArgument() / n);
        return createComplex(magnitude * sc.cos(), magnitude * sc.sin());
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
     * Returns {@link Complex#NaN} if either real or imaginary part of the
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
    public Complex tan() {
        if (isNaN || Double.isInfinite(real)) {
            return NaN;
        }
        if (imaginary > 20.0) {
            return I;
        }
        if (imaginary < -20.0) {
            return MINUS_I;
        }

        final SinCos sc2r = FastMath.sinCos(2.0 * real);
        double imaginary2 = 2.0 * imaginary;
        double d = sc2r.cos() + FastMath.cosh(imaginary2);

        return createComplex(sc2r.sin() / d, FastMath.sinh(imaginary2) / d);

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
     * Returns {@link Complex#NaN} if either real or imaginary part of the
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
    public Complex tanh() {
        if (isNaN || Double.isInfinite(imaginary)) {
            return NaN;
        }
        if (real > 20.0) {
            return ONE;
        }
        if (real < -20.0) {
            return MINUS_ONE;
        }
        double real2 = 2.0 * real;
        final SinCos sc2i = FastMath.sinCos(2.0 * imaginary);
        double d = FastMath.cosh(real2) + sc2i.cos();

        return createComplex(FastMath.sinh(real2) / d, sc2i.sin() / d);
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
    public double getArgument() {
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
     * one element, {@link #NaN} is returned.
     * if neither part is NaN, but at least one part is infinite, the result
     * is a one-element list containing {@link #INF}.
     *
     * @param n Degree of root.
     * @return a List of all {@code n}-th roots of {@code this}.
     * @throws MathIllegalArgumentException if {@code n <= 0}.
     */
    public List<Complex> nthRoot(int n) throws MathIllegalArgumentException {

        if (n <= 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.CANNOT_COMPUTE_NTH_ROOT_FOR_NEGATIVE_N,
                                                   n);
        }

        final List<Complex> result = new ArrayList<>();

        if (isNaN) {
            result.add(NaN);
            return result;
        }
        if (isInfinite()) {
            result.add(INF);
            return result;
        }

        // nth root of abs -- faster / more accurate to use a solver here?
        final double nthRootOfAbs = FastMath.pow(FastMath.hypot(real, imaginary), 1.0 / n);

        // Compute nth roots of complex number with k = 0, 1, ... n-1
        final double nthPhi = getArgument() / n;
        final double slice = 2 * FastMath.PI / n;
        double innerPart = nthPhi;
        for (int k = 0; k < n ; k++) {
            // inner part
            final SinCos scInner = FastMath.sinCos(innerPart);
            final double realPart = nthRootOfAbs *  scInner.cos();
            final double imaginaryPart = nthRootOfAbs *  scInner.sin();
            result.add(createComplex(realPart, imaginaryPart));
            innerPart += slice;
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
     * @see #valueOf(double, double)
     */
    protected Complex createComplex(double realPart,
                                    double imaginaryPart) {
        return new Complex(realPart, imaginaryPart);
    }

    /**
     * Create a complex number given the real and imaginary parts.
     *
     * @param realPart Real part.
     * @param imaginaryPart Imaginary part.
     * @return a Complex instance.
     */
    public static Complex valueOf(double realPart,
                                  double imaginaryPart) {
        if (Double.isNaN(realPart) ||
            Double.isNaN(imaginaryPart)) {
            return NaN;
        }
        return new Complex(realPart, imaginaryPart);
    }

    /**
     * Create a complex number given only the real part.
     *
     * @param realPart Real part.
     * @return a Complex instance.
     */
    public static Complex valueOf(double realPart) {
        if (Double.isNaN(realPart)) {
            return NaN;
        }
        return new Complex(realPart);
    }

    /** {@inheritDoc} */
    @Override
    public Complex newInstance(double realPart) {
        return valueOf(realPart);
    }

    /**
     * Resolve the transient fields in a deserialized Complex Object.
     * Subclasses will need to override {@link #createComplex} to
     * deserialize properly.
     *
     * @return A Complex instance with all fields resolved.
     */
    protected final Object readResolve() {
        return createComplex(real, imaginary);
    }

    /** {@inheritDoc} */
    @Override
    public ComplexField getField() {
        return ComplexField.getInstance();
    }

    /** {@inheritDoc} */
    @Override
    public String toString() {
        return "(" + real + ", " + imaginary + ")";
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex scalb(int n) {
        return createComplex(FastMath.scalb(real, n), FastMath.scalb(imaginary, n));
    }

    /** {@inheritDoc}
     */
    @Override
    public Complex ulp() {
        return createComplex(FastMath.ulp(real), FastMath.ulp(imaginary));
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex hypot(Complex y) {
        if (isInfinite() || y.isInfinite()) {
            return INF;
        } else if (isNaN() || y.isNaN()) {
            return NaN;
        } else {
            return multiply(this).add(y.multiply(y)).sqrt();
        }
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex linearCombination(final Complex[] a, final Complex[] b)
        throws MathIllegalArgumentException {
        final int n = 2 * a.length;
        final double[] realA      = new double[n];
        final double[] realB      = new double[n];
        final double[] imaginaryA = new double[n];
        final double[] imaginaryB = new double[n];
        for (int i = 0; i < a.length; ++i)  {
            final Complex ai = a[i];
            final Complex bi = b[i];
            realA[2 * i    ]      = +ai.real;
            realA[2 * i + 1]      = -ai.imaginary;
            realB[2 * i    ]      = +bi.real;
            realB[2 * i + 1]      = +bi.imaginary;
            imaginaryA[2 * i    ] = +ai.real;
            imaginaryA[2 * i + 1] = +ai.imaginary;
            imaginaryB[2 * i    ] = +bi.imaginary;
            imaginaryB[2 * i + 1] = +bi.real;
        }
        return createComplex(MathArrays.linearCombination(realA,  realB),
                             MathArrays.linearCombination(imaginaryA, imaginaryB));
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex linearCombination(final double[] a, final Complex[] b)
        throws MathIllegalArgumentException {
        final int n = a.length;
        final double[] realB      = new double[n];
        final double[] imaginaryB = new double[n];
        for (int i = 0; i < a.length; ++i)  {
            final Complex bi = b[i];
            realB[i]      = +bi.real;
            imaginaryB[i] = +bi.imaginary;
        }
        return createComplex(MathArrays.linearCombination(a,  realB),
                             MathArrays.linearCombination(a, imaginaryB));
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex linearCombination(final Complex a1, final Complex b1, final Complex a2, final Complex b2) {
        return createComplex(MathArrays.linearCombination(+a1.real, b1.real,
                                                          -a1.imaginary, b1.imaginary,
                                                          +a2.real, b2.real,
                                                          -a2.imaginary, b2.imaginary),
                             MathArrays.linearCombination(+a1.real, b1.imaginary,
                                                          +a1.imaginary, b1.real,
                                                          +a2.real, b2.imaginary,
                                                          +a2.imaginary, b2.real));
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex linearCombination(final double a1, final Complex b1, final double a2, final Complex b2) {
        return createComplex(MathArrays.linearCombination(a1, b1.real,
                                                          a2, b2.real),
                             MathArrays.linearCombination(a1, b1.imaginary,
                                                          a2, b2.imaginary));
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex linearCombination(final Complex a1, final Complex b1,
                                     final Complex a2, final Complex b2,
                                     final Complex a3, final Complex b3) {
        return linearCombination(new Complex[] { a1, a2, a3 },
                                 new Complex[] { b1, b2, b3 });
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex linearCombination(final double a1, final Complex b1,
                                     final double a2, final Complex b2,
                                     final double a3, final Complex b3) {
        return linearCombination(new double[]  { a1, a2, a3 },
                                 new Complex[] { b1, b2, b3 });
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex linearCombination(final Complex a1, final Complex b1,
                                     final Complex a2, final Complex b2,
                                     final Complex a3, final Complex b3,
                                     final Complex a4, final Complex b4) {
        return linearCombination(new Complex[] { a1, a2, a3, a4 },
                                 new Complex[] { b1, b2, b3, b4 });
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex linearCombination(final double a1, final Complex b1,
                                     final double a2, final Complex b2,
                                     final double a3, final Complex b3,
                                     final double a4, final Complex b4) {
        return linearCombination(new double[]  { a1, a2, a3, a4 },
                                 new Complex[] { b1, b2, b3, b4 });
    }

    /** {@inheritDoc} */
    @Override
    public Complex getPi() {
        return PI;
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex ceil() {
        return createComplex(FastMath.ceil(getRealPart()), FastMath.ceil(getImaginaryPart()));
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex floor() {
        return createComplex(FastMath.floor(getRealPart()), FastMath.floor(getImaginaryPart()));
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex rint() {
        return createComplex(FastMath.rint(getRealPart()), FastMath.rint(getImaginaryPart()));
    }

    /** {@inheritDoc}
     * <p>
     * for complex numbers, the integer n corresponding to {@code this.subtract(remainder(a)).divide(a)}
     * is a <a href="https://en.wikipedia.org/wiki/Gaussian_integer">Wikipedia - Gaussian integer</a>.
     * </p>
     * @since 1.7
     */
    @Override
    public Complex remainder(final double a) {
        return createComplex(FastMath.IEEEremainder(getRealPart(), a), FastMath.IEEEremainder(getImaginaryPart(), a));
    }

    /** {@inheritDoc}
     * <p>
     * for complex numbers, the integer n corresponding to {@code this.subtract(remainder(a)).divide(a)}
     * is a <a href="https://en.wikipedia.org/wiki/Gaussian_integer">Wikipedia - Gaussian integer</a>.
     * </p>
     * @since 1.7
     */
    @Override
    public Complex remainder(final Complex a) {
        final Complex complexQuotient = divide(a);
        final double  qRInt           = FastMath.rint(complexQuotient.real);
        final double  qIInt           = FastMath.rint(complexQuotient.imaginary);
        return createComplex(real - qRInt * a.real + qIInt * a.imaginary,
                             imaginary - qRInt * a.imaginary - qIInt * a.real);
    }

    /** {@inheritDoc}
     * @since 2.0
     */
    @Override
    public Complex sign() {
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
     * @since 1.7
     */
    @Override
    public Complex copySign(final Complex z) {
        return createComplex(FastMath.copySign(getRealPart(), z.getRealPart()),
                             FastMath.copySign(getImaginaryPart(), z.getImaginaryPart()));
    }

    /** {@inheritDoc}
     * @since 1.7
     */
    @Override
    public Complex copySign(double r) {
        return createComplex(FastMath.copySign(getRealPart(), r), FastMath.copySign(getImaginaryPart(), r));
    }

    /** {@inheritDoc} */
    @Override
    public Complex toDegrees() {
        return createComplex(FastMath.toDegrees(getRealPart()), FastMath.toDegrees(getImaginaryPart()));
    }

    /** {@inheritDoc} */
    @Override
    public Complex toRadians() {
        return createComplex(FastMath.toRadians(getRealPart()), FastMath.toRadians(getImaginaryPart()));
    }

    /** {@inheritDoc}
     * <p>
     * Comparison us performed using real ordering as the primary sort order and
     * imaginary ordering as the secondary sort order.
     * </p>
     * @since 3.0
     */
    @Override
    public int compareTo(final Complex o) {
        final int cR = Double.compare(getReal(), o.getReal());
        if (cR == 0) {
            return Double.compare(getImaginary(),o.getImaginary());
        } else {
            return cR;
        }
    }

}
