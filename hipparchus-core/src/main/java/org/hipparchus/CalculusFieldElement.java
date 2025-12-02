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
package org.hipparchus;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.FieldSinhCosh;
import org.hipparchus.util.MathArrays;

/**
 * Interface representing a <a href="http://mathworld.wolfram.com/Field.html">field</a>
 * with calculus capabilities (sin, cos, ...).
 * @param <T> the type of the field elements
 * @see FieldElement
 * @since 1.7
 */
public interface CalculusFieldElement<T extends FieldElement<T>> extends FieldElement<T> {

    /** Get the addendum to the real value of the number.
     * <p>
     * The addendum is considered to be the part that when added back to
     * the {@link #getReal() real part} recovers the instance. This means
     * that when {@code e.getReal()} is finite (i.e. neither infinite
     * nor NaN), then {@code e.getAddendum().add(e.getReal())} is {@code e}
     * and {@code e.subtract(e.getReal())} is {@code e.getAddendum()}.
     * Beware that for non-finite numbers, these two equalities may not hold.
     * The first equality (with the addition), always holds even for infinity
     * and NaNs if the real part is independent of the addendum (this is the
     * case for all derivatives types, as well as for complex and Dfp, but it
     * is not the case for Tuple and FieldTuple). The second equality (with
     * the subtraction), generally doesn't hold for non-finite numbers, because
     * the subtraction generates NaNs.
     * </p>
     * @return real value
     * @since 4.0
     */
    T getAddendum();

    /** Get the Archimedes constant π.
     * <p>
     * Archimedes constant is the ratio of a circle's circumference to its diameter.
     * </p>
     * @return Archimedes constant π
     * @since 2.0
     */
    default T getPi() {
        return newInstance(FastMath.PI);
    }

    /** Create an instance corresponding to a constant real value.
     * @param value constant real value
     * @return instance corresponding to a constant real value
     */
    T newInstance(double value);

    /** '+' operator.
     * @param a right hand side parameter of the operator
     * @return this+a
     */
    default T add(double a) {
        return add(newInstance(a));
    }

    /** '-' operator.
     * @param a right hand side parameter of the operator
     * @return this-a
     */
    default T subtract(double a) {
        return subtract(newInstance(a));
    }

    /** {@inheritDoc} */
    @Override
    default T subtract(T a) {
        return add(a.negate());
    }

    /** '&times;' operator.
     * @param a right hand side parameter of the operator
     * @return this&times;a
     */
    default T multiply(double a) {
        return multiply(newInstance(a));
    }

    /** {@inheritDoc} */
    @Override
    default T multiply(int n) {
        return multiply((double) n);
    }

    /** '&divide;' operator.
     * @param a right hand side parameter of the operator
     * @return this&divide;a
     */
    default T divide(double a) {
        return divide(newInstance(a));
    }

    /**
     * Return the exponent of the instance, removing the bias.
     * <p>
     * For double numbers of the form 2<sup>x</sup>, the unbiased
     * exponent is exactly x.
     * </p>
     * @return exponent for the instance, without bias
     */
    default int getExponent() {
        return FastMath.getExponent(getReal());
    }

    /**
     * Multiply the instance by a power of 2.
     * @param n power of 2
     * @return this &times; 2<sup>n</sup>
     */
    T scalb(int n);

    /**
     * Compute least significant bit (Unit in Last Position) for a number.
     * @return ulp(this)
     * @since 2.0
     */
    default T ulp() {
        return newInstance(FastMath.ulp(getReal()));
    }

    /**
     * Returns the hypotenuse of a triangle with sides {@code this} and {@code y}
     * - sqrt(<i>this</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
     * avoiding intermediate overflow or underflow.
     *
     * <ul>
     * <li> If either argument is infinite, then the result is positive infinity.</li>
     * <li> else, if either argument is NaN then the result is NaN.</li>
     * </ul>
     *
     * @param y a value
     * @return sqrt(<i>this</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
     * @exception MathIllegalArgumentException if number of free parameters or orders are inconsistent
     */
    T hypot(T y)
        throws MathIllegalArgumentException;

    /** Compute this divided by two.
     * @return a new element representing it halved
     * @since 4.1
     */
    default T half() {
        return divide(2);
    }

    /** Compute this &times; two.
     * @return a new element representing it doubled
     * @since 4.1
     */
    default T twice() {
        return multiply(2);
    }

    /** {@inheritDoc} */
    @Override
    default T divide(T a) {
        return multiply(a.reciprocal());
    }

    /** Square root.
     * @return square root of the instance
     */
    default T sqrt() {
        return rootN(2);
    }

    /** Cubic root.
     * @return cubic root of the instance
     */
    default T cbrt() {
        return rootN(3);
    }

    /** N<sup>th</sup> root.
     * @param n order of the root
     * @return n<sup>th</sup> root of the instance
     */
    default T rootN(int n) {
        return pow(1. / n);
    }

    /** Compute this &times; this.
     * @return a new element representing this &times; this
     * @since 3.1
     */
    default T square() {
        return pow(2);
    }

    /** Power operation.
     * @param p power to apply
     * @return this<sup>p</sup>
     */
    default T pow(double p) {
        return pow(newInstance(p));
    }

    /** Integer power operation.
     * @param n power to apply
     * @return this<sup>n</sup>
     */
    default T pow(int n) {
        return pow((double) n);
    }

    /** Power operation.
     * @param e exponent
     * @return this<sup>e</sup>
     * @exception MathIllegalArgumentException if number of free parameters or orders are inconsistent
     */
    T pow(T e)
        throws MathIllegalArgumentException;

    /** Exponential.
     * @return exponential of the instance
     */
    T exp();

    /** Exponential minus 1.
     * @return exponential minus one of the instance
     */
    T expm1();

    /** Natural logarithm.
     * @return logarithm of the instance
     */
    T log();

    /** Shifted natural logarithm.
     * @return logarithm of one plus the instance
     */
    T log1p();

    /** Base 10 logarithm.
     * @return base 10 logarithm of the instance
     */
    T log10();

    /** Cosine operation.
     * @return cos(this)
     */
    T cos();

    /** Sine operation.
     * @return sin(this)
     */
    T sin();

    /** Combined Sine and Cosine operation.
     * @return [sin(this), cos(this)]
     * @since 1.4
     */
    default FieldSinCos<T> sinCos() {
        return new FieldSinCos<>(sin(), cos());
    }

    /** Tangent operation.
     * @return tan(this)
     */
    default T tan() {
        return sin().divide(cos());
    }

    /** Arc cosine operation.
     * @return acos(this)
     */
    T acos();

    /** Arc sine operation.
     * @return asin(this)
     */
    T asin();

    /** Arc tangent operation.
     * @return atan(this)
     */
    T atan();

    /** Two arguments arc tangent operation.
     * <p>
     * Beware of the order or arguments! As this is based on a
     * two-arguments functions, in order to be consistent with
     * arguments order, the instance is the <em>first</em> argument
     * and the single provided argument is the <em>second</em> argument.
     * In order to be consistent with programming languages {@code atan2},
     * this method computes {@code atan2(this, x)}, i.e. the instance
     * represents the {@code y} argument and the {@code x} argument is
     * the one passed as a single argument. This may seem confusing especially
     * for users of Wolfram alpha, as this site is <em>not</em> consistent
     * with programming languages {@code atan2} two-arguments arc tangent
     * and puts {@code x} as its first argument.
     * </p>
     * @param x second argument of the arc tangent
     * @return atan2(this, x)
     * @exception MathIllegalArgumentException if number of free parameters or orders are inconsistent
     */
    T atan2(T x)
        throws MathIllegalArgumentException;

    /** Hyperbolic cosine operation.
     * @return cosh(this)
     */
    T cosh();

    /** Hyperbolic sine operation.
     * @return sinh(this)
     */
    T sinh();

    /** Combined hyperbolic sine and cosine operation.
     * @return [sinh(this), cosh(this)]
     * @since 2.0
     */
    default FieldSinhCosh<T> sinhCosh() {
        return new FieldSinhCosh<>(sinh(), cosh());
    }

    /** Hyperbolic tangent operation.
     * @return tanh(this)
     */
    default T tanh() {
        return sinh().divide(cosh());
    }

    /** Inverse hyperbolic cosine operation.
     * @return acosh(this)
     */
    T acosh();

    /** Inverse hyperbolic sine operation.
     * @return asin(this)
     */
    T asinh();

    /** Inverse hyperbolic  tangent operation.
     * @return atanh(this)
     */
    T atanh();

    /** Convert radians to degrees, with error of less than 0.5 ULP
     *  @return instance converted into degrees
     */
    default T toDegrees() {
        return multiply(FastMath.toDegrees(1.));
    }

    /** Convert degrees to radians, with error of less than 0.5 ULP
     *  @return instance converted into radians
     */
    default T toRadians() {
        return multiply(FastMath.toRadians(1.));
    }

    /**
     * Compute a linear combination.
     * @param a Factors.
     * @param b Factors.
     * @return <code>&Sigma;<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
     * @throws MathIllegalArgumentException if arrays dimensions don't match
     */
    T linearCombination(T[] a, T[] b)
        throws MathIllegalArgumentException;

    /**
     * Compute a linear combination.
     * @param a Factors.
     * @param b Factors.
     * @return <code>&Sigma;<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
     * @throws MathIllegalArgumentException if arrays dimensions don't match
     */
    default T linearCombination(double[] a, T[] b) throws MathIllegalArgumentException {
        final T[] newInstances = MathArrays.buildArray(getField(), a.length);
        for (int i = 0; i < a.length; i++) {
            newInstances[i] = newInstance(a[i]);
        }
        return linearCombination(newInstances, b);
    }

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @return a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub>
     * @see #linearCombination(FieldElement, FieldElement, FieldElement, FieldElement, FieldElement, FieldElement)
     * @see #linearCombination(FieldElement, FieldElement, FieldElement, FieldElement, FieldElement, FieldElement, FieldElement, FieldElement)
     */
    T linearCombination(T a1, T b1, T a2, T b2);

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @return a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub>
     * @see #linearCombination(double, FieldElement, double, FieldElement, double, FieldElement)
     * @see #linearCombination(double, FieldElement, double, FieldElement, double, FieldElement, double, FieldElement)
     */
    default T linearCombination(double a1, T b1, double a2, T b2) {
        return linearCombination(newInstance(a1), b1, newInstance(a2), b2);
    }

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @param a3 first factor of the third term
     * @param b3 second factor of the third term
     * @return a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub>
     * @see #linearCombination(FieldElement, FieldElement, FieldElement, FieldElement)
     * @see #linearCombination(FieldElement, FieldElement, FieldElement, FieldElement, FieldElement, FieldElement, FieldElement, FieldElement)
     */
    T linearCombination(T a1, T b1, T a2, T b2, T a3, T b3);

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @param a3 first factor of the third term
     * @param b3 second factor of the third term
     * @return a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub>
     * @see #linearCombination(double, FieldElement, double, FieldElement)
     * @see #linearCombination(double, FieldElement, double, FieldElement, double, FieldElement, double, FieldElement)
     */
    default T linearCombination(double a1, T b1, double a2, T b2, double a3, T b3) {
        return linearCombination(newInstance(a1), b1, newInstance(a2), b2, newInstance(a3), b3);
    }

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @param a3 first factor of the third term
     * @param b3 second factor of the third term
     * @param a4 first factor of the fourth term
     * @param b4 second factor of the fourth term
     * @return a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub> +
     * a<sub>4</sub>&times;b<sub>4</sub>
     * @see #linearCombination(FieldElement, FieldElement, FieldElement, FieldElement)
     * @see #linearCombination(FieldElement, FieldElement, FieldElement, FieldElement, FieldElement, FieldElement)
     */
    T linearCombination(T a1, T b1, T a2, T b2, T a3, T b3, T a4, T b4);

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @param a3 first factor of the third term
     * @param b3 second factor of the third term
     * @param a4 first factor of the fourth term
     * @param b4 second factor of the fourth term
     * @return a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub> +
     * a<sub>4</sub>&times;b<sub>4</sub>
     * @see #linearCombination(double, FieldElement, double, FieldElement)
     * @see #linearCombination(double, FieldElement, double, FieldElement, double, FieldElement)
     */
    default T linearCombination(double a1, T b1, double a2, T b2, double a3, T b3, double a4, T b4) {
        return linearCombination(newInstance(a1), b1, newInstance(a2), b2, newInstance(a3), b3,
                newInstance(a4), b4);
    }

    /** Get the smallest whole number larger than instance.
     * @return ceil(this)
     */
    default T ceil() {
        return newInstance(FastMath.ceil(getReal()));
    }

    /** Get the largest whole number smaller than instance.
     * @return floor(this)
     */
    default T floor() {
        return newInstance(FastMath.floor(getReal()));
    }

    /** Get the whole number that is the nearest to the instance, or the even one if x is exactly half way between two integers.
     * @return a double number r such that r is an integer r - 0.5 &le; this &le; r + 0.5
     */
    default T rint() {
        return newInstance(FastMath.rint(getReal()));
    }

    /** IEEE remainder operator.
     * @param a right hand side parameter of the operator
     * @return this - n &times; a where n is the closest integer to this/a
     */
    T remainder(double a);

    /** IEEE remainder operator.
     * @param a right hand side parameter of the operator
     * @return this - n &times; a where n is the closest integer to this/a
     */
    T remainder(T a);

    /** Compute the sign of the instance.
     * The sign is -1 for negative numbers, +1 for positive numbers and 0 otherwise,
     * for Complex number, it is extended on the unit circle (equivalent to z/|z|,
     * with special handling for 0 and NaN)
     * @return -1.0, -0.0, +0.0, +1.0 or NaN depending on sign of a
     */
    default T sign() {
        return newInstance(FastMath.signum(getReal()));
    }

    /**
     * Returns the instance with the sign of the argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param sign the sign for the returned value
     * @return the instance with the same sign as the {@code sign} argument
     */
    T copySign(T sign);

    /**
     * Returns the instance with the sign of the argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param sign the sign for the returned value
     * @return the instance with the same sign as the {@code sign} argument
     */
    default T copySign(double sign) {
        return copySign(newInstance(sign));
    }

    /**
     * Check if the instance is infinite.
     * @return true if the instance is infinite
     */
    default boolean isInfinite() {
        return Double.isInfinite(getReal());
    }

    /**
     * Check if the instance is finite (neither infinite nor NaN).
     * @return true if the instance is finite (neither infinite nor NaN)
     * @since 2.0
     */
    default boolean isFinite() {
        return Double.isFinite(getReal());
    }

    /**
     * Check if the instance is Not a Number.
     * @return true if the instance is Not a Number
     */
    default boolean isNaN() {
        return Double.isNaN(getReal());
    }

    /** norm.
     * @return norm(this)
     * @since 2.0
     */
    default double norm() {
        return abs().getReal();
    }

    /** absolute value.
     * @return abs(this)
     */
    T abs();

    /** Get the closest long to instance real value.
     * @return closest long to {@link #getReal()}
     */
    default long round() {
        return FastMath.round(getReal());
    }

}
