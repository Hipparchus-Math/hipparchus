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

/**
 * Interface representing a <a href="http://mathworld.wolfram.com/Field.html">field</a>
 * with calculus capabilities (sin, cos, ...).
 * @param <T> the type of the field elements
 * @see FieldElement
 * @since 1.7
 */
public interface CalculusFieldElement<T extends FieldElement<T>> extends FieldElement<T> {

    /** Get the Archimedes constant π.
     * <p>
     * Archimedes constant is the ratio of a circle's circumference to its diameter.
     * </p>
     * @return Archimedes constant π
     * @since 2.0
     */
    T getPi();

    /** Create an instance corresponding to a constant real value.
     * @param value constant real value
     * @return instance corresponding to a constant real value
     */
    T newInstance(double value);

    /** '+' operator.
     * @param a right hand side parameter of the operator
     * @return this+a
     */
    T add(double a);

    /** '-' operator.
     * @param a right hand side parameter of the operator
     * @return this-a
     */
    T subtract(double a);

    /** '&times;' operator.
     * @param a right hand side parameter of the operator
     * @return this&times;a
     */
    T multiply(double a);

    /** '&divide;' operator.
     * @param a right hand side parameter of the operator
     * @return this&divide;a
     */
    T divide(double a);

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
    T ulp();

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

    /** {@inheritDoc} */
    @Override
    T reciprocal();

    /** Square root.
     * @return square root of the instance
     */
    T sqrt();

    /** Cubic root.
     * @return cubic root of the instance
     */
    T cbrt();

    /** N<sup>th</sup> root.
     * @param n order of the root
     * @return n<sup>th</sup> root of the instance
     */
    T rootN(int n);

    /** Power operation.
     * @param p power to apply
     * @return this<sup>p</sup>
     */
    T pow(double p);

    /** Integer power operation.
     * @param n power to apply
     * @return this<sup>n</sup>
     */
    T pow(int n);

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
    FieldSinCos<T> sinCos();

    /** Tangent operation.
     * @return tan(this)
     */
    T tan();

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

    /** Combined hyperbolic sine and sosine operation.
     * @return [sinh(this), cosh(this)]
     * @since 2.0
     */
    FieldSinhCosh<T> sinhCosh();

    /** Hyperbolic tangent operation.
     * @return tanh(this)
     */
    T tanh();

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
    T toDegrees();

    /** Convert degrees to radians, with error of less than 0.5 ULP
     *  @return instance converted into radians
     */
    T toRadians();

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
    T linearCombination(double[] a, T[] b)
        throws MathIllegalArgumentException;

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
    T linearCombination(double a1, T b1, double a2, T b2);

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
    T linearCombination(double a1, T b1, double a2, T b2, double a3, T b3);

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
    T linearCombination(double a1, T b1, double a2, T b2, double a3, T b3, double a4, T b4);

    /** Get the smallest whole number larger than instance.
     * @return ceil(this)
     */
    T ceil();

    /** Get the largest whole number smaller than instance.
     * @return floor(this)
     */
    T floor();

    /** Get the whole number that is the nearest to the instance, or the even one if x is exactly half way between two integers.
     * @return a double number r such that r is an integer r - 0.5 &le; this &le; r + 0.5
     */
    T rint();

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
    T sign();

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
    T copySign(double sign);

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
     * <p>
     * Just another name for {@link #norm()}
     * </p>
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
