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

package org.hipparchus.util;

import java.util.Arrays;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.FieldElement;
import org.hipparchus.exception.Localizable;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.NullArgumentException;

/**
 * Miscellaneous utility functions.
 *
 * @see ArithmeticUtils
 * @see Precision
 * @see MathArrays
 */
public final class MathUtils {
    /** \(2\pi\) */
    public static final double TWO_PI = 2 * FastMath.PI;

    /** \(\pi^2\) */
    public static final double PI_SQUARED = FastMath.PI * FastMath.PI;

    /** \(\pi/2\). */
    public static final double SEMI_PI = 0.5 * FastMath.PI;

    /**
     * Class contains only static methods.
     */
    private MathUtils() {}


    /**
     * Returns an integer hash code representing the given double value.
     *
     * @param value the value to be hashed
     * @return the hash code
     */
    public static int hash(double value) {
        return Double.hashCode(value);
    }

    /**
     * Returns {@code true} if the values are equal according to semantics of
     * {@link Double#equals(Object)}.
     *
     * @param x Value
     * @param y Value
     * @return {@code Double.valueOf(x).equals(Double.valueOf(y))}
     */
    public static boolean equals(double x, double y) {
        return Double.valueOf(x).equals(y);
    }

    /**
     * Returns an integer hash code representing the given double array.
     *
     * @param value the value to be hashed (may be null)
     * @return the hash code
     */
    public static int hash(double[] value) {
        return Arrays.hashCode(value);
    }

    /**
     * Normalize an angle in a 2&pi; wide interval around a center value.
     * <p>This method has three main uses:</p>
     * <ul>
     *   <li>normalize an angle between 0 and 2&pi;:<br>
     *       {@code a = MathUtils.normalizeAngle(a, FastMath.PI);}</li>
     *   <li>normalize an angle between -&pi; and +&pi;<br>
     *       {@code a = MathUtils.normalizeAngle(a, 0.0);}</li>
     *   <li>compute the angle between two defining angular positions:<br>
     *       {@code angle = MathUtils.normalizeAngle(end, start) - start;}</li>
     * </ul>
     * <p>Note that due to numerical accuracy and since &pi; cannot be represented
     * exactly, the result interval is <em>closed</em>, it cannot be half-closed
     * as would be more satisfactory in a purely mathematical view.</p>
     * @param a angle to normalize
     * @param center center of the desired 2&pi; interval for the result
     * @return a-2k&pi; with integer k and center-&pi; &lt;= a-2k&pi; &lt;= center+&pi;
     */
     public static double normalizeAngle(double a, double center) {
         return a - TWO_PI * FastMath.floor((a + FastMath.PI - center) / TWO_PI);
     }

     /**
      * Normalize an angle in a 2&pi; wide interval around a center value.
      * <p>This method has three main uses:</p>
      * <ul>
      *   <li>normalize an angle between 0 and 2&pi;:<br>
      *       {@code a = MathUtils.normalizeAngle(a, FastMath.PI);}</li>
      *   <li>normalize an angle between -&pi; and +&pi;<br>
      *       {@code a = MathUtils.normalizeAngle(a, zero);}</li>
      *   <li>compute the angle between two defining angular positions:<br>
      *       {@code angle = MathUtils.normalizeAngle(end, start).subtract(start);}</li>
      * </ul>
      * <p>Note that due to numerical accuracy and since &pi; cannot be represented
      * exactly, the result interval is <em>closed</em>, it cannot be half-closed
      * as would be more satisfactory in a purely mathematical view.</p>
      * @param <T> the type of the field elements
      * @param a angle to normalize
      * @param center center of the desired 2&pi; interval for the result
      * @return a-2k&pi; with integer k and center-&pi; &lt;= a-2k&pi; &lt;= center+&pi;
      */
      public static <T extends CalculusFieldElement<T>> T normalizeAngle(T a, T center) {
          return a.subtract(FastMath.floor(a.add(FastMath.PI).subtract(center).divide(TWO_PI)).multiply(TWO_PI));
      }

     /** Find the maximum of two field elements.
      * @param <T> the type of the field elements
      * @param e1 first element
      * @param e2 second element
      * @return max(a1, e2)
      */
     public static <T extends CalculusFieldElement<T>> T max(final T e1, final T e2) {
         return e1.subtract(e2).getReal() >= 0 ? e1 : e2;
     }

     /** Find the minimum of two field elements.
      * @param <T> the type of the field elements
      * @param e1 first element
      * @param e2 second element
      * @return min(a1, e2)
      */
     public static <T extends CalculusFieldElement<T>> T min(final T e1, final T e2) {
         return e1.subtract(e2).getReal() >= 0 ? e2 : e1;
     }

    /**
     * <p>Reduce {@code |a - offset|} to the primary interval
     * {@code [0, |period|)}.</p>
     *
     * <p>Specifically, the value returned is <br>
     * {@code a - |period| * floor((a - offset) / |period|) - offset}.</p>
     *
     * <p>If any of the parameters are {@code NaN} or infinite, the result is
     * {@code NaN}.</p>
     *
     * @param a Value to reduce.
     * @param period Period.
     * @param offset Value that will be mapped to {@code 0}.
     * @return the value, within the interval {@code [0 |period|)},
     * that corresponds to {@code a}.
     */
    public static double reduce(double a,
                                double period,
                                double offset) {
        final double p = FastMath.abs(period);
        return a - p * FastMath.floor((a - offset) / p) - offset;
    }

    /**
     * Returns the first argument with the sign of the second argument.
     *
     * @param magnitude Magnitude of the returned value.
     * @param sign Sign of the returned value.
     * @return a value with magnitude equal to {@code magnitude} and with the
     * same sign as the {@code sign} argument.
     * @throws MathRuntimeException if {@code magnitude == Byte.MIN_VALUE}
     * and {@code sign >= 0}.
     */
    public static byte copySign(byte magnitude, byte sign)
        throws MathRuntimeException {
        if ((magnitude >= 0 && sign >= 0) ||
            (magnitude < 0 && sign < 0)) { // Sign is OK.
            return magnitude;
        } else if (sign >= 0 &&
                   magnitude == Byte.MIN_VALUE) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW);
        } else {
            return (byte) -magnitude; // Flip sign.
        }
    }

    /**
     * Returns the first argument with the sign of the second argument.
     *
     * @param magnitude Magnitude of the returned value.
     * @param sign Sign of the returned value.
     * @return a value with magnitude equal to {@code magnitude} and with the
     * same sign as the {@code sign} argument.
     * @throws MathRuntimeException if {@code magnitude == Short.MIN_VALUE}
     * and {@code sign >= 0}.
     */
    public static short copySign(short magnitude, short sign)
            throws MathRuntimeException {
        if ((magnitude >= 0 && sign >= 0) ||
            (magnitude < 0 && sign < 0)) { // Sign is OK.
            return magnitude;
        } else if (sign >= 0 &&
                   magnitude == Short.MIN_VALUE) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW);
        } else {
            return (short) -magnitude; // Flip sign.
        }
    }

    /**
     * Returns the first argument with the sign of the second argument.
     *
     * @param magnitude Magnitude of the returned value.
     * @param sign Sign of the returned value.
     * @return a value with magnitude equal to {@code magnitude} and with the
     * same sign as the {@code sign} argument.
     * @throws MathRuntimeException if {@code magnitude == Integer.MIN_VALUE}
     * and {@code sign >= 0}.
     */
    public static int copySign(int magnitude, int sign)
            throws MathRuntimeException {
        if ((magnitude >= 0 && sign >= 0) ||
            (magnitude < 0 && sign < 0)) { // Sign is OK.
            return magnitude;
        } else if (sign >= 0 &&
                   magnitude == Integer.MIN_VALUE) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW);
        } else {
            return -magnitude; // Flip sign.
        }
    }

    /**
     * Returns the first argument with the sign of the second argument.
     *
     * @param magnitude Magnitude of the returned value.
     * @param sign Sign of the returned value.
     * @return a value with magnitude equal to {@code magnitude} and with the
     * same sign as the {@code sign} argument.
     * @throws MathRuntimeException if {@code magnitude == Long.MIN_VALUE}
     * and {@code sign >= 0}.
     */
    public static long copySign(long magnitude, long sign)
        throws MathRuntimeException {
        if ((magnitude >= 0 && sign >= 0) ||
            (magnitude < 0 && sign < 0)) { // Sign is OK.
            return magnitude;
        } else if (sign >= 0 &&
                   magnitude == Long.MIN_VALUE) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW);
        } else {
            return -magnitude; // Flip sign.
        }
    }
    /**
     * Check that the argument is a real number.
     *
     * @param x Argument.
     * @throws MathIllegalArgumentException if {@code x} is not a
     * finite real number.
     */
    public static void checkFinite(final double x)
        throws MathIllegalArgumentException {
        if (Double.isInfinite(x) || Double.isNaN(x)) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NOT_FINITE_NUMBER, x);
        }
    }

    /**
     * Check that all the elements are real numbers.
     *
     * @param val Arguments.
     * @throws MathIllegalArgumentException if any values of the array is not a
     * finite real number.
     */
    public static void checkFinite(final double[] val)
        throws MathIllegalArgumentException {
        for (final double x : val) {
            if (Double.isInfinite(x) || Double.isNaN(x)) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.NOT_FINITE_NUMBER, x);
            }
        }
    }

    /**
     * Checks that an object is not null.
     *
     * @param o Object to be checked.
     * @param pattern Message pattern.
     * @param args Arguments to replace the placeholders in {@code pattern}.
     * @throws NullArgumentException if {@code o} is {@code null}.
     */
    public static void checkNotNull(Object o,
                                    Localizable pattern,
                                    Object ... args)
        throws NullArgumentException {
        if (o == null) {
            throw new NullArgumentException(pattern, args);
        }
    }

    /**
     * Checks that an object is not null.
     *
     * @param o Object to be checked.
     * @throws NullArgumentException if {@code o} is {@code null}.
     */
    public static void checkNotNull(Object o)
        throws NullArgumentException {
        if (o == null) {
            throw new NullArgumentException(LocalizedCoreFormats.NULL_NOT_ALLOWED);
        }
    }

    /**
     * Checks that the given value is strictly within the range [lo, hi].
     *
     * @param value value to be checked.
     * @param lo the lower bound (inclusive).
     * @param hi the upper bound (inclusive).
     * @throws MathIllegalArgumentException if {@code value} is strictly outside [lo, hi].
     */
    public static void checkRangeInclusive(long value, long lo, long hi) {
        if (value < lo || value > hi) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE,
                                                   value, lo, hi);
        }
    }

    /**
     * Checks that the given value is strictly within the range [lo, hi].
     *
     * @param value value to be checked.
     * @param lo the lower bound (inclusive).
     * @param hi the upper bound (inclusive).
     * @throws MathIllegalArgumentException if {@code value} is strictly outside [lo, hi].
     */
    public static void checkRangeInclusive(double value, double lo, double hi) {
        if (value < lo || value > hi) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE,
                                                   value, lo, hi);
        }
    }

    /**
     * Checks that the given dimensions match.
     *
     * @param dimension the first dimension.
     * @param otherDimension the second dimension.
     * @throws MathIllegalArgumentException if length != otherLength.
     */
    public static void checkDimension(int dimension, int otherDimension) {
        if (dimension != otherDimension) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   dimension, otherDimension);
        }
    }

    /**
     * Sums {@code a} and {@code b} using Møller's 2Sum algorithm.
     * <p>
     * References:
     * <ul>
     * <li>Møller, Ole. "Quasi double-precision in floating point addition." BIT
     * 5, 37–50 (1965).</li>
     * <li>Shewchuk, Richard J. "Adaptive Precision Floating-Point Arithmetic
     * and Fast Robust Geometric Predicates." Discrete &amp; Computational Geometry
     * 18, 305–363 (1997).</li>
     * <li><a href=
     * "https://en.wikipedia.org/wiki/2Sum">https://en.wikipedia.org/wiki/2Sum</a></li>
     * </ul>
     * @param a first summand
     * @param b second summand
     * @return sum and residual error in the sum
     */
    public static SumAndResidual twoSum(final double a, final double b) {
        final double s = a + b;
        final double aPrime = s - b;
        final double bPrime = s - aPrime;
        final double deltaA = a - aPrime;
        final double deltaB = b - bPrime;
        final double t = deltaA + deltaB;
        return new SumAndResidual(s, t);
    }

    /**
     * Sums {@code a} and {@code b} using Møller's 2Sum algorithm.
     * <p>
     * References:
     * <ul>
     * <li>Møller, Ole. "Quasi double-precision in floating point addition." BIT
     * 5, 37–50 (1965).</li>
     * <li>Shewchuk, Richard J. "Adaptive Precision Floating-Point Arithmetic
     * and Fast Robust Geometric Predicates." Discrete &amp; Computational Geometry
     * 18, 305–363 (1997).</li>
     * <li><a href=
     * "https://en.wikipedia.org/wiki/2Sum">https://en.wikipedia.org/wiki/2Sum</a></li>
     * </ul>
     * @param <T> field element type
     * @param a first summand
     * @param b second summand
     * @return sum and residual error in the sum
     */
    public static <T extends FieldElement<T>> FieldSumAndResidual<T> twoSum(final T a, final T b) {
        final T s = a.add(b);
        final T aPrime = s.subtract(b);
        final T bPrime = s.subtract(aPrime);
        final T deltaA = a.subtract(aPrime);
        final T deltaB = b.subtract(bPrime);
        final T t = deltaA.add(deltaB);
        return new FieldSumAndResidual<>(s, t);
    }

    /**
     * Result class for {@link MathUtils#twoSum(double, double)} containing the
     * sum and the residual error in the sum.
     */
    public static final class SumAndResidual {

        /** Sum. */
        private final double sum;
        /** Residual error in the sum. */
        private final double residual;

        /**
         * Constructs a {@link SumAndResidual} instance.
         * @param sum sum
         * @param residual residual error in the sum
         */
        private SumAndResidual(final double sum, final double residual) {
            this.sum = sum;
            this.residual = residual;
        }

        /**
         * Returns the sum.
         * @return sum
         */
        public double getSum() {
            return sum;
        }

        /**
         * Returns the residual error in the sum.
         * @return residual error in the sum
         */
        public double getResidual() {
            return residual;
        }

    }

    /**
     * Result class for
     * {@link MathUtils#twoSum(FieldElement, FieldElement)} containing
     * the sum and the residual error in the sum.
     * @param <T> field element type
     */
    public static final class FieldSumAndResidual<T extends FieldElement<T>> {

        /** Sum. */
        private final T sum;
        /** Residual error in the sum. */
        private final T residual;

        /**
         * Constructs a {@link FieldSumAndResidual} instance.
         * @param sum sum
         * @param residual residual error in the sum
         */
        private FieldSumAndResidual(final T sum, final T residual) {
            this.sum = sum;
            this.residual = residual;
        }

        /**
         * Returns the sum.
         * @return sum
         */
        public T getSum() {
            return sum;
        }

        /**
         * Returns the residual error in the sum.
         * @return residual error in the sum
         */
        public T getResidual() {
            return residual;
        }

    }

}
