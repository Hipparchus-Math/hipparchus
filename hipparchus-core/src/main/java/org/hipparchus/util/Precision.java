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

import java.math.BigDecimal;
import java.math.RoundingMode;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * Utilities for comparing numbers.
 */
public class Precision {
    /**
     * Largest double-precision floating-point number such that
     * {@code 1 + EPSILON} is numerically equal to 1. This value is an upper
     * bound on the relative error due to rounding real numbers to double
     * precision floating-point numbers.
     * <p>
     * In IEEE 754 arithmetic, this is 2<sup>-53</sup>.
     *
     * @see <a href="http://en.wikipedia.org/wiki/Machine_epsilon">Machine epsilon</a>
     */
    public static final double EPSILON;

    /**
     * Safe minimum, such that {@code 1 / SAFE_MIN} does not overflow.
     * <br/>
     * In IEEE 754 arithmetic, this is also the smallest normalized
     * number 2<sup>-1022</sup>.
     */
    public static final double SAFE_MIN;

    /** Exponent offset in IEEE754 representation. */
    private static final long EXPONENT_OFFSET = 1023l;

    /** Offset to order signed double numbers lexicographically. */
    private static final long SGN_MASK = 0x8000000000000000L;
    /** Offset to order signed double numbers lexicographically. */
    private static final int SGN_MASK_FLOAT = 0x80000000;
    /** Positive zero. */
    private static final double POSITIVE_ZERO = 0d;
    /** Positive zero bits. */
    private static final long POSITIVE_ZERO_DOUBLE_BITS = Double.doubleToRawLongBits(+0.0);
    /** Negative zero bits. */
    private static final long NEGATIVE_ZERO_DOUBLE_BITS = Double.doubleToRawLongBits(-0.0);
    /** Positive zero bits. */
    private static final int POSITIVE_ZERO_FLOAT_BITS   = Float.floatToRawIntBits(+0.0f);
    /** Negative zero bits. */
    private static final int NEGATIVE_ZERO_FLOAT_BITS   = Float.floatToRawIntBits(-0.0f);
    /** Mask used to extract exponent from double bits. */
    private static final long MASK_DOUBLE_EXPONENT = 0x7ff0000000000000L;
    /** Mask used to extract mantissa from double bits. */
    private static final long MASK_DOUBLE_MANTISSA = 0x000fffffffffffffL;
    /** Mask used to add implicit high order bit for normalized double. */
    private static final long IMPLICIT_DOUBLE_HIGH_BIT = 0x0010000000000000L;
    /** Mask used to extract exponent from float bits. */
    private static final int MASK_FLOAT_EXPONENT = 0x7f800000;
    /** Mask used to extract mantissa from float bits. */
    private static final int MASK_FLOAT_MANTISSA = 0x007fffff;
    /** Mask used to add implicit high order bit for normalized float. */
    private static final int IMPLICIT_FLOAT_HIGH_BIT = 0x00800000;

    static {
        /*
         *  This was previously expressed as = 0x1.0p-53;
         *  However, OpenJDK (Sparc Solaris) cannot handle such small
         *  constants: MATH-721
         */
        EPSILON = Double.longBitsToDouble((EXPONENT_OFFSET - 53l) << 52);

        /*
         * This was previously expressed as = 0x1.0p-1022;
         * However, OpenJDK (Sparc Solaris) cannot handle such small
         * constants: MATH-721
         */
        SAFE_MIN = Double.longBitsToDouble((EXPONENT_OFFSET - 1022l) << 52);
    }

    /**
     * Private constructor.
     */
    private Precision() {}

    /**
     * Compares two numbers given some amount of allowed error.
     *
     * @param x the first number
     * @param y the second number
     * @param eps the amount of error to allow when checking for equality
     * @return <ul><li>0 if  {@link #equals(double, double, double) equals(x, y, eps)}</li>
     *       <li>&lt; 0 if !{@link #equals(double, double, double) equals(x, y, eps)} &amp;&amp; x &lt; y</li>
     *       <li>&gt; 0 if !{@link #equals(double, double, double) equals(x, y, eps)} &amp;&amp; x &gt; y or
     *       either argument is NaN</li></ul>
     */
    public static int compareTo(double x, double y, double eps) {
        if (equals(x, y, eps)) {
            return 0;
        } else if (x < y) {
            return -1;
        }
        return 1;
    }

    /**
     * Compares two numbers given some amount of allowed error.
     * Two float numbers are considered equal if there are {@code (maxUlps - 1)}
     * (or fewer) floating point numbers between them, i.e. two adjacent floating
     * point numbers are considered equal.
     * Adapted from <a
     * href="http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/">
     * Bruce Dawson</a>. Returns {@code false} if either of the arguments is NaN.
     *
     * @param x first value
     * @param y second value
     * @param maxUlps {@code (maxUlps - 1)} is the number of floating point
     * values between {@code x} and {@code y}.
     * @return <ul><li>0 if  {@link #equals(double, double, int) equals(x, y, maxUlps)}</li>
     *       <li>&lt; 0 if !{@link #equals(double, double, int) equals(x, y, maxUlps)} &amp;&amp; x &lt; y</li>
     *       <li>&gt; 0 if !{@link #equals(double, double, int) equals(x, y, maxUlps)} &amp;&amp; x &gt; y
     *       or either argument is NaN</li></ul>
     */
    public static int compareTo(final double x, final double y, final int maxUlps) {
        if (equals(x, y, maxUlps)) {
            return 0;
        } else if (x < y) {
            return -1;
        }
        return 1;
    }

    /**
     * Returns true iff they are equal as defined by
     * {@link #equals(float,float,int) equals(x, y, 1)}.
     *
     * @param x first value
     * @param y second value
     * @return {@code true} if the values are equal.
     */
    public static boolean equals(float x, float y) {
        return equals(x, y, 1);
    }

    /**
     * Returns true if both arguments are NaN or they are
     * equal as defined by {@link #equals(float,float) equals(x, y, 1)}.
     *
     * @param x first value
     * @param y second value
     * @return {@code true} if the values are equal or both are NaN.
     */
    public static boolean equalsIncludingNaN(float x, float y) {
        return (x != x || y != y) ? !(x != x ^ y != y) : equals(x, y, 1);
    }

    /**
     * Returns true if the arguments are equal or within the range of allowed
     * error (inclusive).  Returns {@code false} if either of the arguments
     * is NaN.
     *
     * @param x first value
     * @param y second value
     * @param eps the amount of absolute error to allow.
     * @return {@code true} if the values are equal or within range of each other.
     */
    public static boolean equals(float x, float y, float eps) {
        return equals(x, y, 1) || FastMath.abs(y - x) <= eps;
    }

    /**
     * Returns true if the arguments are both NaN, are equal, or are within the range
     * of allowed error (inclusive).
     *
     * @param x first value
     * @param y second value
     * @param eps the amount of absolute error to allow.
     * @return {@code true} if the values are equal or within range of each other,
     * or both are NaN.
     */
    public static boolean equalsIncludingNaN(float x, float y, float eps) {
        return equalsIncludingNaN(x, y) || (FastMath.abs(y - x) <= eps);
    }

    /**
     * Returns true if the arguments are equal or within the range of allowed
     * error (inclusive).
     * Two float numbers are considered equal if there are {@code (maxUlps - 1)}
     * (or fewer) floating point numbers between them, i.e. two adjacent floating
     * point numbers are considered equal.
     * Adapted from <a
     * href="http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/">
     * Bruce Dawson</a>.  Returns {@code false} if either of the arguments is NaN.
     *
     * @param x first value
     * @param y second value
     * @param maxUlps {@code (maxUlps - 1)} is the number of floating point
     * values between {@code x} and {@code y}.
     * @return {@code true} if there are fewer than {@code maxUlps} floating
     * point values between {@code x} and {@code y}.
     */
    public static boolean equals(final float x, final float y, final int maxUlps) {

        final int xInt = Float.floatToRawIntBits(x);
        final int yInt = Float.floatToRawIntBits(y);

        final boolean isEqual;
        if (((xInt ^ yInt) & SGN_MASK_FLOAT) == 0) {
            // number have same sign, there is no risk of overflow
            isEqual = FastMath.abs(xInt - yInt) <= maxUlps;
        } else {
            // number have opposite signs, take care of overflow
            final int deltaPlus;
            final int deltaMinus;
            if (xInt < yInt) {
                deltaPlus  = yInt - POSITIVE_ZERO_FLOAT_BITS;
                deltaMinus = xInt - NEGATIVE_ZERO_FLOAT_BITS;
            } else {
                deltaPlus  = xInt - POSITIVE_ZERO_FLOAT_BITS;
                deltaMinus = yInt - NEGATIVE_ZERO_FLOAT_BITS;
            }

            if (deltaPlus > maxUlps) {
                isEqual = false;
            } else {
                isEqual = deltaMinus <= (maxUlps - deltaPlus);
            }

        }

        return isEqual && !Float.isNaN(x) && !Float.isNaN(y);

    }

    /**
     * Returns true if the arguments are both NaN or if they are equal as defined
     * by {@link #equals(float,float,int) equals(x, y, maxUlps)}.
     *
     * @param x first value
     * @param y second value
     * @param maxUlps {@code (maxUlps - 1)} is the number of floating point
     * values between {@code x} and {@code y}.
     * @return {@code true} if both arguments are NaN or if there are less than
     * {@code maxUlps} floating point values between {@code x} and {@code y}.
     */
    public static boolean equalsIncludingNaN(float x, float y, int maxUlps) {
        return (x != x || y != y) ? !(x != x ^ y != y) : equals(x, y, maxUlps);
    }

    /**
     * Returns true iff they are equal as defined by
     * {@link #equals(double,double,int) equals(x, y, 1)}.
     *
     * @param x first value
     * @param y second value
     * @return {@code true} if the values are equal.
     */
    public static boolean equals(double x, double y) {
        return equals(x, y, 1);
    }

    /**
     * Returns true if the arguments are both NaN or they are
     * equal as defined by {@link #equals(double,double) equals(x, y, 1)}.
     *
     * @param x first value
     * @param y second value
     * @return {@code true} if the values are equal or both are NaN.
     */
    public static boolean equalsIncludingNaN(double x, double y) {
        return (x != x || y != y) ? !(x != x ^ y != y) : equals(x, y, 1);
    }

    /**
     * Returns {@code true} if there is no double value strictly between the
     * arguments or the difference between them is within the range of allowed
     * error (inclusive). Returns {@code false} if either of the arguments
     * is NaN.
     *
     * @param x First value.
     * @param y Second value.
     * @param eps Amount of allowed absolute error.
     * @return {@code true} if the values are two adjacent floating point
     * numbers or they are within range of each other.
     */
    public static boolean equals(double x, double y, double eps) {
        return equals(x, y, 1) || FastMath.abs(y - x) <= eps;
    }

    /**
     * Returns {@code true} if there is no double value strictly between the
     * arguments or the relative difference between them is less than or equal
     * to the given tolerance. Returns {@code false} if either of the arguments
     * is NaN.
     *
     * @param x First value.
     * @param y Second value.
     * @param eps Amount of allowed relative error.
     * @return {@code true} if the values are two adjacent floating point
     * numbers or they are within range of each other.
     */
    public static boolean equalsWithRelativeTolerance(double x, double y, double eps) {
        if (equals(x, y, 1)) {
            return true;
        }

        final double absoluteMax = FastMath.max(FastMath.abs(x), FastMath.abs(y));
        final double relativeDifference = FastMath.abs((x - y) / absoluteMax);

        return relativeDifference <= eps;
    }

    /**
     * Returns true if the arguments are both NaN, are equal or are within the range
     * of allowed error (inclusive).
     *
     * @param x first value
     * @param y second value
     * @param eps the amount of absolute error to allow.
     * @return {@code true} if the values are equal or within range of each other,
     * or both are NaN.
     */
    public static boolean equalsIncludingNaN(double x, double y, double eps) {
        return equalsIncludingNaN(x, y) || (FastMath.abs(y - x) <= eps);
    }

    /**
     * Returns true if the arguments are equal or within the range of allowed
     * error (inclusive).
     * <p>
     * Two float numbers are considered equal if there are {@code (maxUlps - 1)}
     * (or fewer) floating point numbers between them, i.e. two adjacent
     * floating point numbers are considered equal.
     * </p>
     * <p>
     * Adapted from <a
     * href="http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/">
     * Bruce Dawson</a>. Returns {@code false} if either of the arguments is NaN.
     * </p>
     *
     * @param x first value
     * @param y second value
     * @param maxUlps {@code (maxUlps - 1)} is the number of floating point
     * values between {@code x} and {@code y}.
     * @return {@code true} if there are fewer than {@code maxUlps} floating
     * point values between {@code x} and {@code y}.
     */
    public static boolean equals(final double x, final double y, final int maxUlps) {

        final long xInt = Double.doubleToRawLongBits(x);
        final long yInt = Double.doubleToRawLongBits(y);

        final boolean isEqual;
        if (((xInt ^ yInt) & SGN_MASK) == 0l) {
            // number have same sign, there is no risk of overflow
            isEqual = FastMath.abs(xInt - yInt) <= maxUlps;
        } else {
            // number have opposite signs, take care of overflow
            final long deltaPlus;
            final long deltaMinus;
            if (xInt < yInt) {
                deltaPlus  = yInt - POSITIVE_ZERO_DOUBLE_BITS;
                deltaMinus = xInt - NEGATIVE_ZERO_DOUBLE_BITS;
            } else {
                deltaPlus  = xInt - POSITIVE_ZERO_DOUBLE_BITS;
                deltaMinus = yInt - NEGATIVE_ZERO_DOUBLE_BITS;
            }

            if (deltaPlus > maxUlps) {
                isEqual = false;
            } else {
                isEqual = deltaMinus <= (maxUlps - deltaPlus);
            }

        }

        return isEqual && !Double.isNaN(x) && !Double.isNaN(y);

    }

    /**
     * Returns true if both arguments are NaN or if they are equal as defined
     * by {@link #equals(double,double,int) equals(x, y, maxUlps)}.
     *
     * @param x first value
     * @param y second value
     * @param maxUlps {@code (maxUlps - 1)} is the number of floating point
     * values between {@code x} and {@code y}.
     * @return {@code true} if both arguments are NaN or if there are less than
     * {@code maxUlps} floating point values between {@code x} and {@code y}.
     */
    public static boolean equalsIncludingNaN(double x, double y, int maxUlps) {
        return (x != x || y != y) ? !(x != x ^ y != y) : equals(x, y, maxUlps);
    }

    /**
     * Rounds the given value to the specified number of decimal places.
     * The value is rounded using the {@link BigDecimal#ROUND_HALF_UP} method.
     *
     * @param x Value to round.
     * @param scale Number of digits to the right of the decimal point.
     * @return the rounded value.
     */
    public static double round(double x, int scale) {
        return round(x, scale, RoundingMode.HALF_UP);
    }

    /**
     * Rounds the given value to the specified number of decimal places.
     * The value is rounded using the given method which is any method defined
     * in {@link BigDecimal}.
     * If {@code x} is infinite or {@code NaN}, then the value of {@code x} is
     * returned unchanged, regardless of the other parameters.
     *
     * @param x Value to round.
     * @param scale Number of digits to the right of the decimal point.
     * @param roundingMethod Rounding method as defined in {@link BigDecimal}.
     * @return the rounded value.
     * @throws ArithmeticException if {@code roundingMethod == ROUND_UNNECESSARY}
     * and the specified scaling operation would require rounding.
     * @throws IllegalArgumentException if {@code roundingMethod} does not
     * represent a valid rounding mode.
     */
    public static double round(double x, int scale, RoundingMode roundingMethod) {
        try {
            final double rounded = (new BigDecimal(Double.toString(x))
                   .setScale(scale, roundingMethod))
                   .doubleValue();
            // MATH-1089: negative values rounded to zero should result in negative zero
            return rounded == POSITIVE_ZERO ? POSITIVE_ZERO * x : rounded;
        } catch (NumberFormatException ex) {
            if (Double.isInfinite(x)) {
                return x;
            } else {
                return Double.NaN;
            }
        }
    }

    /**
     * Rounds the given value to the specified number of decimal places.
     * The value is rounded using the {@link BigDecimal#ROUND_HALF_UP} method.
     *
     * @param x Value to round.
     * @param scale Number of digits to the right of the decimal point.
     * @return the rounded value.
     */
    public static float round(float x, int scale) {
        return round(x, scale, RoundingMode.HALF_UP);
    }

    /**
     * Rounds the given value to the specified number of decimal places.
     * The value is rounded using the given method which is any method defined
     * in {@link BigDecimal}.
     *
     * @param x Value to round.
     * @param scale Number of digits to the right of the decimal point.
     * @param roundingMethod Rounding method as defined in {@link BigDecimal}.
     * @return the rounded value.
     * @throws MathRuntimeException if an exact operation is required but result is not exact
     * @throws MathIllegalArgumentException if {@code roundingMethod} is not a valid rounding method.
     */
    public static float round(float x, int scale, RoundingMode roundingMethod)
        throws MathRuntimeException, MathIllegalArgumentException {
        final float sign = FastMath.copySign(1f, x);
        final float factor = (float) FastMath.pow(10.0f, scale) * sign;
        return (float) roundUnscaled(x * factor, sign, roundingMethod) / factor;
    }

    /**
     * Rounds the given non-negative value to the "nearest" integer. Nearest is
     * determined by the rounding method specified. Rounding methods are defined
     * in {@link BigDecimal}.
     *
     * @param unscaled Value to round.
     * @param sign Sign of the original, scaled value.
     * @param roundingMethod Rounding method, as defined in {@link BigDecimal}.
     * @return the rounded value.
     * @throws MathRuntimeException if an exact operation is required but result is not exact
     * @throws MathIllegalArgumentException if {@code roundingMethod} is not a valid rounding method.
     */
    private static double roundUnscaled(double unscaled,
                                        double sign,
                                        RoundingMode roundingMethod)
        throws MathRuntimeException, MathIllegalArgumentException {
        switch (roundingMethod) {
        case CEILING :
            if (sign == -1) {
                unscaled = FastMath.floor(FastMath.nextAfter(unscaled, Double.NEGATIVE_INFINITY));
            } else {
                unscaled = FastMath.ceil(FastMath.nextAfter(unscaled, Double.POSITIVE_INFINITY));
            }
            break;
        case DOWN :
            unscaled = FastMath.floor(FastMath.nextAfter(unscaled, Double.NEGATIVE_INFINITY));
            break;
        case FLOOR :
            if (sign == -1) {
                unscaled = FastMath.ceil(FastMath.nextAfter(unscaled, Double.POSITIVE_INFINITY));
            } else {
                unscaled = FastMath.floor(FastMath.nextAfter(unscaled, Double.NEGATIVE_INFINITY));
            }
            break;
        case HALF_DOWN : {
            unscaled = FastMath.nextAfter(unscaled, Double.NEGATIVE_INFINITY);
            double fraction = unscaled - FastMath.floor(unscaled);
            if (fraction > 0.5) {
                unscaled = FastMath.ceil(unscaled);
            } else {
                unscaled = FastMath.floor(unscaled);
            }
            break;
        }
        case HALF_EVEN : {
            double fraction = unscaled - FastMath.floor(unscaled);
            if (fraction > 0.5) {
                unscaled = FastMath.ceil(unscaled);
            } else if (fraction < 0.5) {
                unscaled = FastMath.floor(unscaled);
            } else {
                // The following equality test is intentional and needed for rounding purposes
                if (FastMath.floor(unscaled) / 2.0 == FastMath.floor(FastMath.floor(unscaled) / 2.0)) { // even
                    unscaled = FastMath.floor(unscaled);
                } else { // odd
                    unscaled = FastMath.ceil(unscaled);
                }
            }
            break;
        }
        case HALF_UP : {
            unscaled = FastMath.nextAfter(unscaled, Double.POSITIVE_INFINITY);
            double fraction = unscaled - FastMath.floor(unscaled);
            if (fraction >= 0.5) {
                unscaled = FastMath.ceil(unscaled);
            } else {
                unscaled = FastMath.floor(unscaled);
            }
            break;
        }
        case UNNECESSARY :
            if (unscaled != FastMath.floor(unscaled)) {
                throw new MathRuntimeException(LocalizedCoreFormats.ARITHMETIC_EXCEPTION);
            }
            break;
        case UP :
            // do not round if the discarded fraction is equal to zero
            if (unscaled != FastMath.floor(unscaled)) {
                unscaled = FastMath.ceil(FastMath.nextAfter(unscaled, Double.POSITIVE_INFINITY));
            }
            break;
        default :
            // this should nerver happen
            throw MathRuntimeException.createInternalError();
        }
        return unscaled;
    }

    /** Check is x is a mathematical integer.
     * @param x number to check
     * @return true if x is a mathematical integer
     * @since 1.7
     */
    public static boolean isMathematicalInteger(final double x) {
        final long bits   = Double.doubleToRawLongBits(x);
        final int  rawExp = (int) ((bits & MASK_DOUBLE_EXPONENT) >> 52);
        if (rawExp == 2047) {
            // NaN or infinite
            return false;
        } else {
            // a double that may have a fractional part
            final long rawMantissa    = bits & MASK_DOUBLE_MANTISSA;
            final long fullMantissa   = rawExp > 0 ? (IMPLICIT_DOUBLE_HIGH_BIT | rawMantissa) : rawMantissa;
            final long fractionalMask = (IMPLICIT_DOUBLE_HIGH_BIT | MASK_DOUBLE_MANTISSA) >> FastMath.min(53, FastMath.max(0, rawExp - 1022));
            return (fullMantissa & fractionalMask) == 0l;
        }
    }

    /** Check is x is a mathematical integer.
     * @param x number to check
     * @return true if x is a mathematical integer
     * @since 1.7
     */
    public static boolean isMathematicalInteger(final float x) {
        final int bits   = Float.floatToRawIntBits(x);
        final int rawExp = (int) ((bits & MASK_FLOAT_EXPONENT) >> 23);
        if (rawExp == 255) {
            // NaN or infinite
            return false;
        } else {
            // a float that may have a fractional part
            final int rawMantissa    = bits & MASK_FLOAT_MANTISSA;
            final int fullMantissa   = rawExp > 0 ? (IMPLICIT_FLOAT_HIGH_BIT | rawMantissa) : rawMantissa;
            final int fractionalMask = (IMPLICIT_FLOAT_HIGH_BIT | MASK_FLOAT_MANTISSA) >> FastMath.min(24, FastMath.max(0, rawExp - 126));
            return (fullMantissa & fractionalMask) == 0;
        }
    }

    /**
     * Computes a number {@code delta} close to {@code originalDelta} with
     * the property that <pre><code>
     *   x + delta - x
     * </code></pre>
     * is exactly machine-representable.
     * This is useful when computing numerical derivatives, in order to reduce
     * roundoff errors.
     *
     * @param x Value.
     * @param originalDelta Offset value.
     * @return a number {@code delta} so that {@code x + delta} and {@code x}
     * differ by a representable floating number.
     */
    public static double representableDelta(double x,
                                            double originalDelta) {
        return x + originalDelta - x;
    }
}
