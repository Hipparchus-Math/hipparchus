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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * This class wraps a {@code double} value in an object. It is similar to the
 * standard class {@link Double}, while also implementing the
 * {@link CalculusFieldElement} interface.
 */
public class Binary64 extends Number implements CalculusFieldElement<Binary64>, Comparable<Binary64> {

    /** The constant value of {@code 0d} as a {@code Binary64}. */
    public static final Binary64 ZERO;

    /** The constant value of {@code 1d} as a {@code Binary64}. */
    public static final Binary64 ONE;

    /** The constant value of π as a {@code Binary64}. */
    public static final Binary64 PI;

    /**
     * The constant value of {@link Double#NEGATIVE_INFINITY} as a
     * {@code Binary64}.
     */
    public static final Binary64 NEGATIVE_INFINITY;

    /**
     * The constant value of {@link Double#POSITIVE_INFINITY} as a
     * {@code Binary64}.
     */
    public static final Binary64 POSITIVE_INFINITY;

    /** The constant value of {@link Double#NaN} as a {@code Binary64}. */
    public static final Binary64 NAN;

    /** */
    private static final long serialVersionUID = 20120227L;

    static {
        ZERO = new Binary64(0d);
        ONE  = new Binary64(1d);
        PI   = new Binary64(FastMath.PI);
        NEGATIVE_INFINITY = new Binary64(Double.NEGATIVE_INFINITY);
        POSITIVE_INFINITY = new Binary64(Double.POSITIVE_INFINITY);
        NAN = new Binary64(Double.NaN);
    }

    /** The primitive {@code double} value of this object. */
    private final double value;

    /**
     * Creates a new instance of this class.
     *
     * @param x the primitive {@code double} value of the object to be created
     */
    public Binary64(final double x) {
        this.value = x;
    }

    /*
     * Methods from the FieldElement interface.
     */

    /** {@inheritDoc} */
    @Override
    public Binary64 newInstance(final double v) {
        return new Binary64(v);
    }

    /** {@inheritDoc} */
    @Override
    public Field<Binary64> getField() {
        return Binary64Field.getInstance();
    }

    /**
     * {@inheritDoc}
     *
     * The current implementation strictly enforces
     * {@code this.add(a).equals(new Binary64(this.doubleValue()
     * + a.doubleValue()))}.
     */
    @Override
    public Binary64 add(final Binary64 a) {
        return new Binary64(this.value + a.value);
    }

    /**
     * {@inheritDoc}
     *
     * The current implementation strictly enforces
     * {@code this.subtract(a).equals(new Binary64(this.doubleValue()
     * - a.doubleValue()))}.
     */
    @Override
    public Binary64 subtract(final Binary64 a) {
        return new Binary64(this.value - a.value);
    }

    /**
     * {@inheritDoc}
     *
     * The current implementation strictly enforces
     * {@code this.negate().equals(new Binary64(-this.doubleValue()))}.
     */
    @Override
    public Binary64 negate() {
        return new Binary64(-this.value);
    }

    /**
     * {@inheritDoc}
     *
     * The current implementation strictly enforces
     * {@code this.multiply(a).equals(new Binary64(this.doubleValue()
     * * a.doubleValue()))}.
     */
    @Override
    public Binary64 multiply(final Binary64 a) {
        return new Binary64(this.value * a.value);
    }

    /**
     * {@inheritDoc}
     *
     * The current implementation strictly enforces
     * {@code this.multiply(n).equals(new Binary64(n * this.doubleValue()))}.
     */
    @Override
    public Binary64 multiply(final int n) {
        return new Binary64(n * this.value);
    }

    /**
     * {@inheritDoc}
     *
     * The current implementation strictly enforces
     * {@code this.divide(a).equals(new Binary64(this.doubleValue()
     * / a.doubleValue()))}.
     *
     */
    @Override
    public Binary64 divide(final Binary64 a) {
        return new Binary64(this.value / a.value);
    }

    /**
     * {@inheritDoc}
     *
     * The current implementation strictly enforces
     * {@code this.reciprocal().equals(new Binary64(1.0
     * / this.doubleValue()))}.
     */
    @Override
    public Binary64 reciprocal() {
        return new Binary64(1.0 / this.value);
    }

    /*
     * Methods from the Number abstract class
     */

    /**
     * {@inheritDoc}
     *
     * The current implementation performs casting to a {@code byte}.
     */
    @Override
    public byte byteValue() {
        return (byte) value;
    }

    /**
     * {@inheritDoc}
     *
     * The current implementation performs casting to a {@code short}.
     */
    @Override
    public short shortValue() {
        return (short) value;
    }

    /**
     * {@inheritDoc}
     *
     * The current implementation performs casting to a {@code int}.
     */
    @Override
    public int intValue() {
        return (int) value;
    }

    /**
     * {@inheritDoc}
     *
     * The current implementation performs casting to a {@code long}.
     */
    @Override
    public long longValue() {
        return (long) value;
    }

    /**
     * {@inheritDoc}
     *
     * The current implementation performs casting to a {@code float}.
     */
    @Override
    public float floatValue() {
        return (float) value;
    }

    /** {@inheritDoc} */
    @Override
    public double doubleValue() {
        return value;
    }

    /*
     * Methods from the Comparable interface.
     */

    /**
     * {@inheritDoc}
     *
     * The current implementation returns the same value as
     * {@code new Double(this.doubleValue()).compareTo(new
     * Double(o.doubleValue()))}
     *
     * @see Double#compareTo(Double)
     */
    @Override
    public int compareTo(final Binary64 o) {
        return Double.compare(this.value, o.value);
    }

    /*
     * Methods from the Object abstract class.
     */

    /** {@inheritDoc} */
    @Override
    public boolean equals(final Object obj) {
        if (obj instanceof Binary64) {
            final Binary64 that = (Binary64) obj;
            return Double.doubleToLongBits(this.value) == Double
                    .doubleToLongBits(that.value);
        }
        return false;
    }

    /** {@inheritDoc}
     * <p>
     * This implementation considers +0.0 and -0.0 to be equal.
     * </p>
     * @since 1.8
     */
    @Override
    public boolean isZero() {
        return value == 0.0;
    }

    /**
     * {@inheritDoc}
     *
     * The current implementation returns the same value as
     * {@code new Double(this.doubleValue()).hashCode()}
     *
     * @see Double#hashCode()
     */
    @Override
    public int hashCode() {
        long v = Double.doubleToLongBits(value);
        return (int) (v ^ (v >>> 32));
    }

    /**
     * {@inheritDoc}
     *
     * The returned {@code String} is equal to
     * {@code Double.toString(this.doubleValue())}
     *
     * @see Double#toString(double)
     */
    @Override
    public String toString() {
        return Double.toString(value);
    }

    /*
     * Methods inspired by the Double class.
     */

    /**
     * Returns {@code true} if {@code this} double precision number is infinite
     * ({@link Double#POSITIVE_INFINITY} or {@link Double#NEGATIVE_INFINITY}).
     *
     * @return {@code true} if {@code this} number is infinite
     */
    @Override
    public boolean isInfinite() {
        return Double.isInfinite(value);
    }

    /**
     * Returns {@code true} if {@code this} double precision number is
     * Not-a-Number ({@code NaN}), false otherwise.
     *
     * @return {@code true} if {@code this} is {@code NaN}
     */
    @Override
    public boolean isNaN() {
        return Double.isNaN(value);
    }

    /** {@inheritDoc} */
    @Override
    public double getReal() {
        return value;
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 add(final double a) {
        return new Binary64(value + a);
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 subtract(final double a) {
        return new Binary64(value - a);
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 multiply(final double a) {
        return new Binary64(value * a);
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 divide(final double a) {
        return new Binary64(value / a);
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 remainder(final double a) {
        return new Binary64(FastMath.IEEEremainder(value, a));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 remainder(final Binary64 a) {
        return new Binary64(FastMath.IEEEremainder(value, a.value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 abs() {
        return new Binary64(FastMath.abs(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 ceil() {
        return new Binary64(FastMath.ceil(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 floor() {
        return new Binary64(FastMath.floor(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 rint() {
        return new Binary64(FastMath.rint(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 sign() {
        return new Binary64(FastMath.signum(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 copySign(final Binary64 sign) {
        return new Binary64(FastMath.copySign(value, sign.value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 copySign(final double sign) {
        return new Binary64(FastMath.copySign(value, sign));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 scalb(final int n) {
        return new Binary64(FastMath.scalb(value, n));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 ulp() {
        return new Binary64(FastMath.ulp(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 hypot(final Binary64 y) {
        return new Binary64(FastMath.hypot(value, y.value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 sqrt() {
        return new Binary64(FastMath.sqrt(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 cbrt() {
        return new Binary64(FastMath.cbrt(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 rootN(final int n) {
        if (value < 0) {
            return (n % 2 == 0) ? Binary64.NAN : new Binary64(-FastMath.pow(-value, 1.0 / n));
        } else {
            return new Binary64(FastMath.pow(value, 1.0 / n));
        }
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 pow(final double p) {
        return new Binary64(FastMath.pow(value, p));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 pow(final int n) {
        return new Binary64(FastMath.pow(value, n));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 pow(final Binary64 e) {
        return new Binary64(FastMath.pow(value, e.value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 exp() {
        return new Binary64(FastMath.exp(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 expm1() {
        return new Binary64(FastMath.expm1(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 log() {
        return new Binary64(FastMath.log(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 log1p() {
        return new Binary64(FastMath.log1p(value));
    }

    /** Base 10 logarithm.
     * @return base 10 logarithm of the instance
     */
    @Override
    public Binary64 log10() {
        return new Binary64(FastMath.log10(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 cos() {
        return new Binary64(FastMath.cos(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 sin() {
        return new Binary64(FastMath.sin(value));
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinCos<Binary64> sinCos() {
        final SinCos sc = FastMath.sinCos(value);
        return new FieldSinCos<>(new Binary64(sc.sin()), new Binary64(sc.cos()));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 tan() {
        return new Binary64(FastMath.tan(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 acos() {
        return new Binary64(FastMath.acos(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 asin() {
        return new Binary64(FastMath.asin(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 atan() {
        return new Binary64(FastMath.atan(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 atan2(final Binary64 x) {
        return new Binary64(FastMath.atan2(value, x.value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 cosh() {
        return new Binary64(FastMath.cosh(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 sinh() {
        return new Binary64(FastMath.sinh(value));
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinhCosh<Binary64> sinhCosh() {
        final SinhCosh sch = FastMath.sinhCosh(value);
        return new FieldSinhCosh<>(new Binary64(sch.sinh()), new Binary64(sch.cosh()));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 tanh() {
        return new Binary64(FastMath.tanh(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 acosh() {
        return new Binary64(FastMath.acosh(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 asinh() {
        return new Binary64(FastMath.asinh(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 atanh() {
        return new Binary64(FastMath.atanh(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 toDegrees() {
        return new Binary64(FastMath.toDegrees(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 toRadians() {
        return new Binary64(FastMath.toRadians(value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 linearCombination(final Binary64[] a, final Binary64[] b)
        throws MathIllegalArgumentException {
        MathUtils.checkDimension(a.length, b.length);
        final double[] aDouble = new double[a.length];
        final double[] bDouble = new double[b.length];
        for (int i = 0; i < a.length; ++i) {
            aDouble[i] = a[i].value;
            bDouble[i] = b[i].value;
        }
        return new Binary64(MathArrays.linearCombination(aDouble, bDouble));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 linearCombination(final double[] a, final Binary64[] b)
        throws MathIllegalArgumentException {
        MathUtils.checkDimension(a.length, b.length);
        final double[] bDouble = new double[b.length];
        for (int i = 0; i < a.length; ++i) {
            bDouble[i] = b[i].value;
        }
        return new Binary64(MathArrays.linearCombination(a, bDouble));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 linearCombination(final Binary64 a1, final Binary64 b1,
                                       final Binary64 a2, final Binary64 b2) {
        return new Binary64(MathArrays.linearCombination(a1.value, b1.value,
                                                          a2.value, b2.value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 linearCombination(final double a1, final Binary64 b1,
                                       final double a2, final Binary64 b2) {
        return new Binary64(MathArrays.linearCombination(a1, b1.value,
                                                          a2, b2.value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 linearCombination(final Binary64 a1, final Binary64 b1,
                                       final Binary64 a2, final Binary64 b2,
                                       final Binary64 a3, final Binary64 b3) {
        return new Binary64(MathArrays.linearCombination(a1.value, b1.value,
                                                          a2.value, b2.value,
                                                          a3.value, b3.value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 linearCombination(final double a1, final Binary64 b1,
                                       final double a2, final Binary64 b2,
                                       final double a3, final Binary64 b3) {
        return new Binary64(MathArrays.linearCombination(a1, b1.value,
                                                          a2, b2.value,
                                                          a3, b3.value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 linearCombination(final Binary64 a1, final Binary64 b1,
                                       final Binary64 a2, final Binary64 b2,
                                       final Binary64 a3, final Binary64 b3,
                                       final Binary64 a4, final Binary64 b4) {
        return new Binary64(MathArrays.linearCombination(a1.value, b1.value,
                                                          a2.value, b2.value,
                                                          a3.value, b3.value,
                                                          a4.value, b4.value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 linearCombination(final double a1, final Binary64 b1,
                                       final double a2, final Binary64 b2,
                                       final double a3, final Binary64 b3,
                                       final double a4, final Binary64 b4) {
        return new Binary64(MathArrays.linearCombination(a1, b1.value,
                                                          a2, b2.value,
                                                          a3, b3.value,
                                                          a4, b4.value));
    }

    /** {@inheritDoc} */
    @Override
    public Binary64 getPi() {
        return PI;
    }

}
