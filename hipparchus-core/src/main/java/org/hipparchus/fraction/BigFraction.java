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
package org.hipparchus.fraction;

import java.io.Serializable;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.function.Function;
import java.util.stream.Stream;

import org.hipparchus.FieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.fraction.ConvergentsIterator.ConvergenceStep;
import org.hipparchus.util.ArithmeticUtils;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Pair;
import org.hipparchus.util.Precision;

/**
 * Representation of a rational number without any overflow. This class is
 * immutable.
 *
 */
public class BigFraction
    extends Number
    implements FieldElement<BigFraction>, Comparable<BigFraction>, Serializable {

    /** A fraction representing "2 / 1". */
    public static final BigFraction TWO = new BigFraction(2);

    /** A fraction representing "1". */
    public static final BigFraction ONE = new BigFraction(1);

    /** A fraction representing "0". */
    public static final BigFraction ZERO = new BigFraction(0);

    /** A fraction representing "-1 / 1". */
    public static final BigFraction MINUS_ONE = new BigFraction(-1);

    /** A fraction representing "4/5". */
    public static final BigFraction FOUR_FIFTHS = new BigFraction(4, 5);

    /** A fraction representing "1/5". */
    public static final BigFraction ONE_FIFTH = new BigFraction(1, 5);

    /** A fraction representing "1/2". */
    public static final BigFraction ONE_HALF = new BigFraction(1, 2);

    /** A fraction representing "1/4". */
    public static final BigFraction ONE_QUARTER = new BigFraction(1, 4);

    /** A fraction representing "1/3". */
    public static final BigFraction ONE_THIRD = new BigFraction(1, 3);

    /** A fraction representing "3/5". */
    public static final BigFraction THREE_FIFTHS = new BigFraction(3, 5);

    /** A fraction representing "3/4". */
    public static final BigFraction THREE_QUARTERS = new BigFraction(3, 4);

    /** A fraction representing "2/5". */
    public static final BigFraction TWO_FIFTHS = new BigFraction(2, 5);

    /** A fraction representing "2/4". */
    public static final BigFraction TWO_QUARTERS = new BigFraction(2, 4);

    /** A fraction representing "2/3". */
    public static final BigFraction TWO_THIRDS = new BigFraction(2, 3);

    /** Serializable version identifier. */
    private static final long serialVersionUID = -5630213147331578515L;

    /** <code>BigInteger</code> representation of 100. */
    private static final BigInteger ONE_HUNDRED = BigInteger.valueOf(100);

    /** Convert a convergence step to the corresponding double fraction. */
    private static final Function<ConvergenceStep, BigFraction> STEP_TO_FRACTION = //
            s -> new BigFraction(s.getNumerator(), s.getDenominator());

    /** The numerator. */
    private final BigInteger numerator;

    /** The denominator. */
    private final BigInteger denominator;

    /**
     * <p>
     * Create a {@link BigFraction} equivalent to the passed {@code BigInteger}, ie
     * "num / 1".
     * </p>
     *
     * @param num
     *            the numerator.
     */
    public BigFraction(final BigInteger num) {
        this(num, BigInteger.ONE);
    }

    /**
     * Create a {@link BigFraction} given the numerator and denominator as
     * {@code BigInteger}. The {@link BigFraction} is reduced to lowest terms.
     *
     * @param num the numerator, must not be {@code null}.
     * @param den the denominator, must not be {@code null}.
     * @throws MathIllegalArgumentException if the denominator is zero.
     * @throws NullArgumentException if either of the arguments is null
     */
    public BigFraction(BigInteger num, BigInteger den) {
        MathUtils.checkNotNull(num, LocalizedCoreFormats.NUMERATOR);
        MathUtils.checkNotNull(den, LocalizedCoreFormats.DENOMINATOR);
        if (den.signum() == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_DENOMINATOR);
        }
        if (num.signum() == 0) {
            numerator   = BigInteger.ZERO;
            denominator = BigInteger.ONE;
        } else {

            // reduce numerator and denominator by greatest common denominator
            final BigInteger gcd = num.gcd(den);
            if (BigInteger.ONE.compareTo(gcd) < 0) {
                num = num.divide(gcd);
                den = den.divide(gcd);
            }

            // move sign to numerator
            if (den.signum() == -1) {
                num = num.negate();
                den = den.negate();
            }

            // store the values in the final fields
            numerator   = num;
            denominator = den;

        }
    }

    /**
     * Create a fraction given the double value.
     * <p>
     * This constructor behaves <em>differently</em> from
     * {@link #BigFraction(double, double, int)}. It converts the double value
     * exactly, considering its internal bits representation. This works for all
     * values except NaN and infinities and does not requires any loop or
     * convergence threshold.
     * </p>
     * <p>
     * Since this conversion is exact and since double numbers are sometimes
     * approximated, the fraction created may seem strange in some cases. For example,
     * calling <code>new BigFraction(1.0 / 3.0)</code> does <em>not</em> create
     * the fraction 1/3, but the fraction 6004799503160661 / 18014398509481984
     * because the double number passed to the constructor is not exactly 1/3
     * (this number cannot be stored exactly in IEEE754).
     * </p>
     * @see #BigFraction(double, double, int)
     * @param value the double value to convert to a fraction.
     * @exception MathIllegalArgumentException if value is NaN or infinite
     */
    public BigFraction(final double value) throws MathIllegalArgumentException {
        if (Double.isNaN(value)) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NAN_VALUE_CONVERSION);
        }
        if (Double.isInfinite(value)) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INFINITE_VALUE_CONVERSION);
        }

        // compute m and k such that value = m * 2^k
        final long bits     = Double.doubleToLongBits(value);
        final long sign     = bits & 0x8000000000000000L;
        final long exponent = bits & 0x7ff0000000000000L;
        long m              = bits & 0x000fffffffffffffL;
        if (exponent != 0) {
            // this was a normalized number, add the implicit most significant bit
            m |= 0x0010000000000000L;
        }
        if (sign != 0) {
            m = -m;
        }
        int k = ((int) (exponent >> 52)) - 1075;
        while (((m & 0x001ffffffffffffeL) != 0) && ((m & 0x1) == 0)) {
            m >>= 1;
            ++k;
        }

        if (k < 0) {
            numerator   = BigInteger.valueOf(m);
            denominator = BigInteger.ZERO.flipBit(-k);
        } else {
            numerator   = BigInteger.valueOf(m).multiply(BigInteger.ZERO.flipBit(k));
            denominator = BigInteger.ONE;
        }

    }

    /**
     * Create a fraction given the double value and maximum error allowed.
     * <p>* References:</p>
     * <ul>
     * <li><a href="http://mathworld.wolfram.com/ContinuedFraction.html">
     * Continued Fraction</a> equations (11) and (22)-(26)</li>
     * </ul>
     *
     * @param value
     *            the double value to convert to a fraction.
     * @param epsilon
     *            maximum error allowed. The resulting fraction is within
     *            <code>epsilon</code> of <code>value</code>, in absolute terms.
     * @param maxIterations
     *            maximum number of convergents.
     * @throws MathIllegalStateException
     *             if the continued fraction failed to converge.
     * @see #BigFraction(double)
     */
    public BigFraction(final double value, final double epsilon,
                       final int maxIterations)
        throws MathIllegalStateException {
        ConvergenceStep converged = ConvergentsIterator.convergent(value, maxIterations, s -> {
            final double quotient = s.getFractionValue();
            return Precision.equals(quotient, value, 1) || FastMath.abs(quotient - value) < epsilon;
        }).getKey();
        if (FastMath.abs(converged.getFractionValue() - value) < epsilon) {
            this.numerator = BigInteger.valueOf(converged.getNumerator());
            this.denominator = BigInteger.valueOf(converged.getDenominator());
        } else {
            throw new MathIllegalStateException(LocalizedCoreFormats.FAILED_FRACTION_CONVERSION,
                                                value, maxIterations);
        }
    }

    /**
     * Create a fraction given the double value and maximum denominator.
     * <p>* References:</p>
     * <ul>
     * <li><a href="http://mathworld.wolfram.com/ContinuedFraction.html">
     * Continued Fraction</a> equations (11) and (22)-(26)</li>
     * </ul>
     *
     * @param value
     *            the double value to convert to a fraction.
     * @param maxDenominator
     *            The maximum allowed value for denominator.
     * @throws MathIllegalStateException
     *             if the continued fraction failed to converge.
     */
    public BigFraction(final double value, final long maxDenominator)
        throws MathIllegalStateException {
        final int maxIterations = 100;
        ConvergenceStep[] lastValid = new ConvergenceStep[1];
        ConvergentsIterator.convergent(value, maxIterations, s -> {
            if (s.getDenominator() < maxDenominator) {
                lastValid[0] = s;
            }
            return Precision.equals(s.getFractionValue(), value, 1);
        });
        if (lastValid[0] != null) {
            this.numerator   = BigInteger.valueOf(lastValid[0].getNumerator());
            this.denominator = BigInteger.valueOf(lastValid[0].getDenominator());
        } else {
            throw new MathIllegalStateException(LocalizedCoreFormats.FAILED_FRACTION_CONVERSION,
                                                value, maxIterations);
        }
    }

    /**
     * <p>
     * Create a {@link BigFraction} equivalent to the passed {@code int}, ie
     * "num / 1".
     * </p>
     *
     * @param num
     *            the numerator.
     */
    public BigFraction(final int num) {
        this(BigInteger.valueOf(num), BigInteger.ONE);
    }

    /**
     * <p>
     * Create a {@link BigFraction} given the numerator and denominator as simple
     * {@code int}. The {@link BigFraction} is reduced to lowest terms.
     * </p>
     *
     * @param num
     *            the numerator.
     * @param den
     *            the denominator.
     */
    public BigFraction(final int num, final int den) {
        this(BigInteger.valueOf(num), BigInteger.valueOf(den));
    }

    /**
     * <p>
     * Create a {@link BigFraction} equivalent to the passed long, ie "num / 1".
     * </p>
     *
     * @param num
     *            the numerator.
     */
    public BigFraction(final long num) {
        this(BigInteger.valueOf(num), BigInteger.ONE);
    }

    /**
     * <p>
     * Create a {@link BigFraction} given the numerator and denominator as simple
     * {@code long}. The {@link BigFraction} is reduced to lowest terms.
     * </p>
     *
     * @param num
     *            the numerator.
     * @param den
     *            the denominator.
     */
    public BigFraction(final long num, final long den) {
        this(BigInteger.valueOf(num), BigInteger.valueOf(den));
    }

    /**
     * A test to determine if a series of fractions has converged.
     */
    @FunctionalInterface
    public interface ConvergenceTest {
        /**
         * Evaluates if the fraction formed by {@code numerator/denominator} satisfies
         * this convergence test.
         *
         * @param numerator   the numerator
         * @param denominator the denominator
         * @return if this convergence test is satisfied
         */
        boolean test(long numerator, long denominator); // NOPMD - this is not a Junit test, PMD false positive here
    }

    /** Generate a {@link Stream stream} of convergents from a real number.
     * @param value value to approximate
     * @param maxConvergents maximum number of convergents.
     * @return stream of {@link BigFraction} convergents approximating  {@code value}
     * @since 2.1
     */
    public static Stream<BigFraction> convergents(final double value, final int maxConvergents) {
        return ConvergentsIterator.convergents(value, maxConvergents).map(STEP_TO_FRACTION);
    }

    /**
     * Returns the last element of the series of convergent-steps to approximate the
     * given value.
     * <p>
     * The series terminates either at the first step that satisfies the given
     * {@code convergenceTest} or after at most {@code maxConvergents} elements. The
     * returned Pair consists of that terminal {@link BigFraction} and a
     * {@link Boolean} that indicates if it satisfies the given convergence tests.
     * If the returned pair's value is {@code false} the element at position
     * {@code maxConvergents} was examined but failed to satisfy the
     * {@code convergenceTest}. A caller can then decide to accept the result
     * nevertheless or to discard it. This method is usually faster than
     * {@link #convergents(double, int)} if only the terminal element is of
     * interest.
     *
     * @param value           value to approximate
     * @param maxConvergents  maximum number of convergents to examine
     * @param convergenceTest the test if the series has converged at a step
     * @return the pair of last element of the series of convergents and a boolean
     *         indicating if that element satisfies the specified convergent test
     */
    public static Pair<BigFraction, Boolean> convergent(double value, int maxConvergents,
            ConvergenceTest convergenceTest) {
        Pair<ConvergenceStep, Boolean> converged = ConvergentsIterator.convergent(value, maxConvergents,
                s -> convergenceTest.test(s.getNumerator(), s.getDenominator()));
        return Pair.create(STEP_TO_FRACTION.apply(converged.getKey()), converged.getValue());
    }

    /** {@inheritDoc} */
    @Override
    public double getReal() {
        return doubleValue();
    }

    /**
     * <p>
     * Creates a {@code BigFraction} instance with the 2 parts of a fraction
     * Y/Z.
     * </p>
     *
     * <p>
     * Any negative signs are resolved to be on the numerator.
     * </p>
     *
     * @param numerator
     *            the numerator, for example the three in 'three sevenths'.
     * @param denominator
     *            the denominator, for example the seven in 'three sevenths'.
     * @return a new fraction instance, with the numerator and denominator
     *         reduced.
     * @throws ArithmeticException
     *             if the denominator is <code>zero</code>.
     */
    public static BigFraction getReducedFraction(final int numerator,
                                                 final int denominator) {
        if (numerator == 0) {
            return ZERO; // normalize zero.
        }

        return new BigFraction(numerator, denominator);
    }

    /**
     * <p>
     * Returns the absolute value of this {@link BigFraction}.
     * </p>
     *
     * @return the absolute value as a {@link BigFraction}.
     */
    public BigFraction abs() {
        return (numerator.signum() == 1) ? this : negate();
    }

    /** Check if a fraction is an integer.
     * @return true of fraction is an integer
     */
    public boolean isInteger() {
        return denominator.equals(BigInteger.ONE);
    }

    /** Returns the signum function of this {@link BigFraction}.
     * <p>
     * The return value is -1 if the specified value is negative;
     * 0 if the specified value is zero; and 1 if the specified value is positive.
     * </p>
     * @return the signum function of this {@link BigFraction}
     * @since 1.7
     */
    public int signum() {
        return numerator.signum();
    }

    /**
     * <p>
     * Adds the value of this fraction to the passed {@link BigInteger},
     * returning the result in reduced form.
     * </p>
     *
     * @param bg
     *            the {@link BigInteger} to add, must'nt be <code>null</code>.
     * @return a {@code BigFraction} instance with the resulting values.
     * @throws NullArgumentException
     *             if the {@link BigInteger} is <code>null</code>.
     */
    public BigFraction add(final BigInteger bg) throws NullArgumentException {
        MathUtils.checkNotNull(bg);

        if (numerator.signum() == 0) {
            return new BigFraction(bg);
        }
        if (bg.signum() == 0) {
            return this;
        }

        return new BigFraction(numerator.add(denominator.multiply(bg)), denominator);
    }

    /**
     * <p>
     * Adds the value of this fraction to the passed {@code integer}, returning
     * the result in reduced form.
     * </p>
     *
     * @param i
     *            the {@code integer} to add.
     * @return a {@code BigFraction} instance with the resulting values.
     */
    public BigFraction add(final int i) {
        return add(BigInteger.valueOf(i));
    }

    /**
     * <p>
     * Adds the value of this fraction to the passed {@code long}, returning
     * the result in reduced form.
     * </p>
     *
     * @param l
     *            the {@code long} to add.
     * @return a {@code BigFraction} instance with the resulting values.
     */
    public BigFraction add(final long l) {
        return add(BigInteger.valueOf(l));
    }

    /**
     * <p>
     * Adds the value of this fraction to another, returning the result in
     * reduced form.
     * </p>
     *
     * @param fraction
     *            the {@link BigFraction} to add, must not be <code>null</code>.
     * @return a {@link BigFraction} instance with the resulting values.
     * @throws NullArgumentException if the {@link BigFraction} is {@code null}.
     */
    @Override
    public BigFraction add(final BigFraction fraction) {
        MathUtils.checkNotNull(fraction, LocalizedCoreFormats.FRACTION);
        if (fraction.numerator.signum() == 0) {
            return this;
        }
        if (numerator.signum() == 0) {
            return fraction;
        }

        BigInteger num;
        BigInteger den;
        if (denominator.equals(fraction.denominator)) {
            num = numerator.add(fraction.numerator);
            den = denominator;
        } else {
            num = (numerator.multiply(fraction.denominator)).add((fraction.numerator).multiply(denominator));
            den = denominator.multiply(fraction.denominator);
        }

        if (num.signum() == 0) {
            return ZERO;
        }

        return new BigFraction(num, den);

    }

    /**
     * <p>
     * Gets the fraction as a <code>BigDecimal</code>. This calculates the
     * fraction as the numerator divided by denominator.
     * </p>
     *
     * @return the fraction as a <code>BigDecimal</code>.
     * @throws ArithmeticException
     *             if the exact quotient does not have a terminating decimal
     *             expansion.
     * @see BigDecimal
     */
    public BigDecimal bigDecimalValue() {
        return new BigDecimal(numerator).divide(new BigDecimal(denominator));
    }

    /**
     * <p>
     * Gets the fraction as a <code>BigDecimal</code> following the passed
     * rounding mode. This calculates the fraction as the numerator divided by
     * denominator.
     * </p>
     *
     * @param roundingMode
     *            rounding mode to apply. see {@link BigDecimal} constants.
     * @return the fraction as a <code>BigDecimal</code>.
     * @throws IllegalArgumentException
     *             if {@code roundingMode} does not represent a valid rounding
     *             mode.
     * @see BigDecimal
     */
    public BigDecimal bigDecimalValue(final RoundingMode roundingMode) {
        return new BigDecimal(numerator).divide(new BigDecimal(denominator), roundingMode);
    }

    /**
     * <p>
     * Gets the fraction as a <code>BigDecimal</code> following the passed scale
     * and rounding mode. This calculates the fraction as the numerator divided
     * by denominator.
     * </p>
     *
     * @param scale
     *            scale of the <code>BigDecimal</code> quotient to be returned.
     *            see {@link BigDecimal} for more information.
     * @param roundingMode
     *            rounding mode to apply. see {@link BigDecimal} constants.
     * @return the fraction as a <code>BigDecimal</code>.
     * @see BigDecimal
     */
    public BigDecimal bigDecimalValue(final int scale, final RoundingMode roundingMode) {
        return new BigDecimal(numerator).divide(new BigDecimal(denominator), scale, roundingMode);
    }

    /**
     * <p>
     * Compares this object to another based on size.
     * </p>
     *
     * @param object
     *            the object to compare to, must not be <code>null</code>.
     * @return -1 if this is less than {@code object}, +1 if this is greater
     *         than {@code object}, 0 if they are equal.
     * @see java.lang.Comparable#compareTo(java.lang.Object)
     */
    @Override
    public int compareTo(final BigFraction object) {
        int lhsSigNum = numerator.signum();
        int rhsSigNum = object.numerator.signum();

        if (lhsSigNum != rhsSigNum) {
            return (lhsSigNum > rhsSigNum) ? 1 : -1;
        }
        if (lhsSigNum == 0) {
            return 0;
        }

        BigInteger nOd = numerator.multiply(object.denominator);
        BigInteger dOn = denominator.multiply(object.numerator);
        return nOd.compareTo(dOn);
    }

    /**
     * <p>
     * Divide the value of this fraction by the passed {@code BigInteger},
     * ie {@code this * 1 / bg}, returning the result in reduced form.
     * </p>
     *
     * @param bg the {@code BigInteger} to divide by, must not be {@code null}
     * @return a {@link BigFraction} instance with the resulting values
     * @throws NullArgumentException if the {@code BigInteger} is {@code null}
     * @throws MathRuntimeException if the fraction to divide by is zero
     */
    public BigFraction divide(final BigInteger bg) {
        MathUtils.checkNotNull(bg);
        if (bg.signum() == 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
        }
        if (numerator.signum() == 0) {
            return ZERO;
        }
        return new BigFraction(numerator, denominator.multiply(bg));
    }

    /**
     * <p>
     * Divide the value of this fraction by the passed {@code int}, ie
     * {@code this * 1 / i}, returning the result in reduced form.
     * </p>
     *
     * @param i the {@code int} to divide by
     * @return a {@link BigFraction} instance with the resulting values
     * @throws MathRuntimeException if the fraction to divide by is zero
     */
    public BigFraction divide(final int i) {
        return divide(BigInteger.valueOf(i));
    }

    /**
     * <p>
     * Divide the value of this fraction by the passed {@code long}, ie
     * {@code this * 1 / l}, returning the result in reduced form.
     * </p>
     *
     * @param l the {@code long} to divide by
     * @return a {@link BigFraction} instance with the resulting values
     * @throws MathRuntimeException if the fraction to divide by is zero
     */
    public BigFraction divide(final long l) {
        return divide(BigInteger.valueOf(l));
    }

    /**
     * <p>
     * Divide the value of this fraction by another, returning the result in
     * reduced form.
     * </p>
     *
     * @param fraction Fraction to divide by, must not be {@code null}.
     * @return a {@link BigFraction} instance with the resulting values.
     * @throws NullArgumentException if the {@code fraction} is {@code null}.
     * @throws MathRuntimeException if the fraction to divide by is zero
     */
    @Override
    public BigFraction divide(final BigFraction fraction) {
        MathUtils.checkNotNull(fraction, LocalizedCoreFormats.FRACTION);
        if (fraction.numerator.signum() == 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR);
        }
        if (numerator.signum() == 0) {
            return ZERO;
        }

        return multiply(fraction.reciprocal());
    }

    /**
     * <p>
     * Gets the fraction as a {@code double}. This calculates the fraction as
     * the numerator divided by denominator.
     * </p>
     *
     * @return the fraction as a {@code double}
     * @see java.lang.Number#doubleValue()
     */
    @Override
    public double doubleValue() {
        double result = numerator.doubleValue() / denominator.doubleValue();
        if (Double.isInfinite(result) || Double.isNaN(result)) {
            // Numerator and/or denominator must be out of range:
            // Calculate how far to shift them to put them in range.
            int shift = FastMath.max(numerator.bitLength(),
                                     denominator.bitLength()) - FastMath.getExponent(Double.MAX_VALUE);
            result = numerator.shiftRight(shift).doubleValue() /
                denominator.shiftRight(shift).doubleValue();
        }
        return result;
    }

    /**
     * <p>
     * Test for the equality of two fractions. If the lowest term numerator and
     * denominators are the same for both fractions, the two fractions are
     * considered to be equal.
     * </p>
     *
     * @param other
     *            fraction to test for equality to this fraction, can be
     *            <code>null</code>.
     * @return true if two fractions are equal, false if object is
     *         <code>null</code>, not an instance of {@link BigFraction}, or not
     *         equal to this fraction instance.
     * @see java.lang.Object#equals(java.lang.Object)
     */
    @Override
    public boolean equals(final Object other) {
        boolean ret = false;

        if (this == other) {
            ret = true;
        } else if (other instanceof BigFraction) {
            BigFraction rhs = (BigFraction) other;
            ret = numerator.equals(rhs.numerator) && denominator.equals(rhs.denominator);
        }

        return ret;
    }

    /**
     * <p>
     * Gets the fraction as a {@code float}. This calculates the fraction as
     * the numerator divided by denominator.
     * </p>
     *
     * @return the fraction as a {@code float}.
     * @see java.lang.Number#floatValue()
     */
    @Override
    public float floatValue() {
        float result = numerator.floatValue() / denominator.floatValue();
        if (Double.isNaN(result)) {
            // Numerator and/or denominator must be out of range:
            // Calculate how far to shift them to put them in range.
            int shift = FastMath.max(numerator.bitLength(),
                                     denominator.bitLength()) - FastMath.getExponent(Float.MAX_VALUE);
            result = numerator.shiftRight(shift).floatValue() /
                denominator.shiftRight(shift).floatValue();
        }
        return result;
    }

    /**
     * <p>
     * Access the denominator as a <code>BigInteger</code>.
     * </p>
     *
     * @return the denominator as a <code>BigInteger</code>.
     */
    public BigInteger getDenominator() {
        return denominator;
    }

    /**
     * <p>
     * Access the denominator as a {@code int}.
     * </p>
     *
     * @return the denominator as a {@code int}.
     */
    public int getDenominatorAsInt() {
        return denominator.intValue();
    }

    /**
     * <p>
     * Access the denominator as a {@code long}.
     * </p>
     *
     * @return the denominator as a {@code long}.
     */
    public long getDenominatorAsLong() {
        return denominator.longValue();
    }

    /**
     * <p>
     * Access the numerator as a <code>BigInteger</code>.
     * </p>
     *
     * @return the numerator as a <code>BigInteger</code>.
     */
    public BigInteger getNumerator() {
        return numerator;
    }

    /**
     * <p>
     * Access the numerator as a {@code int}.
     * </p>
     *
     * @return the numerator as a {@code int}.
     */
    public int getNumeratorAsInt() {
        return numerator.intValue();
    }

    /**
     * <p>
     * Access the numerator as a {@code long}.
     * </p>
     *
     * @return the numerator as a {@code long}.
     */
    public long getNumeratorAsLong() {
        return numerator.longValue();
    }

    /**
     * <p>
     * Gets a hashCode for the fraction.
     * </p>
     *
     * @return a hash code value for this object.
     * @see java.lang.Object#hashCode()
     */
    @Override
    public int hashCode() {
        return 37 * (37 * 17 + numerator.hashCode()) + denominator.hashCode();
    }

    /**
     * <p>
     * Gets the fraction as an {@code int}. This returns the whole number part
     * of the fraction.
     * </p>
     *
     * @return the whole number fraction part.
     * @see java.lang.Number#intValue()
     */
    @Override
    public int intValue() {
        return numerator.divide(denominator).intValue();
    }

    /**
     * <p>
     * Gets the fraction as a {@code long}. This returns the whole number part
     * of the fraction.
     * </p>
     *
     * @return the whole number fraction part.
     * @see java.lang.Number#longValue()
     */
    @Override
    public long longValue() {
        return numerator.divide(denominator).longValue();
    }

    /**
     * <p>
     * Multiplies the value of this fraction by the passed
     * <code>BigInteger</code>, returning the result in reduced form.
     * </p>
     *
     * @param bg the {@code BigInteger} to multiply by.
     * @return a {@code BigFraction} instance with the resulting values.
     * @throws NullArgumentException if {@code bg} is {@code null}.
     */
    public BigFraction multiply(final BigInteger bg) {
        MathUtils.checkNotNull(bg);
        if (numerator.signum() == 0 || bg.signum() == 0) {
            return ZERO;
        }
        return new BigFraction(bg.multiply(numerator), denominator);
    }

    /**
     * <p>
     * Multiply the value of this fraction by the passed {@code int}, returning
     * the result in reduced form.
     * </p>
     *
     * @param i
     *            the {@code int} to multiply by.
     * @return a {@link BigFraction} instance with the resulting values.
     */
    @Override
    public BigFraction multiply(final int i) {
        if (i == 0 || numerator.signum() == 0) {
            return ZERO;
        }

        return multiply(BigInteger.valueOf(i));
    }

    /**
     * <p>
     * Multiply the value of this fraction by the passed {@code long},
     * returning the result in reduced form.
     * </p>
     *
     * @param l
     *            the {@code long} to multiply by.
     * @return a {@link BigFraction} instance with the resulting values.
     */
    public BigFraction multiply(final long l) {
        if (l == 0 || numerator.signum() == 0) {
            return ZERO;
        }

        return multiply(BigInteger.valueOf(l));
    }

    /**
     * <p>
     * Multiplies the value of this fraction by another, returning the result in
     * reduced form.
     * </p>
     *
     * @param fraction Fraction to multiply by, must not be {@code null}.
     * @return a {@link BigFraction} instance with the resulting values.
     * @throws NullArgumentException if {@code fraction} is {@code null}.
     */
    @Override
    public BigFraction multiply(final BigFraction fraction) {
        MathUtils.checkNotNull(fraction, LocalizedCoreFormats.FRACTION);
        if (numerator.signum() == 0 ||
            fraction.numerator.signum() == 0) {
            return ZERO;
        }
        return new BigFraction(numerator.multiply(fraction.numerator),
                               denominator.multiply(fraction.denominator));
    }

    /**
     * <p>
     * Return the additive inverse of this fraction, returning the result in
     * reduced form.
     * </p>
     *
     * @return the negation of this fraction.
     */
    @Override
    public BigFraction negate() {
        return new BigFraction(numerator.negate(), denominator);
    }

    /**
     * <p>
     * Gets the fraction percentage as a {@code double}. This calculates the
     * fraction as the numerator divided by denominator multiplied by 100.
     * </p>
     *
     * @return the fraction percentage as a {@code double}.
     */
    public double percentageValue() {
        return multiply(ONE_HUNDRED).doubleValue();
    }

    /**
     * <p>
     * Returns a {@code BigFraction} whose value is
     * {@code (this<sup>exponent</sup>)}, returning the result in reduced form.
     * </p>
     *
     * @param exponent
     *            exponent to which this {@code BigFraction} is to be
     *            raised.
     * @return this<sup>exponent</sup>
     */
    public BigFraction pow(final int exponent) {
        if (exponent == 0) {
            return ONE;
        }
        if (numerator.signum() == 0) {
            return this;
        }

        if (exponent < 0) {
            return new BigFraction(denominator.pow(-exponent), numerator.pow(-exponent));
        }
        return new BigFraction(numerator.pow(exponent), denominator.pow(exponent));
    }

    /**
     * <p>
     * Returns a {@code BigFraction} whose value is
     * this<sup>exponent</sup>, returning the result in reduced form.
     * </p>
     *
     * @param exponent
     *            exponent to which this {@code BigFraction} is to be raised.
     * @return this<sup>exponent</sup> as a {@code BigFraction}.
     */
    public BigFraction pow(final long exponent) {
        if (exponent == 0) {
            return ONE;
        }
        if (numerator.signum() == 0) {
            return this;
        }

        if (exponent < 0) {
            return new BigFraction(ArithmeticUtils.pow(denominator, -exponent),
                                   ArithmeticUtils.pow(numerator,   -exponent));
        }
        return new BigFraction(ArithmeticUtils.pow(numerator,   exponent),
                               ArithmeticUtils.pow(denominator, exponent));
    }

    /**
     * <p>
     * Returns a {@code BigFraction} whose value is
     * this<sup>exponent</sup>, returning the result in reduced form.
     * </p>
     *
     * @param exponent
     *            exponent to which this {@code BigFraction} is to be raised.
     * @return this<sup>exponent</sup> as a {@code BigFraction}.
     */
    public BigFraction pow(final BigInteger exponent) {
        if (exponent.signum() == 0) {
            return ONE;
        }
        if (numerator.signum() == 0) {
            return this;
        }

        if (exponent.signum() == -1) {
            final BigInteger eNeg = exponent.negate();
            return new BigFraction(ArithmeticUtils.pow(denominator, eNeg),
                                   ArithmeticUtils.pow(numerator,   eNeg));
        }
        return new BigFraction(ArithmeticUtils.pow(numerator,   exponent),
                               ArithmeticUtils.pow(denominator, exponent));
    }

    /**
     * <p>
     * Returns a <code>double</code> whose value is
     * this<sup>exponent</sup>, returning the result in reduced form.
     * </p>
     *
     * @param exponent
     *            exponent to which this {@code BigFraction} is to be raised.
     * @return this<sup>exponent</sup>
     */
    public double pow(final double exponent) {
        return FastMath.pow(numerator.doubleValue(),   exponent) /
               FastMath.pow(denominator.doubleValue(), exponent);
    }

    /**
     * <p>
     * Return the multiplicative inverse of this fraction.
     * </p>
     *
     * @return the reciprocal fraction.
     */
    @Override
    public BigFraction reciprocal() {
        return new BigFraction(denominator, numerator);
    }

    /**
     * <p>
     * Reduce this {@code BigFraction} to its lowest terms.
     * </p>
     *
     * @return the reduced {@code BigFraction}. It doesn't change anything if
     *         the fraction can be reduced.
     */
    public BigFraction reduce() {
        final BigInteger gcd = numerator.gcd(denominator);

        if (BigInteger.ONE.compareTo(gcd) < 0) {
            return new BigFraction(numerator.divide(gcd), denominator.divide(gcd));
        } else {
            return this;
        }
    }

    /**
     * <p>
     * Subtracts the value of an {@link BigInteger} from the value of this
     * {@code BigFraction}, returning the result in reduced form.
     * </p>
     *
     * @param bg the {@link BigInteger} to subtract, cannot be {@code null}.
     * @return a {@code BigFraction} instance with the resulting values.
     * @throws NullArgumentException if the {@link BigInteger} is {@code null}.
     */
    public BigFraction subtract(final BigInteger bg) {
        MathUtils.checkNotNull(bg);
        if (bg.signum() == 0) {
            return this;
        }
        if (numerator.signum() == 0) {
            return new BigFraction(bg.negate());
        }

        return new BigFraction(numerator.subtract(denominator.multiply(bg)), denominator);
    }

    /**
     * <p>
     * Subtracts the value of an {@code integer} from the value of this
     * {@code BigFraction}, returning the result in reduced form.
     * </p>
     *
     * @param i the {@code integer} to subtract.
     * @return a {@code BigFraction} instance with the resulting values.
     */
    public BigFraction subtract(final int i) {
        return subtract(BigInteger.valueOf(i));
    }

    /**
     * <p>
     * Subtracts the value of a {@code long} from the value of this
     * {@code BigFraction}, returning the result in reduced form.
     * </p>
     *
     * @param l the {@code long} to subtract.
     * @return a {@code BigFraction} instance with the resulting values.
     */
    public BigFraction subtract(final long l) {
        return subtract(BigInteger.valueOf(l));
    }

    /**
     * <p>
     * Subtracts the value of another fraction from the value of this one,
     * returning the result in reduced form.
     * </p>
     *
     * @param fraction {@link BigFraction} to subtract, must not be {@code null}.
     * @return a {@link BigFraction} instance with the resulting values
     * @throws NullArgumentException if the {@code fraction} is {@code null}.
     */
    @Override
    public BigFraction subtract(final BigFraction fraction) {
        MathUtils.checkNotNull(fraction, LocalizedCoreFormats.FRACTION);
        if (fraction.numerator.signum() == 0) {
            return this;
        }
        if (numerator.signum() == 0) {
            return fraction.negate();
        }

        BigInteger num;
        BigInteger den;
        if (denominator.equals(fraction.denominator)) {
            num = numerator.subtract(fraction.numerator);
            den = denominator;
        } else {
            num = (numerator.multiply(fraction.denominator)).subtract((fraction.numerator).multiply(denominator));
            den = denominator.multiply(fraction.denominator);
        }
        return new BigFraction(num, den);

    }

    /**
     * <p>
     * Returns the <code>String</code> representing this fraction, ie
     * "num / dem" or just "num" if the denominator is one.
     * </p>
     *
     * @return a string representation of the fraction.
     * @see java.lang.Object#toString()
     */
    @Override
    public String toString() {
        if (BigInteger.ONE.equals(denominator)) {
            return numerator.toString();
        } else if (BigInteger.ZERO.equals(numerator)) {
            return "0";
        } else {
            return numerator + " / " + denominator;
        }
    }

    /** {@inheritDoc} */
    @Override
    public BigFractionField getField() {
        return BigFractionField.getInstance();
    }

}
