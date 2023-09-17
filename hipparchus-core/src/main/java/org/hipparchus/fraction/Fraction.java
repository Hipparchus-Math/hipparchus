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
import java.math.BigInteger;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Stream;

import org.hipparchus.FieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.fraction.ConvergentsIterator.ConvergenceStep;
import org.hipparchus.util.ArithmeticUtils;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Pair;
import org.hipparchus.util.Precision;

/**
 * Representation of a rational number.
 */
public class Fraction
    extends Number
    implements FieldElement<Fraction>, Comparable<Fraction>, Serializable {

    /** A fraction representing "2 / 1". */
    public static final Fraction TWO = new Fraction(2, 1);

    /** A fraction representing "1". */
    public static final Fraction ONE = new Fraction(1, 1);

    /** A fraction representing "0". */
    public static final Fraction ZERO = new Fraction(0, 1);

    /** A fraction representing "4/5". */
    public static final Fraction FOUR_FIFTHS = new Fraction(4, 5);

    /** A fraction representing "1/5". */
    public static final Fraction ONE_FIFTH = new Fraction(1, 5);

    /** A fraction representing "1/2". */
    public static final Fraction ONE_HALF = new Fraction(1, 2);

    /** A fraction representing "1/4". */
    public static final Fraction ONE_QUARTER = new Fraction(1, 4);

    /** A fraction representing "1/3". */
    public static final Fraction ONE_THIRD = new Fraction(1, 3);

    /** A fraction representing "3/5". */
    public static final Fraction THREE_FIFTHS = new Fraction(3, 5);

    /** A fraction representing "3/4". */
    public static final Fraction THREE_QUARTERS = new Fraction(3, 4);

    /** A fraction representing "2/5". */
    public static final Fraction TWO_FIFTHS = new Fraction(2, 5);

    /** A fraction representing "2/4". */
    public static final Fraction TWO_QUARTERS = new Fraction(2, 4);

    /** A fraction representing "2/3". */
    public static final Fraction TWO_THIRDS = new Fraction(2, 3);

    /** A fraction representing "-1 / 1". */
    public static final Fraction MINUS_ONE = new Fraction(-1, 1);

    /** Serializable version identifier */
    private static final long serialVersionUID = 3698073679419233275L;

    /** The default epsilon used for convergence. */
    private static final double DEFAULT_EPSILON = 1e-5;

    /** Convert a convergence step to the corresponding double fraction. */
    private static final Function<ConvergenceStep, Fraction> STEP_TO_FRACTION = //
                    s -> new Fraction((int) s.getNumerator(), (int) s.getDenominator());

    /** The denominator. */
    private final int denominator;

    /** The numerator. */
    private final int numerator;

    /**
     * Create a fraction given the double value.
     * @param value the double value to convert to a fraction.
     * @throws MathIllegalStateException if the continued fraction failed to
     *         converge.
     */
    public Fraction(double value) throws MathIllegalStateException {
        this(value, DEFAULT_EPSILON, 100);
    }

    /**
     * Create a fraction given the double value and maximum error allowed.
     * <p>
     * References:
     * <ul>
     * <li><a href="http://mathworld.wolfram.com/ContinuedFraction.html">
     * Continued Fraction</a> equations (11) and (22)-(26)</li>
     * </ul>
     *
     * @param value the double value to convert to a fraction.
     * @param epsilon maximum error allowed.  The resulting fraction is within
     *        {@code epsilon} of {@code value}, in absolute terms.
     * @param maxIterations maximum number of convergents
     * @throws MathIllegalStateException if the continued fraction failed to
     *         converge.
     */
    public Fraction(double value, double epsilon, int maxIterations)
        throws MathIllegalStateException {
        ConvergenceStep converged = convergent(value, maxIterations, s -> {
            double quotient = s.getFractionValue();
            return Precision.equals(quotient, value, 1) || FastMath.abs(quotient - value) < epsilon;
        }).getKey();
        if (FastMath.abs(converged.getFractionValue() - value) < epsilon) {
            this.numerator   = (int) converged.getNumerator();
            this.denominator = (int) converged.getDenominator();
        } else {
            throw new MathIllegalStateException(LocalizedCoreFormats.FAILED_FRACTION_CONVERSION,
                                                value, maxIterations);
        }
    }

    /**
     * Create a fraction given the double value and maximum denominator.
     * <p>
     * References:
     * <ul>
     * <li><a href="http://mathworld.wolfram.com/ContinuedFraction.html">
     * Continued Fraction</a> equations (11) and (22)-(26)</li>
     * </ul>
     *
     * @param value the double value to convert to a fraction.
     * @param maxDenominator The maximum allowed value for denominator
     * @throws MathIllegalStateException if the continued fraction failed to
     *         converge
     */
    public Fraction(double value, int maxDenominator)
        throws MathIllegalStateException {
        final int maxIterations = 100;
        ConvergenceStep[] lastValid = new ConvergenceStep[1];
        try {
            convergent(value, maxIterations, s -> {
                if (s.getDenominator() < maxDenominator) {
                    lastValid[0] = s;
                }
                return Precision.equals(s.getFractionValue(), value, 1);
            });
        } catch (MathIllegalStateException e) { // NOPMD - ignore overflows and just take the last valid result
        }
        if (lastValid[0] != null) {
            this.numerator   = (int) lastValid[0].getNumerator();
            this.denominator = (int) lastValid[0].getDenominator();
        } else {
            throw new MathIllegalStateException(LocalizedCoreFormats.FAILED_FRACTION_CONVERSION,
                                                value, maxIterations);
        }
    }

    /**
     * Create a fraction from an int.
     * The fraction is num / 1.
     * @param num the numerator.
     */
    public Fraction(int num) {
        this(num, 1);
    }

    /**
     * Create a fraction given the numerator and denominator.  The fraction is
     * reduced to lowest terms.
     * @param num the numerator.
     * @param den the denominator.
     * @throws MathRuntimeException if the denominator is {@code zero}
     */
    public Fraction(int num, int den) {
        if (den == 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR_IN_FRACTION,
                                           num, den);
        }
        if (den < 0) {
            if (num == Integer.MIN_VALUE ||
                den == Integer.MIN_VALUE) {
                throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_FRACTION,
                                               num, den);
            }
            num = -num;
            den = -den;
        }
        // reduce numerator and denominator by greatest common denominator.
        final int d = ArithmeticUtils.gcd(num, den);
        if (d > 1) {
            num /= d;
            den /= d;
        }

        // move sign to numerator.
        if (den < 0) {
            num = -num;
            den = -den;
        }
        this.numerator   = num;
        this.denominator = den;
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
        boolean test(int numerator, int denominator); // NOPMD - this is not a Junit test, PMD false positive here
    }

    /** Generate a {@link Stream stream} of convergents from a real number.
     * @param value value to approximate
     * @param maxConvergents maximum number of convergents.
     * @return stream of {@link Fraction} convergents approximating  {@code value}
     * @since 2.1
     */
    public static Stream<Fraction> convergents(final double value, final int maxConvergents) {
        if (FastMath.abs(value) > Integer.MAX_VALUE) {
            throw new MathIllegalStateException(LocalizedCoreFormats.FRACTION_CONVERSION_OVERFLOW, value, value, 1l);
        }
        return ConvergentsIterator.convergents(value, maxConvergents).map(STEP_TO_FRACTION);
    }

    /**
     * Returns the last element of the series of convergent-steps to approximate the
     * given value.
     * <p>
     * The series terminates either at the first step that satisfies the given
     * {@code convergenceTest} or after at most {@code maxConvergents} elements. The
     * returned Pair consists of that terminal {@link Fraction} and a
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
    public static Pair<Fraction, Boolean> convergent(double value, int maxConvergents, ConvergenceTest convergenceTest) {
        Pair<ConvergenceStep, Boolean> converged = convergent(value, maxConvergents, s -> {
            assertNoIntegerOverflow(s, value);
            return convergenceTest.test((int) s.getNumerator(), (int) s.getDenominator());
        });
        return Pair.create(STEP_TO_FRACTION.apply(converged.getKey()), converged.getValue());
    }

    /** Create a convergent-steps to approximate the given value.
     * @param value           value to approximate
     * @param maxConvergents  maximum number of convergents to examine
     * @param convergenceTests the test if the series has converged at a step
     * @return the pair of last element of the series of convergents and a boolean
     *         indicating if that element satisfies the specified convergent test
     */
    private static Pair<ConvergenceStep, Boolean> convergent(double value, int maxConvergents,
                                                             Predicate<ConvergenceStep> convergenceTests) {
        if (FastMath.abs(value) > Integer.MAX_VALUE) {
            throw new MathIllegalStateException(LocalizedCoreFormats.FRACTION_CONVERSION_OVERFLOW, value, value, 1l);
        }
        return ConvergentsIterator.convergent(value, maxConvergents, s -> {
            assertNoIntegerOverflow(s, value);
            return convergenceTests.test(s);
        });
    }

    /** Check no overflow occurred.
     * @param s convergent
     * @param value corresponding value
     */
    private static void assertNoIntegerOverflow(ConvergenceStep s, double value) {
        if (s.getNumerator() > Integer.MAX_VALUE || s.getDenominator() > Integer.MAX_VALUE) {
            throw new MathIllegalStateException(LocalizedCoreFormats.FRACTION_CONVERSION_OVERFLOW, value,
                    s.getNumerator(), s.getDenominator());
        }
    }

    /** {@inheritDoc} */
    @Override
    public double getReal() {
        return doubleValue();
    }

    /** Check if a fraction is an integer.
     * @return true of fraction is an integer
     */
    public boolean isInteger() {
        return denominator == 1;
    }

    /** Returns the signum function of this fraction.
     * <p>
     * The return value is -1 if the specified value is negative;
     * 0 if the specified value is zero; and 1 if the specified value is positive.
     * </p>
     * @return the signum function of this fraction
     * @since 1.7
     */
    public int signum() {
        return Integer.signum(numerator);
    }

    /**
     * Returns the absolute value of this fraction.
     * @return the absolute value.
     */
    public Fraction abs() {
        Fraction ret;
        if (numerator >= 0) {
            ret = this;
        } else {
            ret = negate();
        }
        return ret;
    }

    /**
     * Compares this object to another based on size.
     * @param object the object to compare to
     * @return -1 if this is less than {@code object}, +1 if this is greater
     *         than {@code object}, 0 if they are equal.
     */
    @Override
    public int compareTo(Fraction object) {
        long nOd = ((long) numerator) * object.denominator;
        long dOn = ((long) denominator) * object.numerator;
        return Long.compare(nOd, dOn);
    }

    /**
     * Gets the fraction as a {@code double}. This calculates the fraction as
     * the numerator divided by denominator.
     * @return the fraction as a {@code double}
     */
    @Override
    public double doubleValue() {
        return (double)numerator / (double)denominator;
    }

    /**
     * Test for the equality of two fractions.  If the lowest term
     * numerator and denominators are the same for both fractions, the two
     * fractions are considered to be equal.
     * @param other fraction to test for equality to this fraction
     * @return true if two fractions are equal, false if object is
     *         {@code null}, not an instance of {@link Fraction}, or not equal
     *         to this fraction instance.
     */
    @Override
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (other instanceof Fraction) {
            // since fractions are always in lowest terms, numerators and
            // denominators can be compared directly for equality.
            Fraction rhs = (Fraction)other;
            return (numerator == rhs.numerator) &&
                (denominator == rhs.denominator);
        }
        return false;
    }

    /**
     * Gets the fraction as a {@code float}. This calculates the fraction as
     * the numerator divided by denominator.
     * @return the fraction as a {@code float}
     */
    @Override
    public float floatValue() {
        return (float)doubleValue();
    }

    /**
     * Access the denominator.
     * @return the denominator.
     */
    public int getDenominator() {
        return denominator;
    }

    /**
     * Access the numerator.
     * @return the numerator.
     */
    public int getNumerator() {
        return numerator;
    }

    /**
     * Gets a hashCode for the fraction.
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        return 37 * (37 * 17 + numerator) + denominator;
    }

    /**
     * Gets the fraction as an {@code int}. This returns the whole number part
     * of the fraction.
     * @return the whole number fraction part
     */
    @Override
    public int intValue() {
        return (int)doubleValue();
    }

    /**
     * Gets the fraction as a {@code long}. This returns the whole number part
     * of the fraction.
     * @return the whole number fraction part
     */
    @Override
    public long longValue() {
        return (long)doubleValue();
    }

    /**
     * Return the additive inverse of this fraction.
     * @return the negation of this fraction.
     */
    @Override
    public Fraction negate() {
        if (numerator==Integer.MIN_VALUE) {
            throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_FRACTION, numerator, denominator);
        }
        return new Fraction(-numerator, denominator);
    }

    /**
     * Return the multiplicative inverse of this fraction.
     * @return the reciprocal fraction
     */
    @Override
    public Fraction reciprocal() {
        return new Fraction(denominator, numerator);
    }

    /**
     * Adds the value of this fraction to another, returning the result in reduced form.
     * The algorithm follows Knuth, 4.5.1.
     *
     * @param fraction  the fraction to add, must not be {@code null}
     * @return a {@code Fraction} instance with the resulting values
     * @throws org.hipparchus.exception.NullArgumentException if the fraction is {@code null}
     * @throws MathRuntimeException if the resulting numerator or denominator exceeds
     *  {@code Integer.MAX_VALUE}
     */
    @Override
    public Fraction add(Fraction fraction) {
        return addSub(fraction, true /* add */);
    }

    /**
     * Add an integer to the fraction.
     * @param i the {@code integer} to add.
     * @return this + i
     */
    public Fraction add(final int i) {
        return new Fraction(numerator + i * denominator, denominator);
    }

    /**
     * Subtracts the value of another fraction from the value of this one,
     * returning the result in reduced form.
     *
     * @param fraction  the fraction to subtract, must not be {@code null}
     * @return a {@code Fraction} instance with the resulting values
     * @throws org.hipparchus.exception.NullArgumentException if the fraction is {@code null}
     * @throws MathRuntimeException if the resulting numerator or denominator
     *   cannot be represented in an {@code int}.
     */
    @Override
    public Fraction subtract(Fraction fraction) {
        return addSub(fraction, false /* subtract */);
    }

    /**
     * Subtract an integer from the fraction.
     * @param i the {@code integer} to subtract.
     * @return this - i
     */
    public Fraction subtract(final int i) {
        return new Fraction(numerator - i * denominator, denominator);
    }

    /**
     * Implement add and subtract using algorithm described in Knuth 4.5.1.
     *
     * @param fraction the fraction to subtract, must not be {@code null}
     * @param isAdd true to add, false to subtract
     * @return a {@code Fraction} instance with the resulting values
     * @throws org.hipparchus.exception.NullArgumentException if the fraction is {@code null}
     * @throws MathRuntimeException if the resulting numerator or denominator
     *   cannot be represented in an {@code int}.
     */
    private Fraction addSub(Fraction fraction, boolean isAdd) {
        MathUtils.checkNotNull(fraction, LocalizedCoreFormats.FRACTION);

        // zero is identity for addition.
        if (numerator == 0) {
            return isAdd ? fraction : fraction.negate();
        }
        if (fraction.numerator == 0) {
            return this;
        }
        // if denominators are randomly distributed, d1 will be 1 about 61%
        // of the time.
        int d1 = ArithmeticUtils.gcd(denominator, fraction.denominator);
        if (d1==1) {
            // result is ( (u*v' +/- u'v) / u'v')
            int uvp = ArithmeticUtils.mulAndCheck(numerator, fraction.denominator);
            int upv = ArithmeticUtils.mulAndCheck(fraction.numerator, denominator);
            return new Fraction
                (isAdd ? ArithmeticUtils.addAndCheck(uvp, upv) :
                 ArithmeticUtils.subAndCheck(uvp, upv),
                 ArithmeticUtils.mulAndCheck(denominator, fraction.denominator));
        }
        // the quantity 't' requires 65 bits of precision; see knuth 4.5.1
        // exercise 7.  we're going to use a BigInteger.
        // t = u(v'/d1) +/- v(u'/d1)
        BigInteger uvp = BigInteger.valueOf(numerator)
                                   .multiply(BigInteger.valueOf(fraction.denominator / d1));
        BigInteger upv = BigInteger.valueOf(fraction.numerator)
                                   .multiply(BigInteger.valueOf(denominator / d1));
        BigInteger t = isAdd ? uvp.add(upv) : uvp.subtract(upv);
        // but d2 doesn't need extra precision because
        // d2 = gcd(t,d1) = gcd(t mod d1, d1)
        int tmodd1 = t.mod(BigInteger.valueOf(d1)).intValue();
        int d2 = (tmodd1==0)?d1:ArithmeticUtils.gcd(tmodd1, d1);

        // result is (t/d2) / (u'/d1)(v'/d2)
        BigInteger w = t.divide(BigInteger.valueOf(d2));
        if (w.bitLength() > 31) {
            throw new MathRuntimeException(LocalizedCoreFormats.NUMERATOR_OVERFLOW_AFTER_MULTIPLY,
                                           w);
        }
        return new Fraction (w.intValue(),
                ArithmeticUtils.mulAndCheck(denominator/d1,
                                            fraction.denominator/d2));
    }

    /**
     * Multiplies the value of this fraction by another, returning the
     * result in reduced form.
     *
     * @param fraction  the fraction to multiply by, must not be {@code null}
     * @return a {@code Fraction} instance with the resulting values
     * @throws org.hipparchus.exception.NullArgumentException if the fraction is {@code null}
     * @throws MathRuntimeException if the resulting numerator or denominator exceeds
     *  {@code Integer.MAX_VALUE}
     */
    @Override
    public Fraction multiply(Fraction fraction) {
        MathUtils.checkNotNull(fraction, LocalizedCoreFormats.FRACTION);
        if (numerator == 0 || fraction.numerator == 0) {
            return ZERO;
        }
        // knuth 4.5.1
        // make sure we don't overflow unless the result *must* overflow.
        int d1 = ArithmeticUtils.gcd(numerator, fraction.denominator);
        int d2 = ArithmeticUtils.gcd(fraction.numerator, denominator);
        return getReducedFraction
                (ArithmeticUtils.mulAndCheck(numerator/d1, fraction.numerator/d2),
                 ArithmeticUtils.mulAndCheck(denominator/d2, fraction.denominator/d1));
    }

    /**
     * Multiply the fraction by an integer.
     * @param i the {@code integer} to multiply by.
     * @return this * i
     */
    @Override
    public Fraction multiply(final int i) {
        return multiply(new Fraction(i));
    }

    /**
     * Divide the value of this fraction by another.
     *
     * @param fraction  the fraction to divide by, must not be {@code null}
     * @return a {@code Fraction} instance with the resulting values
     * @throws IllegalArgumentException if the fraction is {@code null}
     * @throws MathRuntimeException if the fraction to divide by is zero
     * @throws MathRuntimeException if the resulting numerator or denominator exceeds
     *  {@code Integer.MAX_VALUE}
     */
    @Override
    public Fraction divide(Fraction fraction) {
        MathUtils.checkNotNull(fraction, LocalizedCoreFormats.FRACTION);
        if (fraction.numerator == 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_FRACTION_TO_DIVIDE_BY,
                                           fraction.numerator, fraction.denominator);
        }
        return multiply(fraction.reciprocal());
    }

    /**
     * Divide the fraction by an integer.
     * @param i the {@code integer} to divide by.
     * @return this * i
     */
    public Fraction divide(final int i) {
        return divide(new Fraction(i));
    }

    /**
     * Gets the fraction percentage as a {@code double}. This calculates the
     * fraction as the numerator divided by denominator multiplied by 100.
     *
     * @return the fraction percentage as a {@code double}.
     */
    public double percentageValue() {
        return 100 * doubleValue();
    }

    /**
     * Creates a {@code Fraction} instance with the 2 parts
     * of a fraction Y/Z.
     * <p>
     * Any negative signs are resolved to be on the numerator.
     *
     * @param numerator  the numerator, for example the three in 'three sevenths'
     * @param denominator  the denominator, for example the seven in 'three sevenths'
     * @return a new fraction instance, with the numerator and denominator reduced
     * @throws MathRuntimeException if the denominator is {@code zero}
     */
    public static Fraction getReducedFraction(int numerator, int denominator) {
        if (denominator == 0) {
            throw new MathRuntimeException(LocalizedCoreFormats.ZERO_DENOMINATOR_IN_FRACTION,
                                           numerator, denominator);
        }
        if (numerator==0) {
            return ZERO; // normalize zero.
        }
        // allow 2^k/-2^31 as a valid fraction (where k>0)
        if (denominator==Integer.MIN_VALUE && (numerator&1)==0) {
            numerator/=2; denominator/=2;
        }
        if (denominator < 0) {
            if (numerator==Integer.MIN_VALUE ||
                denominator==Integer.MIN_VALUE) {
                throw new MathRuntimeException(LocalizedCoreFormats.OVERFLOW_IN_FRACTION,
                                               numerator, denominator);
            }
            numerator = -numerator;
            denominator = -denominator;
        }
        // simplify fraction.
        int gcd = ArithmeticUtils.gcd(numerator, denominator);
        numerator /= gcd;
        denominator /= gcd;
        return new Fraction(numerator, denominator);
    }

    /**
     * Returns the {@code String} representing this fraction, ie
     * "num / dem" or just "num" if the denominator is one.
     *
     * @return a string representation of the fraction.
     * @see java.lang.Object#toString()
     */
    @Override
    public String toString() {
        if (denominator == 1) {
            return Integer.toString(numerator);
        } else if (numerator == 0) {
            return "0";
        } else {
            return numerator + " / " + denominator;
        }
    }

    /** {@inheritDoc} */
    @Override
    public FractionField getField() {
        return FractionField.getInstance();
    }

}
