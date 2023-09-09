/* Copyright 2018 Ulf Adams
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

/*
 * This is not the original file distributed by Ulf Adams in project
 * https://github.com/ulfjack/ryu
 * It has been modified by the Hipparchus project.
 */
package org.hipparchus.util;

import java.math.BigInteger;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * An implementation of Ryū for double.
 * <p>
 * Ryū generates the shortest decimal representation of a floating point number
 * that maintains round-trip safety. That is, a correct parser can recover the
 * exact original number. Ryū is very fast (about 10 time faster than {@code
 * Double.toString()}).
 * </p>
 * @see <a href="https://dl.acm.org/citation.cfm?doid=3296979.3192369">Ryū: fast float-to-string conversion</a>
 */
public final class RyuDouble {

    /** Default low switch level to scientific notation. */
    public static final int DEFAULT_LOW_EXP = -3;

    /** Default high switch level to scientific notation. */
    public static final int DEFAULT_HIGH_EXP = 7;

    /** Number of bits in a double mantissa. */
    private static final int DOUBLE_MANTISSA_BITS = 52;

    /** Bit mask for retrieving mantissa. */
    private static final long DOUBLE_MANTISSA_MASK = (1L << DOUBLE_MANTISSA_BITS) - 1;

    /** Number of bits in a double exponant. */
    private static final int DOUBLE_EXPONENT_BITS = 11;

    /** Bit mask for retrieving exponent. */
    private static final int DOUBLE_EXPONENT_MASK = (1 << DOUBLE_EXPONENT_BITS) - 1;

    /** Bias of the exponent. */
    private static final int DOUBLE_EXPONENT_BIAS = (1 << (DOUBLE_EXPONENT_BITS - 1)) - 1;

    /** Size of the factors table for positive exponents. */
    private static final int POS_TABLE_SIZE = 326;

    /** Size of the factors table for negative exponents. */
    private static final int NEG_TABLE_SIZE = 291;

    /** Bit count for complete entries in the positive exponent tables. */
    private static final int POW5_BITCOUNT = 121; // max 3*31 = 124

    /** Bit count for split entries in the positive exponent tables. */
    private static final int POW5_QUARTER_BITCOUNT = 31;

    /** Split table for positive exponents. */
    private static final int[][] POW5_SPLIT = new int[POS_TABLE_SIZE][4];

    /** Bit count for complete entries in the negative exponent tables. */
    private static final int POW5_INV_BITCOUNT = 122; // max 3*31 = 124

    /** Bit count for split entries in the negative exponent tables. */
    private static final int POW5_INV_QUARTER_BITCOUNT = 31;

    /** Split table for negative exponents. */
    private static final int[][] POW5_INV_SPLIT = new int[NEG_TABLE_SIZE][4];

    /** Create the tables. */
    static {
        final BigInteger mask = BigInteger.valueOf(1).shiftLeft(POW5_QUARTER_BITCOUNT).subtract(BigInteger.ONE);
        final BigInteger invMask = BigInteger.valueOf(1).shiftLeft(POW5_INV_QUARTER_BITCOUNT).subtract(BigInteger.ONE);
        for (int i = 0; i < FastMath.max(POS_TABLE_SIZE, NEG_TABLE_SIZE); i++) {
            final BigInteger pow = BigInteger.valueOf(5).pow(i);
            final int pow5len = pow.bitLength();
            if (i < POW5_SPLIT.length) {
                for (int j = 0; j < 4; j++) {
                    POW5_SPLIT[i][j] = pow.
                                       shiftRight(pow5len - POW5_BITCOUNT + (3 - j) * POW5_QUARTER_BITCOUNT).
                                       and(mask).
                                       intValueExact();
                }
            }

            if (i < POW5_INV_SPLIT.length) {
                // We want floor(log_2 5^q) here, which is pow5len - 1.
                final int j = pow5len - 1 + POW5_INV_BITCOUNT;
                final BigInteger inv = BigInteger.ONE.shiftLeft(j).divide(pow).add(BigInteger.ONE);
                for (int k = 0; k < 4; k++) {
                    if (k == 0) {
                        POW5_INV_SPLIT[i][k] = inv.shiftRight((3 - k) * POW5_INV_QUARTER_BITCOUNT).intValueExact();
                    } else {
                        POW5_INV_SPLIT[i][k] = inv.shiftRight((3 - k) * POW5_INV_QUARTER_BITCOUNT).and(invMask).intValueExact();
                    }
                }
            }
        }
    }

    /** Private constructor for a utility class.
     */
    private RyuDouble() {
        // nothing to do
    }

    /** Convert a double to shortest string representation, preserving full accuracy.
     * <p>
     * This implementation uses the same specifications as {@code Double.toString()},
     * i.e. it uses scientific notation if for numbers smaller than 10⁻³ or larger
     * than 10⁺⁷, and decimal notion in between. That is it call {@link #doubleToString(double,
     * int, int) doubleToString(value, -3, 7)}.
     * </p>
     * @param value double number to convert
     * @return shortest string representation
     * @see #doubleToString(double, int, int)
     * @see #DEFAULT_LOW_EXP
     * @see #DEFAULT_HIGH_EXP
     */
    public static String doubleToString(double value) {
        return doubleToString(value, DEFAULT_LOW_EXP, DEFAULT_HIGH_EXP);
    }

    /** Convert a double to shortest string representation, preserving full accuracy.
     * <p>
     * Number inside of the interval [10<sup>lowExp</sup>, 10<sup>highExp</sup>]
     * are represented using decimal notation, numbers outside of this
     * range are represented using scientific notation.
     * </p>
     * @param value double number to convert
     * @param lowExp lowest decimal exponent for which decimal notation can be used
     * @param highExp highest decimal exponent for which decimal notation can be used
     * @return shortest string representation
     * @see #doubleToString(double)
     * @see #DEFAULT_LOW_EXP
     * @see #DEFAULT_HIGH_EXP
     */
    public static String doubleToString(double value, int lowExp, int highExp) {
        // Step 1: Decode the floating point number, and unify normalized and subnormal cases.
        // First, handle all the trivial cases.
        if (Double.isNaN(value)) {
            return "NaN";
        }
        if (value == Double.POSITIVE_INFINITY) {
            return "Infinity";
        }
        if (value == Double.NEGATIVE_INFINITY) {
            return "-Infinity";
        }
        long bits = Double.doubleToLongBits(value);
        if (bits == 0) {
            return "0.0";
        }
        if (bits == 0x8000000000000000L) {
            return "-0.0";
        }

        // Otherwise extract the mantissa and exponent bits and run the full algorithm.
        final int ieeeExponent = (int) ((bits >>> DOUBLE_MANTISSA_BITS) & DOUBLE_EXPONENT_MASK);
        final long ieeeMantissa = bits & DOUBLE_MANTISSA_MASK;
        int e2;
        final long m2;
        if (ieeeExponent == 0) {
            // Denormal number - no implicit leading 1, and the exponent is 1, not 0.
            e2 = 1 - DOUBLE_EXPONENT_BIAS - DOUBLE_MANTISSA_BITS;
            m2 = ieeeMantissa;
        } else {
            // Add implicit leading 1.
            e2 = ieeeExponent - DOUBLE_EXPONENT_BIAS - DOUBLE_MANTISSA_BITS;
            m2 = ieeeMantissa | (1L << DOUBLE_MANTISSA_BITS);
        }

        final boolean sign = bits < 0;

        // Step 2: Determine the interval of legal decimal representations.
        final boolean even = (m2 & 1) == 0;
        final long mv = 4 * m2;
        final long mp = 4 * m2 + 2;
        final int mmShift = ((m2 != (1L << DOUBLE_MANTISSA_BITS)) || (ieeeExponent <= 1)) ? 1 : 0;
        final long mm = 4 * m2 - 1 - mmShift;
        e2 -= 2;

        // Step 3: Convert to a decimal power base using 128-bit arithmetic.
        // -1077 = 1 - 1023 - 53 - 2 <= e_2 - 2 <= 2046 - 1023 - 53 - 2 = 968
        long dv;
        long dp;
        long dm;
        final int e10;
        boolean dmIsTrailingZeros = false;
        boolean dvIsTrailingZeros = false;
        if (e2 >= 0) {
            final int q = FastMath.max(0, ((e2 * 78913) >>> 18) - 1);
            // k = constant + floor(log_2(5^q))
            final int k = POW5_INV_BITCOUNT + pow5bits(q) - 1;
            final int i = -e2 + q + k;
            dv = mulPow5InvDivPow2(mv, q, i);
            dp = mulPow5InvDivPow2(mp, q, i);
            dm = mulPow5InvDivPow2(mm, q, i);
            e10 = q;

            if (q <= 21) {
                if (mv % 5 == 0) {
                    dvIsTrailingZeros = multipleOfPowerOf5(mv, q);
                } else if (even) {
                    dmIsTrailingZeros = multipleOfPowerOf5(mm, q);
                } else if (multipleOfPowerOf5(mp, q)) {
                    dp--;
                }
            }
        } else {
            final int q = FastMath.max(0, ((-e2 * 732923) >>> 20) - 1);
            final int i = -e2 - q;
            final int k = pow5bits(i) - POW5_BITCOUNT;
            final int j = q - k;
            dv = mulPow5divPow2(mv, i, j);
            dp = mulPow5divPow2(mp, i, j);
            dm = mulPow5divPow2(mm, i, j);
            e10 = q + e2;
            if (q <= 1) {
                dvIsTrailingZeros = true;
                if (even) {
                    dmIsTrailingZeros = mmShift == 1;
                } else {
                    dp--;
                }
            } else if (q < 63) {
                dvIsTrailingZeros = (mv & ((1L << (q - 1)) - 1)) == 0;
            }
        }

        // Step 4: Find the shortest decimal representation in the interval of legal representations.
        //
        // We do some extra work here in order to follow Float/Double.toString semantics. In particular,
        // that requires printing in scientific format if and only if the exponent is between lowExp and highExp,
        // and it requires printing at least two decimal digits.
        //
        // Above, we moved the decimal dot all the way to the right, so now we need to count digits to
        // figure out the correct exponent for scientific notation.
        final int vplength = decimalLength(dp);
        int exp = e10 + vplength - 1;

        // use scientific notation if and only if outside this range.
        final boolean scientificNotation = !((exp >= lowExp) && (exp < highExp));

        int removed = 0;

        int lastRemovedDigit = 0;
        long output;
        if (dmIsTrailingZeros || dvIsTrailingZeros) {
            while (dp / 10 > dm / 10) {
                if ((dp < 100) && scientificNotation) {
                    // Double.toString semantics requires printing at least two digits.
                    break;
                }
                dmIsTrailingZeros &= dm % 10 == 0;
                dvIsTrailingZeros &= lastRemovedDigit == 0;
                lastRemovedDigit = (int) (dv % 10);
                dp /= 10;
                dv /= 10;
                dm /= 10;
                removed++;
            }
            if (dmIsTrailingZeros && even) {
                while (dm % 10 == 0) {
                    if ((dp < 100) && scientificNotation) {
                        // Double.toString semantics requires printing at least two digits.
                        break;
                    }
                    dvIsTrailingZeros &= lastRemovedDigit == 0;
                    lastRemovedDigit = (int) (dv % 10);
                    dp /= 10;
                    dv /= 10;
                    dm /= 10;
                    removed++;
                }
            }
            if (dvIsTrailingZeros && (lastRemovedDigit == 5) && (dv % 2 == 0)) {
                // Round even if the exact numbers is .....50..0.
                lastRemovedDigit = 4;
            }
            output = dv +
                            ((dv == dm && !(dmIsTrailingZeros && even)) || (lastRemovedDigit >= 5) ? 1 : 0);
        } else {
            while (dp / 10 > dm / 10) {
                if ((dp < 100) && scientificNotation) {
                    // Double.toString semantics requires printing at least two digits.
                    break;
                }
                lastRemovedDigit = (int) (dv % 10);
                dp /= 10;
                dv /= 10;
                dm /= 10;
                removed++;
            }
            output = dv + ((dv == dm || (lastRemovedDigit >= 5)) ? 1 : 0);
        }
        int olength = vplength - removed;

        // Step 5: Print the decimal representation.
        // We follow Double.toString semantics here,
        // but adjusting the boundaries at which we switch to scientific notation
        char[] result = new char[14 - lowExp + highExp];
        int index = 0;
        if (sign) {
            result[index++] = '-';
        }

        // Values in the interval [10^lowExp, 10^highExp) are special.
        if (scientificNotation) {
            // Print in the format x.xxxxxE-yy.
            for (int i = 0; i < olength - 1; i++) {
                int c = (int) (output % 10);
                output /= 10;
                result[index + olength - i] = (char) ('0' + c);
            }
            result[index] = (char) ('0' + output % 10);
            result[index + 1] = '.';
            index += olength + 1;
            if (olength == 1) {
                result[index++] = '0';
            }

            // Print 'E', the exponent sign, and the exponent, which has at most three digits.
            result[index++] = 'E';
            if (exp < 0) {
                result[index++] = '-';
                exp = -exp;
            }
            if (exp >= 100) {
                result[index++] = (char) ('0' + exp / 100);
                exp %= 100;
                result[index++] = (char) ('0' + exp / 10);
            } else if (exp >= 10) {
                result[index++] = (char) ('0' + exp / 10);
            }
            result[index++] = (char) ('0' + exp % 10);
            return String.valueOf(result, 0, index);
        } else {
            // Otherwise follow the Java spec for values in the interval [10^lowExp, 10^highExp).
            if (exp < 0) {
                // Decimal dot is before any of the digits.
                result[index++] = '0';
                result[index++] = '.';
                for (int i = -1; i > exp; i--) {
                    result[index++] = '0';
                }
                int current = index;
                for (int i = 0; i < olength; i++) {
                    result[current + olength - i - 1] = (char) ('0' + output % 10);
                    output /= 10;
                    index++;
                }
            } else if (exp + 1 >= olength) {
                // Decimal dot is after any of the digits.
                for (int i = 0; i < olength; i++) {
                    result[index + olength - i - 1] = (char) ('0' + output % 10);
                    output /= 10;
                }
                index += olength;
                for (int i = olength; i < exp + 1; i++) {
                    result[index++] = '0';
                }
                result[index++] = '.';
                result[index++] = '0';
            } else {
                // Decimal dot is somewhere between the digits.
                int current = index + 1;
                for (int i = 0; i < olength; i++) {
                    if (olength - i - 1 == exp) {
                        result[current + olength - i - 1] = '.';
                        current--;
                    }
                    result[current + olength - i - 1] = (char) ('0' + output % 10);
                    output /= 10;
                }
                index += olength + 1;
            }
            return String.valueOf(result, 0, index);
        }
    }

    /** Get the number of bits of 5<sup>e</sup>.
     * @param e exponent
     * @return number of bits of 5<sup>e</sup>
     */
    private static int pow5bits(int e) {
        return ((e * 1217359) >>> 19) + 1;
    }

    /** Compute decimal length of an integer.
     * @param v integer to check
     * @return decimal length of {@code v}
     */
    private static int decimalLength(long v) {
        if (v >= 1000000000000000000L) {
            return 19;
        }
        if (v >= 100000000000000000L) {
            return 18;
        }
        if (v >= 10000000000000000L) {
            return 17;
        }
        if (v >= 1000000000000000L) {
            return 16;
        }
        if (v >= 100000000000000L) {
            return 15;
        }
        if (v >= 10000000000000L) {
            return 14;
        }
        if (v >= 1000000000000L) {
            return 13;
        }
        if (v >= 100000000000L) {
            return 12;
        }
        if (v >= 10000000000L) {
            return 11;
        }
        if (v >= 1000000000L) {
            return 10;
        }
        if (v >= 100000000L) {
            return 9;
        }
        if (v >= 10000000L) {
            return 8;
        }
        if (v >= 1000000L) {
            return 7;
        }
        if (v >= 100000L) {
            return 6;
        }
        if (v >= 10000L) {
            return 5;
        }
        if (v >= 1000L) {
            return 4;
        }
        if (v >= 100L) {
            return 3;
        }
        if (v >= 10L) {
            return 2;
        }
        return 1;
    }

    private static boolean multipleOfPowerOf5(long value, int q) {
        return pow5Factor(value) >= q;
    }

    /** Compute largest power of 5 that divides the value.
     * @param value value to check
     * @return largest power of 5 that divides the value
     */
    private static int pow5Factor(long value) {
        // We want to find the largest power of 5 that divides value.
        if ((value % 5) != 0) {
            return 0;
        }
        if ((value % 25) != 0) {
            return 1;
        }
        if ((value % 125) != 0) {
            return 2;
        }
        if ((value % 625) != 0) {
            return 3;
        }
        int count = 4;
        value /= 625;
        while (value > 0) {
            if (value % 5 != 0) {
                return count;
            }
            value /= 5;
            count++;
        }
        throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL, value, 0);
    }

    /**
     * Compute the high digits of m * 5^p / 10^q = m * 5^(p - q) / 2^q = m * 5^i / 2^j, with q chosen
     * such that m * 5^i / 2^j has sufficiently many decimal digits to represent the original floating
     * point number.
     * @param m mantissa
     * @param i power of 5
     * @param j power of 2
     * @return high digits of m * 5^i / 2^j
     */
    private static long mulPow5divPow2(final long m, final int i, final int j) {
        // m has at most 55 bits.
        long mHigh = m >>> 31;
        long mLow = m & 0x7fffffff;
        long bits13 = mHigh * POW5_SPLIT[i][0]; // 124
        long bits03 = mLow * POW5_SPLIT[i][0];  // 93
        long bits12 = mHigh * POW5_SPLIT[i][1]; // 93
        long bits02 = mLow * POW5_SPLIT[i][1];  // 62
        long bits11 = mHigh * POW5_SPLIT[i][2]; // 62
        long bits01 = mLow * POW5_SPLIT[i][2];  // 31
        long bits10 = mHigh * POW5_SPLIT[i][3]; // 31
        long bits00 = mLow * POW5_SPLIT[i][3];  // 0
        int actualShift = j - 3 * 31 - 21;
        if (actualShift < 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL,
                                                   j, 3 * 31 + 21);
        }
        return ((((((((bits00 >>> 31) + bits01 + bits10) >>> 31) +
                        bits02 + bits11) >>> 31) +
                        bits03 + bits12) >>> 21) +
                        (bits13 << 10)) >>> actualShift;
    }

    /**
     * Compute the high digits of m / 5^i / 2^j such that the result is accurate to at least 9
     * decimal digits. i and j are already chosen appropriately.
     * @param m mantissa
     * @param i power of 5
     * @param j power of 2
     * @return high digits of m / 5^i / 2^j
     */
    private static long mulPow5InvDivPow2(final long m, final int i, final int j) {
        // m has at most 55 bits.
        final long mHigh = m >>> 31;
        final long mLow = m & 0x7fffffff;
        final long bits13 = mHigh * POW5_INV_SPLIT[i][0];
        final long bits03 = mLow * POW5_INV_SPLIT[i][0];
        final long bits12 = mHigh * POW5_INV_SPLIT[i][1];
        final long bits02 = mLow * POW5_INV_SPLIT[i][1];
        final long bits11 = mHigh * POW5_INV_SPLIT[i][2];
        final long bits01 = mLow * POW5_INV_SPLIT[i][2];
        final long bits10 = mHigh * POW5_INV_SPLIT[i][3];
        final long bits00 = mLow * POW5_INV_SPLIT[i][3];

        final int actualShift = j - 3 * 31 - 21;
        if (actualShift < 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL,
                                                   j, 3 * 31 + 21);
        }
        return ((((((((bits00 >>> 31) + bits01 + bits10) >>> 31) +
                        bits02 + bits11) >>> 31) +
                        bits03 + bits12) >>> 21) +
                        (bits13 << 10)) >>> actualShift;
    }

}
