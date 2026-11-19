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
package org.hipparchus.util;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * Utility class to manage roman numerals.
 * @since 4.1
 */
public final class RomanNumeral {

    /** Components in decreased value order. */
    private static final Component[] COMPONENTS = {
        new Component( "M", 1000),
        new Component("CM",  900),
        new Component( "D",  500),
        new Component("CD",  400),
        new Component( "C",  100),
        new Component("XC",   90),
        new Component( "L",   50),
        new Component("XL",   40),
        new Component( "X",   10),
        new Component("IX",    9),
        new Component( "V",    5),
        new Component("IV",    4),
        new Component( "I",    1)
    };

    /** Private constructor for a utility class. */
    private RomanNumeral() {
        // nothing to do
    }

    /** Convert an integer to roman numeral.
     * <p>
     * The classical subtractive form is generated here,
     * i.e. 4 is represented as IV and not IIII.
     * </p>
     * @param n number to convert
     * @return string representation as a roman numeral
     */
    public static String toRoman(final int n) {

        checkRange(n);

        // build the string representation from left to right
        final StringBuilder builder = new StringBuilder();
        int remaining = n;
        int componentIndex = -1;
        while (remaining > 0) {
            Component c = COMPONENTS[++componentIndex];
            while (remaining >= c.value) {
                builder.append(c.roman);
                remaining -= c.value;
            }
        }

        return builder.toString();

    }

    /** Convert roman numeral to an integer.
     * <p>
     * We accept here the additive form (IIII instead of IV).
     * </p>
     * @param r roman numeral to convert
     * @return value of the roman numeral
     */
    public static int parse(final String r) {

        // special handling
        if (r.isEmpty()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INVALID_ROMAN_NUMERAL, r);
        }

        int n = 0;

        // parse components
        int stringIndex    = 0;
        int componentIndex = -1;
        while (stringIndex < r.length() && componentIndex < COMPONENTS.length - 1) {
            Component c = COMPONENTS[++componentIndex];
            while (r.substring(stringIndex).startsWith(c.roman)) {
                n           += c.value;
                stringIndex += c.roman.length();
                if (c.roman.length() > 1) {
                    // subtractive components like IV cannot be repeated
                    break;
                }
            }
        }

        if (stringIndex < r.length()) {
            // we were not able to parse the complete string
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INVALID_ROMAN_NUMERAL, r);
        }

        // a posteriori check
        checkRange(n);

        return n;

    }

    /** Check the range of a number.
     * @param n number to check
     */
    private static void checkRange(final int n) {
        if (n < 1 || n > 3999) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.ROMAN_NUMERAL_RANGE, n);
        }
    }

    /** Container for roman numerals. */
    private static class Component {

        /** Roman representation. */
        private final String roman;

        /** Integer value. */
        private final int value;

        /** Simple constructor.
         * @param roman roman representation
         * @param value integer value
         */
        private Component(final String roman, final int value) {
            this.roman = roman;
            this.value = value;
        }

    }

}
