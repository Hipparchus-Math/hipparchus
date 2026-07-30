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
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

class RomanNumeralTest {

    @Test
    public void testZero() {
        doTestWrongRange(0);
    }

    @Test
    public void testNegative() {
        doTestWrongRange(-4);
    }

    @Test
    public void testTooLargeValue() {
        doTestWrongRange(4000);
    }

    @Test
    public void testTooLargeRoman() {
        try {
            RomanNumeral.parse("MMMM");
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.ROMAN_NUMERAL_RANGE, e.getSpecifier());
            Assertions.assertEquals(4000, (Integer) e.getParts()[0]);
        }
    }

    @Test
    public void testComponents() {
        Assertions.assertEquals( "I", RomanNumeral.toRoman(   1));
        Assertions.assertEquals("IV", RomanNumeral.toRoman(   4));
        Assertions.assertEquals( "V", RomanNumeral.toRoman(   5));
        Assertions.assertEquals("IX", RomanNumeral.toRoman(   9));
        Assertions.assertEquals( "X", RomanNumeral.toRoman(  10));
        Assertions.assertEquals("XL", RomanNumeral.toRoman(  40));
        Assertions.assertEquals( "L", RomanNumeral.toRoman(  50));
        Assertions.assertEquals("XC", RomanNumeral.toRoman(  90));
        Assertions.assertEquals( "C", RomanNumeral.toRoman( 100));
        Assertions.assertEquals("CD", RomanNumeral.toRoman( 400));
        Assertions.assertEquals( "D", RomanNumeral.toRoman( 500));
        Assertions.assertEquals("CM", RomanNumeral.toRoman( 900));
        Assertions.assertEquals( "M", RomanNumeral.toRoman(1000));
    }

    @Test
    public void testStandardForm() {
        // reference values from the Wikipedia page
        Assertions.assertEquals(     "XXXIX", RomanNumeral.toRoman(  39));
        Assertions.assertEquals(    "CCXLVI", RomanNumeral.toRoman( 246));
        Assertions.assertEquals( "DCCLXXXIX", RomanNumeral.toRoman( 789));
        Assertions.assertEquals(   "MMCDXXI", RomanNumeral.toRoman(2421));
        Assertions.assertEquals(       "CLX", RomanNumeral.toRoman( 160));
        Assertions.assertEquals(     "CCVII", RomanNumeral.toRoman( 207));
        Assertions.assertEquals(       "MIX", RomanNumeral.toRoman(1009));
        Assertions.assertEquals(     "MLXVI", RomanNumeral.toRoman(1066));
        Assertions.assertEquals( "MMMCMXCIX", RomanNumeral.toRoman(3999));
        Assertions.assertEquals(  "MCMXVIII", RomanNumeral.toRoman(1918));
        Assertions.assertEquals(   "MCMXLIV", RomanNumeral.toRoman(1944));
        Assertions.assertEquals(    "MMXXVI", RomanNumeral.toRoman(2026));
    }

    @Test
    public void testOtherAdditiveForm() {
        // reference values from the Wikipedia page
        Assertions.assertEquals(   4, RomanNumeral.parse(     "IIII"));
        Assertions.assertEquals(  40, RomanNumeral.parse(     "XXXX"));
        Assertions.assertEquals( 400, RomanNumeral.parse(     "CCCC"));
        Assertions.assertEquals(  24, RomanNumeral.parse(   "XXIIII"));
        Assertions.assertEquals(  74, RomanNumeral.parse(  "LXXIIII"));
        Assertions.assertEquals( 490, RomanNumeral.parse("CCCCLXXXX"));
        Assertions.assertEquals(   9, RomanNumeral.parse(    "VIIII"));
        Assertions.assertEquals(  90, RomanNumeral.parse(    "LXXXX"));
        Assertions.assertEquals( 900, RomanNumeral.parse(    "DCCCC"));
        Assertions.assertEquals(  44, RomanNumeral.parse(   "XLIIII"));
        Assertions.assertEquals(1910, RomanNumeral.parse(  "MDCCCCX"));
        Assertions.assertEquals(1903, RomanNumeral.parse(  "MDCDIII"));
    }

    @Test
    public void testEmpty() {
        doTestInvalid("");
    }

    @Test
    public void testNoRepeatedSubtraction() {
        doTestInvalid("IVIV");
    }

    @Test
    public void testRoundTrip() {
        for (int i = 1; i < 4000; ++i) {
            Assertions.assertEquals(i, RomanNumeral.parse(RomanNumeral.toRoman(i)));
        }
    }

    private void doTestWrongRange(final int n) {
        try {
            RomanNumeral.toRoman(n);
            Assertions.fail("an exception should have been thrown");
        } catch (final MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.ROMAN_NUMERAL_RANGE, e.getSpecifier());
            Assertions.assertEquals(n, (Integer) e.getParts()[0]);
        }
    }

    private void doTestInvalid(final String r) {
        try {
            RomanNumeral.parse(r);
            Assertions.fail("an exception should have been thrown");
        } catch (final MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.INVALID_ROMAN_NUMERAL, e.getSpecifier());
            Assertions.assertEquals(r, e.getParts()[0]);
        }
    }

}
