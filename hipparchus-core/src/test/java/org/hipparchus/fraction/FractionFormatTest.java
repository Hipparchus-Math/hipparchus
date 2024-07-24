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

import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Locale;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.fail;


class FractionFormatTest {

    FractionFormat properFormat = null;
    FractionFormat improperFormat = null;

    protected Locale getLocale() {
        return Locale.getDefault();
    }

    @BeforeEach
    void setUp() {
        properFormat = FractionFormat.getProperInstance(getLocale());
        improperFormat = FractionFormat.getImproperInstance(getLocale());
    }

    @Test
    void testFormat() {
        Fraction c = new Fraction(1, 2);
        String expected = "1 / 2";

        String actual = properFormat.format(c);
        assertEquals(expected, actual);

        actual = improperFormat.format(c);
        assertEquals(expected, actual);
    }

    @Test
    void testFormatNegative() {
        Fraction c = new Fraction(-1, 2);
        String expected = "-1 / 2";

        String actual = properFormat.format(c);
        assertEquals(expected, actual);

        actual = improperFormat.format(c);
        assertEquals(expected, actual);
    }

    @Test
    void testFormatZero() {
        Fraction c = new Fraction(0, 1);
        String expected = "0 / 1";

        String actual = properFormat.format(c);
        assertEquals(expected, actual);

        actual = improperFormat.format(c);
        assertEquals(expected, actual);
    }

    @Test
    void testFormatImproper() {
        Fraction c = new Fraction(5, 3);

        String actual = properFormat.format(c);
        assertEquals("1 2 / 3", actual);

        actual = improperFormat.format(c);
        assertEquals("5 / 3", actual);
    }

    @Test
    void testFormatImproperNegative() {
        Fraction c = new Fraction(-5, 3);

        String actual = properFormat.format(c);
        assertEquals("-1 2 / 3", actual);

        actual = improperFormat.format(c);
        assertEquals("-5 / 3", actual);
    }

    @Test
    void testParse() {
        String source = "1 / 2";

        try {
            Fraction c = properFormat.parse(source);
            assertNotNull(c);
            assertEquals(1, c.getNumerator());
            assertEquals(2, c.getDenominator());

            c = improperFormat.parse(source);
            assertNotNull(c);
            assertEquals(1, c.getNumerator());
            assertEquals(2, c.getDenominator());
        } catch (MathIllegalStateException ex) {
            fail(ex.getMessage());
        }
    }

    @Test
    void testParseInteger() {
        String source = "10";
        {
            Fraction c = properFormat.parse(source);
            assertNotNull(c);
            assertEquals(10, c.getNumerator());
            assertEquals(1, c.getDenominator());
        }
        {
            Fraction c = improperFormat.parse(source);
            assertNotNull(c);
            assertEquals(10, c.getNumerator());
            assertEquals(1, c.getDenominator());
        }
    }

    @Test
    void testParseOne1() {
        String source = "1 / 1";
        Fraction c = properFormat.parse(source);
        assertNotNull(c);
        assertEquals(1, c.getNumerator());
        assertEquals(1, c.getDenominator());
    }

    @Test
    void testParseOne2() {
        String source = "10 / 10";
        Fraction c = properFormat.parse(source);
        assertNotNull(c);
        assertEquals(1, c.getNumerator());
        assertEquals(1, c.getDenominator());
    }

    @Test
    void testParseZero1() {
        String source = "0 / 1";
        Fraction c = properFormat.parse(source);
        assertNotNull(c);
        assertEquals(0, c.getNumerator());
        assertEquals(1, c.getDenominator());
    }

    @Test
    void testParseZero2() {
        String source = "-0 / 1";
        Fraction c = properFormat.parse(source);
        assertNotNull(c);
        assertEquals(0, c.getNumerator());
        assertEquals(1, c.getDenominator());
        // This test shows that the sign is not preserved.
        assertEquals(Double.POSITIVE_INFINITY, 1d / c.doubleValue(), 0);
    }

    @Test
    void testParseInvalid() {
        String source = "a";
        String msg = "should not be able to parse '10 / a'.";
        try {
            properFormat.parse(source);
            fail(msg);
        } catch (MathIllegalStateException ex) {
            // success
        }
        try {
            improperFormat.parse(source);
            fail(msg);
        } catch (MathIllegalStateException ex) {
            // success
        }
    }

    @Test
    void testParseInvalidDenominator() {
        String source = "10 / a";
        String msg = "should not be able to parse '10 / a'.";
        try {
            properFormat.parse(source);
            fail(msg);
        } catch (MathIllegalStateException ex) {
            // success
        }
        try {
            improperFormat.parse(source);
            fail(msg);
        } catch (MathIllegalStateException ex) {
            // success
        }
    }

    @Test
    void testParseNegative() {

        {
            String source = "-1 / 2";
            Fraction c = properFormat.parse(source);
            assertNotNull(c);
            assertEquals(-1, c.getNumerator());
            assertEquals(2, c.getDenominator());

            c = improperFormat.parse(source);
            assertNotNull(c);
            assertEquals(-1, c.getNumerator());
            assertEquals(2, c.getDenominator());

            source = "1 / -2";
            c = properFormat.parse(source);
            assertNotNull(c);
            assertEquals(-1, c.getNumerator());
            assertEquals(2, c.getDenominator());

            c = improperFormat.parse(source);
            assertNotNull(c);
            assertEquals(-1, c.getNumerator());
            assertEquals(2, c.getDenominator());
        }
    }

    @Test
    void testParseProper() {
        String source = "1 2 / 3";

        {
            Fraction c = properFormat.parse(source);
            assertNotNull(c);
            assertEquals(5, c.getNumerator());
            assertEquals(3, c.getDenominator());
        }

        try {
            improperFormat.parse(source);
            fail("invalid improper fraction.");
        } catch (MathIllegalStateException ex) {
            // success
        }
    }

    @Test
    void testParseProperNegative() {
        String source = "-1 2 / 3";
        {
            Fraction c = properFormat.parse(source);
            assertNotNull(c);
            assertEquals(-5, c.getNumerator());
            assertEquals(3, c.getDenominator());
        }

        try {
            improperFormat.parse(source);
            fail("invalid improper fraction.");
        } catch (MathIllegalStateException ex) {
            // success
        }
    }

    @Test
    void testParseProperInvalidMinus() {
        String source = "2 -2 / 3";
        try {
            properFormat.parse(source);
            fail("invalid minus in improper fraction.");
        } catch (MathIllegalStateException ex) {
            // expected
        }
        source = "2 2 / -3";
        try {
            properFormat.parse(source);
            fail("invalid minus in improper fraction.");
        } catch (MathIllegalStateException ex) {
            // expected
        }
    }

    @Test
    void testLongFormat() {
        assertEquals("10 / 1", improperFormat.format(10l));
    }

    @Test
    void testDoubleFormat() {
        assertEquals("355 / 113", improperFormat.format(FastMath.PI));
    }
}
