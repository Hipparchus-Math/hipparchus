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
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Locale;


public class FractionFormatTest {

    FractionFormat properFormat = null;
    FractionFormat improperFormat = null;

    protected Locale getLocale() {
        return Locale.getDefault();
    }

    @BeforeEach
    public void setUp() {
        properFormat = FractionFormat.getProperInstance(getLocale());
        improperFormat = FractionFormat.getImproperInstance(getLocale());
    }

    @Test
    public void testFormat() {
        Fraction c = new Fraction(1, 2);
        String expected = "1 / 2";

        String actual = properFormat.format(c);
        Assertions.assertEquals(expected, actual);

        actual = improperFormat.format(c);
        Assertions.assertEquals(expected, actual);
    }

    @Test
    public void testFormatNegative() {
        Fraction c = new Fraction(-1, 2);
        String expected = "-1 / 2";

        String actual = properFormat.format(c);
        Assertions.assertEquals(expected, actual);

        actual = improperFormat.format(c);
        Assertions.assertEquals(expected, actual);
    }

    @Test
    public void testFormatZero() {
        Fraction c = new Fraction(0, 1);
        String expected = "0 / 1";

        String actual = properFormat.format(c);
        Assertions.assertEquals(expected, actual);

        actual = improperFormat.format(c);
        Assertions.assertEquals(expected, actual);
    }

    @Test
    public void testFormatImproper() {
        Fraction c = new Fraction(5, 3);

        String actual = properFormat.format(c);
        Assertions.assertEquals("1 2 / 3", actual);

        actual = improperFormat.format(c);
        Assertions.assertEquals("5 / 3", actual);
    }

    @Test
    public void testFormatImproperNegative() {
        Fraction c = new Fraction(-5, 3);

        String actual = properFormat.format(c);
        Assertions.assertEquals("-1 2 / 3", actual);

        actual = improperFormat.format(c);
        Assertions.assertEquals("-5 / 3", actual);
    }

    @Test
    public void testParse() {
        String source = "1 / 2";

        try {
            Fraction c = properFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(1, c.getNumerator());
            Assertions.assertEquals(2, c.getDenominator());

            c = improperFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(1, c.getNumerator());
            Assertions.assertEquals(2, c.getDenominator());
        } catch (MathIllegalStateException ex) {
            Assertions.fail(ex.getMessage());
        }
    }

    @Test
    public void testParseInteger() {
        String source = "10";
        {
            Fraction c = properFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(10, c.getNumerator());
            Assertions.assertEquals(1, c.getDenominator());
        }
        {
            Fraction c = improperFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(10, c.getNumerator());
            Assertions.assertEquals(1, c.getDenominator());
        }
    }

    @Test
    public void testParseOne1() {
        String source = "1 / 1";
        Fraction c = properFormat.parse(source);
        Assertions.assertNotNull(c);
        Assertions.assertEquals(1, c.getNumerator());
        Assertions.assertEquals(1, c.getDenominator());
    }

    @Test
    public void testParseOne2() {
        String source = "10 / 10";
        Fraction c = properFormat.parse(source);
        Assertions.assertNotNull(c);
        Assertions.assertEquals(1, c.getNumerator());
        Assertions.assertEquals(1, c.getDenominator());
    }

    @Test
    public void testParseZero1() {
        String source = "0 / 1";
        Fraction c = properFormat.parse(source);
        Assertions.assertNotNull(c);
        Assertions.assertEquals(0, c.getNumerator());
        Assertions.assertEquals(1, c.getDenominator());
    }

    @Test
    public void testParseZero2() {
        String source = "-0 / 1";
        Fraction c = properFormat.parse(source);
        Assertions.assertNotNull(c);
        Assertions.assertEquals(0, c.getNumerator());
        Assertions.assertEquals(1, c.getDenominator());
        // This test shows that the sign is not preserved.
        Assertions.assertEquals(Double.POSITIVE_INFINITY, 1d / c.doubleValue(), 0);
    }

    @Test
    public void testParseInvalid() {
        String source = "a";
        String msg = "should not be able to parse '10 / a'.";
        try {
            properFormat.parse(source);
            Assertions.fail(msg);
        } catch (MathIllegalStateException ex) {
            // success
        }
        try {
            improperFormat.parse(source);
            Assertions.fail(msg);
        } catch (MathIllegalStateException ex) {
            // success
        }
    }

    @Test
    public void testParseInvalidDenominator() {
        String source = "10 / a";
        String msg = "should not be able to parse '10 / a'.";
        try {
            properFormat.parse(source);
            Assertions.fail(msg);
        } catch (MathIllegalStateException ex) {
            // success
        }
        try {
            improperFormat.parse(source);
            Assertions.fail(msg);
        } catch (MathIllegalStateException ex) {
            // success
        }
    }

    @Test
    public void testParseNegative() {

        {
            String source = "-1 / 2";
            Fraction c = properFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(-1, c.getNumerator());
            Assertions.assertEquals(2, c.getDenominator());

            c = improperFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(-1, c.getNumerator());
            Assertions.assertEquals(2, c.getDenominator());

            source = "1 / -2";
            c = properFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(-1, c.getNumerator());
            Assertions.assertEquals(2, c.getDenominator());

            c = improperFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(-1, c.getNumerator());
            Assertions.assertEquals(2, c.getDenominator());
        }
    }

    @Test
    public void testParseProper() {
        String source = "1 2 / 3";

        {
            Fraction c = properFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(5, c.getNumerator());
            Assertions.assertEquals(3, c.getDenominator());
        }

        try {
            improperFormat.parse(source);
            Assertions.fail("invalid improper fraction.");
        } catch (MathIllegalStateException ex) {
            // success
        }
    }

    @Test
    public void testParseProperNegative() {
        String source = "-1 2 / 3";
        {
            Fraction c = properFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(-5, c.getNumerator());
            Assertions.assertEquals(3, c.getDenominator());
        }

        try {
            improperFormat.parse(source);
            Assertions.fail("invalid improper fraction.");
        } catch (MathIllegalStateException ex) {
            // success
        }
    }

    @Test
    public void testParseProperInvalidMinus() {
        String source = "2 -2 / 3";
        try {
            properFormat.parse(source);
            Assertions.fail("invalid minus in improper fraction.");
        } catch (MathIllegalStateException ex) {
            // expected
        }
        source = "2 2 / -3";
        try {
            properFormat.parse(source);
            Assertions.fail("invalid minus in improper fraction.");
        } catch (MathIllegalStateException ex) {
            // expected
        }
    }

    @Test
    public void testLongFormat() {
        Assertions.assertEquals("10 / 1", improperFormat.format(10l));
    }

    @Test
    public void testDoubleFormat() {
        Assertions.assertEquals("355 / 113", improperFormat.format(FastMath.PI));
    }
}
