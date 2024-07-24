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

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.Locale;


public class BigFractionFormatTest {

    BigFractionFormat properFormat = null;
    BigFractionFormat improperFormat = null;

    protected Locale getLocale() {
        return Locale.getDefault();
    }

    @BeforeEach
    public void setUp() {
        properFormat = BigFractionFormat.getProperInstance(getLocale());
        improperFormat = BigFractionFormat.getImproperInstance(getLocale());
    }

    @Test
    public void testFormat() {
        BigFraction c = new BigFraction(1, 2);
        String expected = "1 / 2";

        String actual = properFormat.format(c);
        Assertions.assertEquals(expected, actual);

        actual = improperFormat.format(c);
        Assertions.assertEquals(expected, actual);
    }

    @Test
    public void testFormatNegative() {
        BigFraction c = new BigFraction(-1, 2);
        String expected = "-1 / 2";

        String actual = properFormat.format(c);
        Assertions.assertEquals(expected, actual);

        actual = improperFormat.format(c);
        Assertions.assertEquals(expected, actual);
    }

    @Test
    public void testFormatZero() {
        BigFraction c = new BigFraction(0, 1);
        String expected = "0 / 1";

        String actual = properFormat.format(c);
        Assertions.assertEquals(expected, actual);

        actual = improperFormat.format(c);
        Assertions.assertEquals(expected, actual);
    }

    @Test
    public void testFormatImproper() {
        BigFraction c = new BigFraction(5, 3);

        String actual = properFormat.format(c);
        Assertions.assertEquals("1 2 / 3", actual);

        actual = improperFormat.format(c);
        Assertions.assertEquals("5 / 3", actual);
    }

    @Test
    public void testFormatImproperNegative() {
        BigFraction c = new BigFraction(-5, 3);

        String actual = properFormat.format(c);
        Assertions.assertEquals("-1 2 / 3", actual);

        actual = improperFormat.format(c);
        Assertions.assertEquals("-5 / 3", actual);
    }

    @Test
    public void testParse() {
        String source = "1 / 2";

        {
            BigFraction c = properFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(BigInteger.ONE, c.getNumerator());
            Assertions.assertEquals(BigInteger.valueOf(2l), c.getDenominator());

            c = improperFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(BigInteger.ONE, c.getNumerator());
            Assertions.assertEquals(BigInteger.valueOf(2l), c.getDenominator());
        }
    }

    @Test
    public void testParseInteger() {
        String source = "10";
        {
            BigFraction c = properFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(BigInteger.TEN, c.getNumerator());
            Assertions.assertEquals(BigInteger.ONE, c.getDenominator());
        }
        {
            BigFraction c = improperFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(BigInteger.TEN, c.getNumerator());
            Assertions.assertEquals(BigInteger.ONE, c.getDenominator());
        }
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
            BigFraction c = properFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(-1, c.getNumeratorAsInt());
            Assertions.assertEquals(2, c.getDenominatorAsInt());

            c = improperFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(-1, c.getNumeratorAsInt());
            Assertions.assertEquals(2, c.getDenominatorAsInt());

            source = "1 / -2";
            c = properFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(-1, c.getNumeratorAsInt());
            Assertions.assertEquals(2, c.getDenominatorAsInt());

            c = improperFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(-1, c.getNumeratorAsInt());
            Assertions.assertEquals(2, c.getDenominatorAsInt());
        }
    }

    @Test
    public void testParseProper() {
        String source = "1 2 / 3";

        {
            BigFraction c = properFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(5, c.getNumeratorAsInt());
            Assertions.assertEquals(3, c.getDenominatorAsInt());
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
            BigFraction c = properFormat.parse(source);
            Assertions.assertNotNull(c);
            Assertions.assertEquals(-5, c.getNumeratorAsInt());
            Assertions.assertEquals(3, c.getDenominatorAsInt());
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
    public void testParseBig() {
        BigFraction f1 =
            improperFormat.parse("167213075789791382630275400487886041651764456874403" +
                                 " / " +
                                 "53225575123090058458126718248444563466137046489291");
        Assertions.assertEquals(FastMath.PI, f1.doubleValue(), 0.0);
        BigFraction f2 =
            properFormat.parse("3 " +
                               "7536350420521207255895245742552351253353317406530" +
                               " / " +
                               "53225575123090058458126718248444563466137046489291");
        Assertions.assertEquals(FastMath.PI, f2.doubleValue(), 0.0);
        Assertions.assertEquals(f1, f2);
        BigDecimal pi =
            new BigDecimal("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068");
        Assertions.assertEquals(pi, f1.bigDecimalValue(99, RoundingMode.HALF_EVEN));
    }

    @Test
    public void testLongFormat() {
        Assertions.assertEquals("10 / 1", improperFormat.format(10l));
    }

    @Test
    public void testDoubleFormat() {
        Assertions.assertEquals("1 / 16", improperFormat.format(0.0625));
    }
}
