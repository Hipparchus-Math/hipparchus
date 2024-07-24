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

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.stream.Collectors;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;


/**
 */
class FractionTest {

    private void customAssertFraction(int expectedNumerator, int expectedDenominator, Fraction actual) {
        assertEquals(expectedNumerator, actual.getNumerator());
        assertEquals(expectedDenominator, actual.getDenominator());
    }

    @Test
    void testConstructor() {
        customAssertFraction(0, 1, new Fraction(0, 1));
        customAssertFraction(0, 1, new Fraction(0, 2));
        customAssertFraction(0, 1, new Fraction(0, -1));
        customAssertFraction(1, 2, new Fraction(1, 2));
        customAssertFraction(1, 2, new Fraction(2, 4));
        customAssertFraction(-1, 2, new Fraction(-1, 2));
        customAssertFraction(-1, 2, new Fraction(1, -2));
        customAssertFraction(-1, 2, new Fraction(-2, 4));
        customAssertFraction(-1, 2, new Fraction(2, -4));

        // overflow
        try {
            new Fraction(Integer.MIN_VALUE, -1);
            fail();
        } catch (MathRuntimeException ex) {
            // success
        }
        try {
            new Fraction(1, Integer.MIN_VALUE);
            fail();
        } catch (MathRuntimeException ex) {
            // success
        }

        customAssertFraction(0, 1, new Fraction(0.00000000000001));
        customAssertFraction(2, 5, new Fraction(0.40000000000001));
        customAssertFraction(15, 1, new Fraction(15.0000000000001));
    }

    @Test
    void testIsInteger() {
        assertTrue(new Fraction(12, 12).isInteger());
        assertTrue(new Fraction(14, 7).isInteger());
        assertFalse(new Fraction(12, 11).isInteger());
    }

    @Test
    void testGoldenRatio() {
        assertThrows(MathIllegalStateException.class, () -> {
            // the golden ratio is notoriously a difficult number for continuous fraction
            new Fraction((1 + FastMath.sqrt(5)) / 2, 1.0e-12, 25);
        });
    }

    // MATH-179
    @Test
    void testDoubleConstructor() {
        customAssertFraction(1, 2, new Fraction((double)1 / (double)2));
        customAssertFraction(1, 3, new Fraction((double)1 / (double)3));
        customAssertFraction(2, 3, new Fraction((double)2 / (double)3));
        customAssertFraction(1, 4, new Fraction((double)1 / (double)4));
        customAssertFraction(3, 4, new Fraction((double)3 / (double)4));
        customAssertFraction(1, 5, new Fraction((double)1 / (double)5));
        customAssertFraction(2, 5, new Fraction((double)2 / (double)5));
        customAssertFraction(3, 5, new Fraction((double)3 / (double)5));
        customAssertFraction(4, 5, new Fraction((double)4 / (double)5));
        customAssertFraction(1, 6, new Fraction((double)1 / (double)6));
        customAssertFraction(5, 6, new Fraction((double)5 / (double)6));
        customAssertFraction(1, 7, new Fraction((double)1 / (double)7));
        customAssertFraction(2, 7, new Fraction((double)2 / (double)7));
        customAssertFraction(3, 7, new Fraction((double)3 / (double)7));
        customAssertFraction(4, 7, new Fraction((double)4 / (double)7));
        customAssertFraction(5, 7, new Fraction((double)5 / (double)7));
        customAssertFraction(6, 7, new Fraction((double)6 / (double)7));
        customAssertFraction(1, 8, new Fraction((double)1 / (double)8));
        customAssertFraction(3, 8, new Fraction((double)3 / (double)8));
        customAssertFraction(5, 8, new Fraction((double)5 / (double)8));
        customAssertFraction(7, 8, new Fraction((double)7 / (double)8));
        customAssertFraction(1, 9, new Fraction((double)1 / (double)9));
        customAssertFraction(2, 9, new Fraction((double)2 / (double)9));
        customAssertFraction(4, 9, new Fraction((double)4 / (double)9));
        customAssertFraction(5, 9, new Fraction((double)5 / (double)9));
        customAssertFraction(7, 9, new Fraction((double)7 / (double)9));
        customAssertFraction(8, 9, new Fraction((double)8 / (double)9));
        customAssertFraction(1, 10, new Fraction((double)1 / (double)10));
        customAssertFraction(3, 10, new Fraction((double)3 / (double)10));
        customAssertFraction(7, 10, new Fraction((double)7 / (double)10));
        customAssertFraction(9, 10, new Fraction((double)9 / (double)10));
        customAssertFraction(1, 11, new Fraction((double)1 / (double)11));
        customAssertFraction(2, 11, new Fraction((double)2 / (double)11));
        customAssertFraction(3, 11, new Fraction((double)3 / (double)11));
        customAssertFraction(4, 11, new Fraction((double)4 / (double)11));
        customAssertFraction(5, 11, new Fraction((double)5 / (double)11));
        customAssertFraction(6, 11, new Fraction((double)6 / (double)11));
        customAssertFraction(7, 11, new Fraction((double)7 / (double)11));
        customAssertFraction(8, 11, new Fraction((double)8 / (double)11));
        customAssertFraction(9, 11, new Fraction((double)9 / (double)11));
        customAssertFraction(10, 11, new Fraction((double)10 / (double)11));
    }

    // MATH-181
    @Test
    void testDigitLimitConstructor() {
        customAssertFraction(2, 5, new Fraction(0.4, 9));
        customAssertFraction(2, 5, new Fraction(0.4, 99));
        customAssertFraction(2, 5, new Fraction(0.4, 999));

        customAssertFraction(3, 5, new Fraction(0.6152, 9));
        customAssertFraction(8, 13, new Fraction(0.6152, 99));
        customAssertFraction(510, 829, new Fraction(0.6152, 999));
        customAssertFraction(769, 1250, new Fraction(0.6152, 9999));

        // MATH-996
        customAssertFraction(1, 2, new Fraction(0.5000000001, 10));
    }

    @Test
    void testIntegerOverflow() {
        checkIntegerOverflow(0.75000000001455192);
        checkIntegerOverflow(1.0e10);
        checkIntegerOverflow(-1.0e10);
    }

    @Test
    void testSignum() {
        assertEquals(-1, new Fraction(4, -5).signum());
        assertEquals(-1, new Fraction(-4, 5).signum());
        assertEquals( 0, new Fraction(0).signum());
        assertEquals(+1, new Fraction(-4, -5).signum());
        assertEquals(+1, new Fraction(4, 5).signum());
    }

    private void checkIntegerOverflow(double a) {
        assertThrows(MathIllegalStateException.class, () -> new Fraction(a, 1.0e-12, 1000));
    }

    @Test
    void testEpsilonLimitConstructor() {
        customAssertFraction(2, 5, new Fraction(0.4, 1.0e-5, 100));

        customAssertFraction(3, 5, new Fraction(0.6152, 0.02, 100));
        customAssertFraction(8, 13, new Fraction(0.6152, 1.0e-3, 100));
        customAssertFraction(251, 408, new Fraction(0.6152, 1.0e-4, 100));
        customAssertFraction(251, 408, new Fraction(0.6152, 1.0e-5, 100));
        customAssertFraction(510, 829, new Fraction(0.6152, 1.0e-6, 100));
        customAssertFraction(769, 1250, new Fraction(0.6152, 1.0e-7, 100));
    }

    @Test
    void testCompareTo() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(1, 3);
        Fraction third = new Fraction(1, 2);

        assertEquals(0, first.compareTo(first));
        assertEquals(0, first.compareTo(third));
        assertEquals(1, first.compareTo(second));
        assertEquals(-1, second.compareTo(first));

        // these two values are different approximations of PI
        // the first  one is approximately PI - 3.07e-18
        // the second one is approximately PI + 1.936e-17
        Fraction pi1 = new Fraction(1068966896, 340262731);
        Fraction pi2 = new Fraction( 411557987, 131002976);
        assertEquals(-1, pi1.compareTo(pi2));
        assertEquals( 1, pi2.compareTo(pi1));
        assertEquals(0.0, pi1.doubleValue() - pi2.doubleValue(), 1.0e-20);
    }

    @Test
    void testDoubleValue() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(1, 3);

        assertEquals(0.5, first.doubleValue(), 0.0);
        assertEquals(1.0 / 3.0, second.doubleValue(), 0.0);
    }

    @Test
    void testFloatValue() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(1, 3);

        assertEquals(0.5f, first.floatValue(), 0.0f);
        assertEquals((float)(1.0 / 3.0), second.floatValue(), 0.0f);
    }

    @Test
    void testIntValue() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(3, 2);

        assertEquals(0, first.intValue());
        assertEquals(1, second.intValue());
    }

    @Test
    void testLongValue() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(3, 2);

        assertEquals(0L, first.longValue());
        assertEquals(1L, second.longValue());
    }

    @Test
    void testConstructorDouble() {
        customAssertFraction(1, 2, new Fraction(0.5));
        customAssertFraction(1, 3, new Fraction(1.0 / 3.0));
        customAssertFraction(17, 100, new Fraction(17.0 / 100.0));
        customAssertFraction(317, 100, new Fraction(317.0 / 100.0));
        customAssertFraction(-1, 2, new Fraction(-0.5));
        customAssertFraction(-1, 3, new Fraction(-1.0 / 3.0));
        customAssertFraction(-17, 100, new Fraction(17.0 / -100.0));
        customAssertFraction(-317, 100, new Fraction(-317.0 / 100.0));
    }

    @Test
    void testAbs() {
        Fraction a = new Fraction(10, 21);
        Fraction b = new Fraction(-10, 21);
        Fraction c = new Fraction(10, -21);

        customAssertFraction(10, 21, a.abs());
        customAssertFraction(10, 21, b.abs());
        customAssertFraction(10, 21, c.abs());
    }

    @Test
    void testPercentage() {
        assertEquals(50.0, new Fraction(1, 2).percentageValue(), 1.0e-15);
    }

    @Test
    void testMath835() {
        final int numer = Integer.MAX_VALUE / 99;
        final int denom = 1;
        final double percentage = 100 * ((double) numer) / denom;
        final Fraction frac = new Fraction(numer, denom);
        // With the implementation that preceded the fix suggested in MATH-835,
        // this test was failing, due to overflow.
        assertEquals(percentage, frac.percentageValue(), Math.ulp(percentage));
    }

    @Test
    void testMath1261() {
        final Fraction a = new Fraction(Integer.MAX_VALUE, 2);
        final Fraction b = a.multiply(2);
        assertEquals(b, new Fraction(Integer.MAX_VALUE));

        final Fraction c = new Fraction(2, Integer.MAX_VALUE);
        final Fraction d = c.divide(2);
        assertEquals(d, new Fraction(1, Integer.MAX_VALUE));
    }

    @Test
    void testReciprocal() {
        Fraction f = null;

        f = new Fraction(50, 75);
        f = f.reciprocal();
        assertEquals(3, f.getNumerator());
        assertEquals(2, f.getDenominator());

        f = new Fraction(4, 3);
        f = f.reciprocal();
        assertEquals(3, f.getNumerator());
        assertEquals(4, f.getDenominator());

        f = new Fraction(-15, 47);
        f = f.reciprocal();
        assertEquals(-47, f.getNumerator());
        assertEquals(15, f.getDenominator());

        f = new Fraction(0, 3);
        try {
            f = f.reciprocal();
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        // large values
        f = new Fraction(Integer.MAX_VALUE, 1);
        f = f.reciprocal();
        assertEquals(1, f.getNumerator());
        assertEquals(Integer.MAX_VALUE, f.getDenominator());
    }

    @Test
    void testNegate() {
        Fraction f = null;

        f = new Fraction(50, 75);
        f = f.negate();
        assertEquals(-2, f.getNumerator());
        assertEquals(3, f.getDenominator());

        f = new Fraction(-50, 75);
        f = f.negate();
        assertEquals(2, f.getNumerator());
        assertEquals(3, f.getDenominator());

        // large values
        f = new Fraction(Integer.MAX_VALUE-1, Integer.MAX_VALUE);
        f = f.negate();
        assertEquals(Integer.MIN_VALUE+2, f.getNumerator());
        assertEquals(Integer.MAX_VALUE, f.getDenominator());

        f = new Fraction(Integer.MIN_VALUE, 1);
        try {
            f = f.negate();
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}
    }

    @Test
    void testAdd() {
        Fraction a = new Fraction(1, 2);
        Fraction b = new Fraction(2, 3);

        customAssertFraction(1, 1, a.add(a));
        customAssertFraction(7, 6, a.add(b));
        customAssertFraction(7, 6, b.add(a));
        customAssertFraction(4, 3, b.add(b));

        Fraction f1 = new Fraction(Integer.MAX_VALUE - 1, 1);
        Fraction f2 = Fraction.ONE;
        Fraction f = f1.add(f2);
        assertEquals(Integer.MAX_VALUE, f.getNumerator());
        assertEquals(1, f.getDenominator());
        f = f1.add(1);
        assertEquals(Integer.MAX_VALUE, f.getNumerator());
        assertEquals(1, f.getDenominator());

        f1 = new Fraction(-1, 13*13*2*2);
        f2 = new Fraction(-2, 13*17*2);
        f = f1.add(f2);
        assertEquals(13*13*17*2*2, f.getDenominator());
        assertEquals(-17 - 2*13*2, f.getNumerator());

        try {
            f.add(null);
            fail("expecting NullArgumentException");
        } catch (NullArgumentException ex) {}

        // if this fraction is added naively, it will overflow.
        // check that it doesn't.
        f1 = new Fraction(1,32768*3);
        f2 = new Fraction(1,59049);
        f = f1.add(f2);
        assertEquals(52451, f.getNumerator());
        assertEquals(1934917632, f.getDenominator());

        f1 = new Fraction(Integer.MIN_VALUE, 3);
        f2 = new Fraction(1,3);
        f = f1.add(f2);
        assertEquals(Integer.MIN_VALUE+1, f.getNumerator());
        assertEquals(3, f.getDenominator());

        f1 = new Fraction(Integer.MAX_VALUE - 1, 1);
        f2 = Fraction.ONE;
        f = f1.add(f2);
        assertEquals(Integer.MAX_VALUE, f.getNumerator());
        assertEquals(1, f.getDenominator());

        try {
            f = f.add(Fraction.ONE); // should overflow
            fail("expecting MathRuntimeException but got: " + f.toString());
        } catch (MathRuntimeException ex) {}

        // denominator should not be a multiple of 2 or 3 to trigger overflow
        f1 = new Fraction(Integer.MIN_VALUE, 5);
        f2 = new Fraction(-1,5);
        try {
            f = f1.add(f2); // should overflow
            fail("expecting MathRuntimeException but got: " + f.toString());
        } catch (MathRuntimeException ex) {}

        try {
            f= new Fraction(-Integer.MAX_VALUE, 1);
            f = f.add(f);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        try {
            f= new Fraction(-Integer.MAX_VALUE, 1);
            f = f.add(f);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        f1 = new Fraction(3,327680);
        f2 = new Fraction(2,59049);
        try {
            f = f1.add(f2); // should overflow
            fail("expecting MathRuntimeException but got: " + f.toString());
        } catch (MathRuntimeException ex) {}
    }

    @Test
    void testDivide() {
        Fraction a = new Fraction(1, 2);
        Fraction b = new Fraction(2, 3);

        customAssertFraction(1, 1, a.divide(a));
        customAssertFraction(3, 4, a.divide(b));
        customAssertFraction(4, 3, b.divide(a));
        customAssertFraction(1, 1, b.divide(b));

        Fraction f1 = new Fraction(3, 5);
        Fraction f2 = Fraction.ZERO;
        try {
            f1.divide(f2);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        f1 = new Fraction(0, 5);
        f2 = new Fraction(2, 7);
        Fraction f = f1.divide(f2);
        assertSame(Fraction.ZERO, f);

        f1 = new Fraction(2, 7);
        f2 = Fraction.ONE;
        f = f1.divide(f2);
        assertEquals(2, f.getNumerator());
        assertEquals(7, f.getDenominator());

        f1 = new Fraction(1, Integer.MAX_VALUE);
        f = f1.divide(f1);
        assertEquals(1, f.getNumerator());
        assertEquals(1, f.getDenominator());

        f1 = new Fraction(Integer.MIN_VALUE, Integer.MAX_VALUE);
        f2 = new Fraction(1, Integer.MAX_VALUE);
        f = f1.divide(f2);
        assertEquals(Integer.MIN_VALUE, f.getNumerator());
        assertEquals(1, f.getDenominator());

        try {
            f.divide(null);
            fail("NullArgumentException");
        } catch (NullArgumentException ex) {}

        try {
            f1 = new Fraction(1, Integer.MAX_VALUE);
            f = f1.divide(f1.reciprocal());  // should overflow
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}
        try {
            f1 = new Fraction(1, -Integer.MAX_VALUE);
            f = f1.divide(f1.reciprocal());  // should overflow
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        f1 = new Fraction(6, 35);
        f  = f1.divide(15);
        assertEquals(2, f.getNumerator());
        assertEquals(175, f.getDenominator());

    }

    @Test
    void testMultiply() {
        Fraction a = new Fraction(1, 2);
        Fraction b = new Fraction(2, 3);

        customAssertFraction(1, 4, a.multiply(a));
        customAssertFraction(1, 3, a.multiply(b));
        customAssertFraction(1, 3, b.multiply(a));
        customAssertFraction(4, 9, b.multiply(b));

        Fraction f1 = new Fraction(Integer.MAX_VALUE, 1);
        Fraction f2 = new Fraction(Integer.MIN_VALUE, Integer.MAX_VALUE);
        Fraction f = f1.multiply(f2);
        assertEquals(Integer.MIN_VALUE, f.getNumerator());
        assertEquals(1, f.getDenominator());

        try {
            f.multiply(null);
            fail("expecting NullArgumentException");
        } catch (NullArgumentException ex) {}

        f1 = new Fraction(6, 35);
        f  = f1.multiply(15);
        assertEquals(18, f.getNumerator());
        assertEquals(7, f.getDenominator());
    }

    @Test
    void testSubtract() {
        Fraction a = new Fraction(1, 2);
        Fraction b = new Fraction(2, 3);

        customAssertFraction(0, 1, a.subtract(a));
        customAssertFraction(-1, 6, a.subtract(b));
        customAssertFraction(1, 6, b.subtract(a));
        customAssertFraction(0, 1, b.subtract(b));

        Fraction f = new Fraction(1,1);
        try {
            f.subtract(null);
            fail("expecting NullArgumentException");
        } catch (NullArgumentException ex) {}

        // if this fraction is subtracted naively, it will overflow.
        // check that it doesn't.
        Fraction f1 = new Fraction(1,32768*3);
        Fraction f2 = new Fraction(1,59049);
        f = f1.subtract(f2);
        assertEquals(-13085, f.getNumerator());
        assertEquals(1934917632, f.getDenominator());

        f1 = new Fraction(Integer.MIN_VALUE, 3);
        f2 = new Fraction(1,3).negate();
        f = f1.subtract(f2);
        assertEquals(Integer.MIN_VALUE+1, f.getNumerator());
        assertEquals(3, f.getDenominator());

        f1 = new Fraction(Integer.MAX_VALUE, 1);
        f2 = Fraction.ONE;
        f = f1.subtract(f2);
        assertEquals(Integer.MAX_VALUE-1, f.getNumerator());
        assertEquals(1, f.getDenominator());
        f = f1.subtract(1);
        assertEquals(Integer.MAX_VALUE-1, f.getNumerator());
        assertEquals(1, f.getDenominator());

        try {
            f1 = new Fraction(1, Integer.MAX_VALUE);
            f2 = new Fraction(1, Integer.MAX_VALUE - 1);
            f = f1.subtract(f2);
            fail("expecting MathRuntimeException");  //should overflow
        } catch (MathRuntimeException ex) {}

        // denominator should not be a multiple of 2 or 3 to trigger overflow
        f1 = new Fraction(Integer.MIN_VALUE, 5);
        f2 = new Fraction(1,5);
        try {
            f = f1.subtract(f2); // should overflow
            fail("expecting MathRuntimeException but got: " + f.toString());
        } catch (MathRuntimeException ex) {}

        try {
            f= new Fraction(Integer.MIN_VALUE, 1);
            f = f.subtract(Fraction.ONE);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        try {
            f= new Fraction(Integer.MAX_VALUE, 1);
            f = f.subtract(Fraction.ONE.negate());
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        f1 = new Fraction(3,327680);
        f2 = new Fraction(2,59049);
        try {
            f = f1.subtract(f2); // should overflow
            fail("expecting MathRuntimeException but got: " + f.toString());
        } catch (MathRuntimeException ex) {}
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testEqualsAndHashCode() {
        Fraction zero  = new Fraction(0,1);
        Fraction nullFraction = null;
        assertEquals(zero, zero);
        assertNotEquals(zero, nullFraction);
        assertNotEquals(zero, Double.valueOf(0));
        Fraction zero2 = new Fraction(0,2);
        assertEquals(zero, zero2);
        assertEquals(zero.hashCode(), zero2.hashCode());
        Fraction one = new Fraction(1,1);
        assertFalse((one.equals(zero) ||zero.equals(one)));
    }

    @Test
    void testGCD() {
      Fraction first = new Fraction(1, 3);
      Fraction second = new Fraction(2, 5);
      Fraction third = new Fraction(3, 7);
      Fraction gcd1 = first.gcd(second);
        assertEquals(gcd1, Fraction.getReducedFraction(1, 15));
      Fraction gcd2 = gcd1.gcd(third);
        assertEquals(gcd2, Fraction.getReducedFraction(1, 105));

      // example from https://math.stackexchange.com/a/151089
      Fraction x = new Fraction(3, 7);
      Fraction y = new Fraction(12, 22);
      Fraction gcd = x.gcd(y);
        assertEquals(gcd, Fraction.getReducedFraction(3, 77));

      x = new Fraction(13, 6);
      y = new Fraction(3, 4);
      gcd = x.gcd(y);
        assertEquals(gcd, Fraction.getReducedFraction(1, 12));

    }

    @Test
    void testLCM() {
      Fraction first = new Fraction(1, 3);
      Fraction second = new Fraction(2, 5);
      Fraction third = new Fraction(3, 7);
      Fraction lcm1 = first.lcm(second);
        assertEquals(lcm1, Fraction.getReducedFraction(2, 1));
      Fraction lcm2 = lcm1.lcm(third);
        assertEquals(lcm2, Fraction.getReducedFraction(6, 1));
    }

    @Test
    void testGetReducedFraction() {
        Fraction threeFourths = new Fraction(3, 4);
        assertEquals(threeFourths, Fraction.getReducedFraction(6, 8));
        assertEquals(Fraction.ZERO, Fraction.getReducedFraction(0, -1));
        try {
            Fraction.getReducedFraction(1, 0);
            fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
            // expected
        }
        assertEquals(-1, Fraction.getReducedFraction
            (2, Integer.MIN_VALUE).getNumerator());
        assertEquals(-1, Fraction.getReducedFraction
            (1, -1).getNumerator());
    }

    @Test
    void testNormalizedEquals() {
        assertEquals(new Fraction(237, -3871), new Fraction(-51, 833));
    }

    @Test
    void testToString() {
        assertEquals("0", new Fraction(0, 3).toString());
        assertEquals("3", new Fraction(6, 2).toString());
        assertEquals("2 / 3", new Fraction(18, 27).toString());
    }

    @Test
    void testSerial() {
        Fraction[] fractions = {
            new Fraction(3, 4), Fraction.ONE, Fraction.ZERO,
            new Fraction(17), new Fraction(FastMath.PI, 1000),
            new Fraction(-5, 2)
        };
        for (Fraction fraction : fractions) {
            assertEquals(fraction, UnitTestUtils.serializeAndRecover(fraction));
        }
    }

    @Test
    void testConvergents() {
        // OEIS A002485, Numerators of convergents to Pi (https://oeis.org/A002485)
        // 0, 1, 3, 22, 333, 355, 103993, 104348, 208341, 312689, 833719, 1146408, 4272943, 5419351, 80143857, 165707065, 245850922
        // OEIS A002486, Apart from two leading terms (which are present by convention), denominators of convergents to Pi (https://oeis.org/A002486)
        // 1, 0, 1,  7, 106, 113,  33102,  33215,  66317,  99532, 265381,  364913, 1360120, 1725033, 25510582,  52746197, 78256779
        List<Fraction> convergents = Fraction.convergents(FastMath.PI, 20).collect(Collectors.toList());
        assertEquals(13, convergents.size());
        assertEquals(new Fraction(       3,        1), convergents.get( 0));
        assertEquals(new Fraction(      22,        7), convergents.get( 1));
        assertEquals(new Fraction(     333,      106), convergents.get( 2));
        assertEquals(new Fraction(     355,      113), convergents.get( 3));
        assertEquals(new Fraction(  103993,    33102), convergents.get( 4));
        assertEquals(new Fraction(  104348,    33215), convergents.get( 5));
        assertEquals(new Fraction(  208341,    66317), convergents.get( 6));
        assertEquals(new Fraction(  312689,    99532), convergents.get( 7));
        assertEquals(new Fraction(  833719,   265381), convergents.get( 8));
        assertEquals(new Fraction( 1146408,   364913), convergents.get( 9));
        assertEquals(new Fraction( 4272943,  1360120), convergents.get(10));
        assertEquals(new Fraction( 5419351,  1725033), convergents.get(11));
        assertEquals(new Fraction(80143857, 25510582), convergents.get(12));
    }

    @Test
    void testLimitedConvergents() {
        double value = FastMath.PI;
        assertEquals(new Fraction(  208341,    66317),
                Fraction.convergent(value, 7, (p, q) -> Precision.equals(p / (double) q, value, 1)).getKey());
    }

    @Test
    void testTruncatedConvergents() {
        final double value = FastMath.PI;
        assertEquals(new Fraction(   355,   113),
                Fraction.convergent(value, 20, (p, q) -> FastMath.abs(p / (double) q - value) < 1.0e-6).getKey());
        assertEquals(new Fraction(312689, 99532),
                Fraction.convergent(value, 20, (p, q) -> FastMath.abs(p / (double) q - value) < 1.0e-10).getKey());
    }

}
