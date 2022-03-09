/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
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

import static org.junit.Assert.assertThrows;

import java.util.List;
import java.util.stream.Collectors;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.Assert;
import org.junit.Test;


/**
 */
public class FractionTest {

    private void assertFraction(int expectedNumerator, int expectedDenominator, Fraction actual) {
        Assert.assertEquals(expectedNumerator, actual.getNumerator());
        Assert.assertEquals(expectedDenominator, actual.getDenominator());
    }

    @Test
    public void testConstructor() {
        assertFraction(0, 1, new Fraction(0, 1));
        assertFraction(0, 1, new Fraction(0, 2));
        assertFraction(0, 1, new Fraction(0, -1));
        assertFraction(1, 2, new Fraction(1, 2));
        assertFraction(1, 2, new Fraction(2, 4));
        assertFraction(-1, 2, new Fraction(-1, 2));
        assertFraction(-1, 2, new Fraction(1, -2));
        assertFraction(-1, 2, new Fraction(-2, 4));
        assertFraction(-1, 2, new Fraction(2, -4));

        // overflow
        try {
            new Fraction(Integer.MIN_VALUE, -1);
            Assert.fail();
        } catch (MathRuntimeException ex) {
            // success
        }
        try {
            new Fraction(1, Integer.MIN_VALUE);
            Assert.fail();
        } catch (MathRuntimeException ex) {
            // success
        }

        assertFraction(0, 1, new Fraction(0.00000000000001));
        assertFraction(2, 5, new Fraction(0.40000000000001));
        assertFraction(15, 1, new Fraction(15.0000000000001));
    }

    @Test
    public void testIsInteger() {
        Assert.assertTrue(new Fraction(12, 12).isInteger());
        Assert.assertTrue(new Fraction(14, 7).isInteger());
        Assert.assertFalse(new Fraction(12, 11).isInteger());
    }

    @Test(expected=MathIllegalStateException.class)
    public void testGoldenRatio() {
        // the golden ratio is notoriously a difficult number for continuous fraction
        new Fraction((1 + FastMath.sqrt(5)) / 2, 1.0e-12, 25);
    }

    // MATH-179
    @Test
    public void testDoubleConstructor() {
        assertFraction(1, 2, new Fraction((double)1 / (double)2));
        assertFraction(1, 3, new Fraction((double)1 / (double)3));
        assertFraction(2, 3, new Fraction((double)2 / (double)3));
        assertFraction(1, 4, new Fraction((double)1 / (double)4));
        assertFraction(3, 4, new Fraction((double)3 / (double)4));
        assertFraction(1, 5, new Fraction((double)1 / (double)5));
        assertFraction(2, 5, new Fraction((double)2 / (double)5));
        assertFraction(3, 5, new Fraction((double)3 / (double)5));
        assertFraction(4, 5, new Fraction((double)4 / (double)5));
        assertFraction(1, 6, new Fraction((double)1 / (double)6));
        assertFraction(5, 6, new Fraction((double)5 / (double)6));
        assertFraction(1, 7, new Fraction((double)1 / (double)7));
        assertFraction(2, 7, new Fraction((double)2 / (double)7));
        assertFraction(3, 7, new Fraction((double)3 / (double)7));
        assertFraction(4, 7, new Fraction((double)4 / (double)7));
        assertFraction(5, 7, new Fraction((double)5 / (double)7));
        assertFraction(6, 7, new Fraction((double)6 / (double)7));
        assertFraction(1, 8, new Fraction((double)1 / (double)8));
        assertFraction(3, 8, new Fraction((double)3 / (double)8));
        assertFraction(5, 8, new Fraction((double)5 / (double)8));
        assertFraction(7, 8, new Fraction((double)7 / (double)8));
        assertFraction(1, 9, new Fraction((double)1 / (double)9));
        assertFraction(2, 9, new Fraction((double)2 / (double)9));
        assertFraction(4, 9, new Fraction((double)4 / (double)9));
        assertFraction(5, 9, new Fraction((double)5 / (double)9));
        assertFraction(7, 9, new Fraction((double)7 / (double)9));
        assertFraction(8, 9, new Fraction((double)8 / (double)9));
        assertFraction(1, 10, new Fraction((double)1 / (double)10));
        assertFraction(3, 10, new Fraction((double)3 / (double)10));
        assertFraction(7, 10, new Fraction((double)7 / (double)10));
        assertFraction(9, 10, new Fraction((double)9 / (double)10));
        assertFraction(1, 11, new Fraction((double)1 / (double)11));
        assertFraction(2, 11, new Fraction((double)2 / (double)11));
        assertFraction(3, 11, new Fraction((double)3 / (double)11));
        assertFraction(4, 11, new Fraction((double)4 / (double)11));
        assertFraction(5, 11, new Fraction((double)5 / (double)11));
        assertFraction(6, 11, new Fraction((double)6 / (double)11));
        assertFraction(7, 11, new Fraction((double)7 / (double)11));
        assertFraction(8, 11, new Fraction((double)8 / (double)11));
        assertFraction(9, 11, new Fraction((double)9 / (double)11));
        assertFraction(10, 11, new Fraction((double)10 / (double)11));
    }

    // MATH-181
    @Test
    public void testDigitLimitConstructor() {
        assertFraction(2, 5, new Fraction(0.4,   9));
        assertFraction(2, 5, new Fraction(0.4,  99));
        assertFraction(2, 5, new Fraction(0.4, 999));

        assertFraction(3, 5,      new Fraction(0.6152,    9));
        assertFraction(8, 13,     new Fraction(0.6152,   99));
        assertFraction(510, 829,  new Fraction(0.6152,  999));
        assertFraction(769, 1250, new Fraction(0.6152, 9999));

        // MATH-996
        assertFraction(1, 2, new Fraction(0.5000000001, 10));
    }

    @Test
    public void testIntegerOverflow() {
        checkIntegerOverflow(0.75000000001455192);
        checkIntegerOverflow(1.0e10);
        checkIntegerOverflow(-1.0e10);
    }

    @Test
    public void testSignum() {
        Assert.assertEquals(-1, new Fraction(4, -5).signum());
        Assert.assertEquals(-1, new Fraction(-4, 5).signum());
        Assert.assertEquals( 0, new Fraction(0).signum());
        Assert.assertEquals(+1, new Fraction(-4, -5).signum());
        Assert.assertEquals(+1, new Fraction(4, 5).signum());
    }

    private void checkIntegerOverflow(double a) {
        assertThrows(MathIllegalStateException.class, () -> new Fraction(a, 1.0e-12, 1000));
    }

    @Test
    public void testEpsilonLimitConstructor() {
        assertFraction(2, 5, new Fraction(0.4, 1.0e-5, 100));

        assertFraction(3, 5,      new Fraction(0.6152, 0.02, 100));
        assertFraction(8, 13,     new Fraction(0.6152, 1.0e-3, 100));
        assertFraction(251, 408,  new Fraction(0.6152, 1.0e-4, 100));
        assertFraction(251, 408,  new Fraction(0.6152, 1.0e-5, 100));
        assertFraction(510, 829,  new Fraction(0.6152, 1.0e-6, 100));
        assertFraction(769, 1250, new Fraction(0.6152, 1.0e-7, 100));
    }

    @Test
    public void testCompareTo() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(1, 3);
        Fraction third = new Fraction(1, 2);

        Assert.assertEquals(0, first.compareTo(first));
        Assert.assertEquals(0, first.compareTo(third));
        Assert.assertEquals(1, first.compareTo(second));
        Assert.assertEquals(-1, second.compareTo(first));

        // these two values are different approximations of PI
        // the first  one is approximately PI - 3.07e-18
        // the second one is approximately PI + 1.936e-17
        Fraction pi1 = new Fraction(1068966896, 340262731);
        Fraction pi2 = new Fraction( 411557987, 131002976);
        Assert.assertEquals(-1, pi1.compareTo(pi2));
        Assert.assertEquals( 1, pi2.compareTo(pi1));
        Assert.assertEquals(0.0, pi1.doubleValue() - pi2.doubleValue(), 1.0e-20);
    }

    @Test
    public void testDoubleValue() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(1, 3);

        Assert.assertEquals(0.5, first.doubleValue(), 0.0);
        Assert.assertEquals(1.0 / 3.0, second.doubleValue(), 0.0);
    }

    @Test
    public void testFloatValue() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(1, 3);

        Assert.assertEquals(0.5f, first.floatValue(), 0.0f);
        Assert.assertEquals((float)(1.0 / 3.0), second.floatValue(), 0.0f);
    }

    @Test
    public void testIntValue() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(3, 2);

        Assert.assertEquals(0, first.intValue());
        Assert.assertEquals(1, second.intValue());
    }

    @Test
    public void testLongValue() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(3, 2);

        Assert.assertEquals(0L, first.longValue());
        Assert.assertEquals(1L, second.longValue());
    }

    @Test
    public void testConstructorDouble() {
        assertFraction(1, 2, new Fraction(0.5));
        assertFraction(1, 3, new Fraction(1.0 / 3.0));
        assertFraction(17, 100, new Fraction(17.0 / 100.0));
        assertFraction(317, 100, new Fraction(317.0 / 100.0));
        assertFraction(-1, 2, new Fraction(-0.5));
        assertFraction(-1, 3, new Fraction(-1.0 / 3.0));
        assertFraction(-17, 100, new Fraction(17.0 / -100.0));
        assertFraction(-317, 100, new Fraction(-317.0 / 100.0));
    }

    @Test
    public void testAbs() {
        Fraction a = new Fraction(10, 21);
        Fraction b = new Fraction(-10, 21);
        Fraction c = new Fraction(10, -21);

        assertFraction(10, 21, a.abs());
        assertFraction(10, 21, b.abs());
        assertFraction(10, 21, c.abs());
    }

    @Test
    public void testPercentage() {
        Assert.assertEquals(50.0, new Fraction(1, 2).percentageValue(), 1.0e-15);
    }

    @Test
    public void testMath835() {
        final int numer = Integer.MAX_VALUE / 99;
        final int denom = 1;
        final double percentage = 100 * ((double) numer) / denom;
        final Fraction frac = new Fraction(numer, denom);
        // With the implementation that preceded the fix suggested in MATH-835,
        // this test was failing, due to overflow.
        Assert.assertEquals(percentage, frac.percentageValue(), Math.ulp(percentage));
    }

    @Test
    public void testMath1261() {
        final Fraction a = new Fraction(Integer.MAX_VALUE, 2);
        final Fraction b = a.multiply(2);
        Assert.assertTrue(b.equals(new Fraction(Integer.MAX_VALUE)));

        final Fraction c = new Fraction(2, Integer.MAX_VALUE);
        final Fraction d = c.divide(2);
        Assert.assertTrue(d.equals(new Fraction(1, Integer.MAX_VALUE)));
    }

    @Test
    public void testReciprocal() {
        Fraction f = null;

        f = new Fraction(50, 75);
        f = f.reciprocal();
        Assert.assertEquals(3, f.getNumerator());
        Assert.assertEquals(2, f.getDenominator());

        f = new Fraction(4, 3);
        f = f.reciprocal();
        Assert.assertEquals(3, f.getNumerator());
        Assert.assertEquals(4, f.getDenominator());

        f = new Fraction(-15, 47);
        f = f.reciprocal();
        Assert.assertEquals(-47, f.getNumerator());
        Assert.assertEquals(15, f.getDenominator());

        f = new Fraction(0, 3);
        try {
            f = f.reciprocal();
            Assert.fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        // large values
        f = new Fraction(Integer.MAX_VALUE, 1);
        f = f.reciprocal();
        Assert.assertEquals(1, f.getNumerator());
        Assert.assertEquals(Integer.MAX_VALUE, f.getDenominator());
    }

    @Test
    public void testNegate() {
        Fraction f = null;

        f = new Fraction(50, 75);
        f = f.negate();
        Assert.assertEquals(-2, f.getNumerator());
        Assert.assertEquals(3, f.getDenominator());

        f = new Fraction(-50, 75);
        f = f.negate();
        Assert.assertEquals(2, f.getNumerator());
        Assert.assertEquals(3, f.getDenominator());

        // large values
        f = new Fraction(Integer.MAX_VALUE-1, Integer.MAX_VALUE);
        f = f.negate();
        Assert.assertEquals(Integer.MIN_VALUE+2, f.getNumerator());
        Assert.assertEquals(Integer.MAX_VALUE, f.getDenominator());

        f = new Fraction(Integer.MIN_VALUE, 1);
        try {
            f = f.negate();
            Assert.fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}
    }

    @Test
    public void testAdd() {
        Fraction a = new Fraction(1, 2);
        Fraction b = new Fraction(2, 3);

        assertFraction(1, 1, a.add(a));
        assertFraction(7, 6, a.add(b));
        assertFraction(7, 6, b.add(a));
        assertFraction(4, 3, b.add(b));

        Fraction f1 = new Fraction(Integer.MAX_VALUE - 1, 1);
        Fraction f2 = Fraction.ONE;
        Fraction f = f1.add(f2);
        Assert.assertEquals(Integer.MAX_VALUE, f.getNumerator());
        Assert.assertEquals(1, f.getDenominator());
        f = f1.add(1);
        Assert.assertEquals(Integer.MAX_VALUE, f.getNumerator());
        Assert.assertEquals(1, f.getDenominator());

        f1 = new Fraction(-1, 13*13*2*2);
        f2 = new Fraction(-2, 13*17*2);
        f = f1.add(f2);
        Assert.assertEquals(13*13*17*2*2, f.getDenominator());
        Assert.assertEquals(-17 - 2*13*2, f.getNumerator());

        try {
            f.add(null);
            Assert.fail("expecting NullArgumentException");
        } catch (NullArgumentException ex) {}

        // if this fraction is added naively, it will overflow.
        // check that it doesn't.
        f1 = new Fraction(1,32768*3);
        f2 = new Fraction(1,59049);
        f = f1.add(f2);
        Assert.assertEquals(52451, f.getNumerator());
        Assert.assertEquals(1934917632, f.getDenominator());

        f1 = new Fraction(Integer.MIN_VALUE, 3);
        f2 = new Fraction(1,3);
        f = f1.add(f2);
        Assert.assertEquals(Integer.MIN_VALUE+1, f.getNumerator());
        Assert.assertEquals(3, f.getDenominator());

        f1 = new Fraction(Integer.MAX_VALUE - 1, 1);
        f2 = Fraction.ONE;
        f = f1.add(f2);
        Assert.assertEquals(Integer.MAX_VALUE, f.getNumerator());
        Assert.assertEquals(1, f.getDenominator());

        try {
            f = f.add(Fraction.ONE); // should overflow
            Assert.fail("expecting MathRuntimeException but got: " + f.toString());
        } catch (MathRuntimeException ex) {}

        // denominator should not be a multiple of 2 or 3 to trigger overflow
        f1 = new Fraction(Integer.MIN_VALUE, 5);
        f2 = new Fraction(-1,5);
        try {
            f = f1.add(f2); // should overflow
            Assert.fail("expecting MathRuntimeException but got: " + f.toString());
        } catch (MathRuntimeException ex) {}

        try {
            f= new Fraction(-Integer.MAX_VALUE, 1);
            f = f.add(f);
            Assert.fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        try {
            f= new Fraction(-Integer.MAX_VALUE, 1);
            f = f.add(f);
            Assert.fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        f1 = new Fraction(3,327680);
        f2 = new Fraction(2,59049);
        try {
            f = f1.add(f2); // should overflow
            Assert.fail("expecting MathRuntimeException but got: " + f.toString());
        } catch (MathRuntimeException ex) {}
    }

    @Test
    public void testDivide() {
        Fraction a = new Fraction(1, 2);
        Fraction b = new Fraction(2, 3);

        assertFraction(1, 1, a.divide(a));
        assertFraction(3, 4, a.divide(b));
        assertFraction(4, 3, b.divide(a));
        assertFraction(1, 1, b.divide(b));

        Fraction f1 = new Fraction(3, 5);
        Fraction f2 = Fraction.ZERO;
        try {
            f1.divide(f2);
            Assert.fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        f1 = new Fraction(0, 5);
        f2 = new Fraction(2, 7);
        Fraction f = f1.divide(f2);
        Assert.assertSame(Fraction.ZERO, f);

        f1 = new Fraction(2, 7);
        f2 = Fraction.ONE;
        f = f1.divide(f2);
        Assert.assertEquals(2, f.getNumerator());
        Assert.assertEquals(7, f.getDenominator());

        f1 = new Fraction(1, Integer.MAX_VALUE);
        f = f1.divide(f1);
        Assert.assertEquals(1, f.getNumerator());
        Assert.assertEquals(1, f.getDenominator());

        f1 = new Fraction(Integer.MIN_VALUE, Integer.MAX_VALUE);
        f2 = new Fraction(1, Integer.MAX_VALUE);
        f = f1.divide(f2);
        Assert.assertEquals(Integer.MIN_VALUE, f.getNumerator());
        Assert.assertEquals(1, f.getDenominator());

        try {
            f.divide(null);
            Assert.fail("NullArgumentException");
        } catch (NullArgumentException ex) {}

        try {
            f1 = new Fraction(1, Integer.MAX_VALUE);
            f = f1.divide(f1.reciprocal());  // should overflow
            Assert.fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}
        try {
            f1 = new Fraction(1, -Integer.MAX_VALUE);
            f = f1.divide(f1.reciprocal());  // should overflow
            Assert.fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        f1 = new Fraction(6, 35);
        f  = f1.divide(15);
        Assert.assertEquals(2, f.getNumerator());
        Assert.assertEquals(175, f.getDenominator());

    }

    @Test
    public void testMultiply() {
        Fraction a = new Fraction(1, 2);
        Fraction b = new Fraction(2, 3);

        assertFraction(1, 4, a.multiply(a));
        assertFraction(1, 3, a.multiply(b));
        assertFraction(1, 3, b.multiply(a));
        assertFraction(4, 9, b.multiply(b));

        Fraction f1 = new Fraction(Integer.MAX_VALUE, 1);
        Fraction f2 = new Fraction(Integer.MIN_VALUE, Integer.MAX_VALUE);
        Fraction f = f1.multiply(f2);
        Assert.assertEquals(Integer.MIN_VALUE, f.getNumerator());
        Assert.assertEquals(1, f.getDenominator());

        try {
            f.multiply(null);
            Assert.fail("expecting NullArgumentException");
        } catch (NullArgumentException ex) {}

        f1 = new Fraction(6, 35);
        f  = f1.multiply(15);
        Assert.assertEquals(18, f.getNumerator());
        Assert.assertEquals(7, f.getDenominator());
    }

    @Test
    public void testSubtract() {
        Fraction a = new Fraction(1, 2);
        Fraction b = new Fraction(2, 3);

        assertFraction(0, 1, a.subtract(a));
        assertFraction(-1, 6, a.subtract(b));
        assertFraction(1, 6, b.subtract(a));
        assertFraction(0, 1, b.subtract(b));

        Fraction f = new Fraction(1,1);
        try {
            f.subtract(null);
            Assert.fail("expecting NullArgumentException");
        } catch (NullArgumentException ex) {}

        // if this fraction is subtracted naively, it will overflow.
        // check that it doesn't.
        Fraction f1 = new Fraction(1,32768*3);
        Fraction f2 = new Fraction(1,59049);
        f = f1.subtract(f2);
        Assert.assertEquals(-13085, f.getNumerator());
        Assert.assertEquals(1934917632, f.getDenominator());

        f1 = new Fraction(Integer.MIN_VALUE, 3);
        f2 = new Fraction(1,3).negate();
        f = f1.subtract(f2);
        Assert.assertEquals(Integer.MIN_VALUE+1, f.getNumerator());
        Assert.assertEquals(3, f.getDenominator());

        f1 = new Fraction(Integer.MAX_VALUE, 1);
        f2 = Fraction.ONE;
        f = f1.subtract(f2);
        Assert.assertEquals(Integer.MAX_VALUE-1, f.getNumerator());
        Assert.assertEquals(1, f.getDenominator());
        f = f1.subtract(1);
        Assert.assertEquals(Integer.MAX_VALUE-1, f.getNumerator());
        Assert.assertEquals(1, f.getDenominator());

        try {
            f1 = new Fraction(1, Integer.MAX_VALUE);
            f2 = new Fraction(1, Integer.MAX_VALUE - 1);
            f = f1.subtract(f2);
            Assert.fail("expecting MathRuntimeException");  //should overflow
        } catch (MathRuntimeException ex) {}

        // denominator should not be a multiple of 2 or 3 to trigger overflow
        f1 = new Fraction(Integer.MIN_VALUE, 5);
        f2 = new Fraction(1,5);
        try {
            f = f1.subtract(f2); // should overflow
            Assert.fail("expecting MathRuntimeException but got: " + f.toString());
        } catch (MathRuntimeException ex) {}

        try {
            f= new Fraction(Integer.MIN_VALUE, 1);
            f = f.subtract(Fraction.ONE);
            Assert.fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        try {
            f= new Fraction(Integer.MAX_VALUE, 1);
            f = f.subtract(Fraction.ONE.negate());
            Assert.fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {}

        f1 = new Fraction(3,327680);
        f2 = new Fraction(2,59049);
        try {
            f = f1.subtract(f2); // should overflow
            Assert.fail("expecting MathRuntimeException but got: " + f.toString());
        } catch (MathRuntimeException ex) {}
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    public void testEqualsAndHashCode() {
        Fraction zero  = new Fraction(0,1);
        Fraction nullFraction = null;
        Assert.assertTrue( zero.equals(zero));
        Assert.assertFalse(zero.equals(nullFraction));
        Assert.assertFalse(zero.equals(Double.valueOf(0)));
        Fraction zero2 = new Fraction(0,2);
        Assert.assertTrue(zero.equals(zero2));
        Assert.assertEquals(zero.hashCode(), zero2.hashCode());
        Fraction one = new Fraction(1,1);
        Assert.assertFalse((one.equals(zero) ||zero.equals(one)));
    }

    @Test
    public void testGetReducedFraction() {
        Fraction threeFourths = new Fraction(3, 4);
        Assert.assertTrue(threeFourths.equals(Fraction.getReducedFraction(6, 8)));
        Assert.assertTrue(Fraction.ZERO.equals(Fraction.getReducedFraction(0, -1)));
        try {
            Fraction.getReducedFraction(1, 0);
            Assert.fail("expecting MathRuntimeException");
        } catch (MathRuntimeException ex) {
            // expected
        }
        Assert.assertEquals(Fraction.getReducedFraction
                (2, Integer.MIN_VALUE).getNumerator(),-1);
        Assert.assertEquals(Fraction.getReducedFraction
                (1, -1).getNumerator(), -1);
    }

    @Test
    public void testNormalizedEquals() {
        Assert.assertEquals(new Fraction(237, -3871), new Fraction(-51, 833));
    }

    @Test
    public void testToString() {
        Assert.assertEquals("0", new Fraction(0, 3).toString());
        Assert.assertEquals("3", new Fraction(6, 2).toString());
        Assert.assertEquals("2 / 3", new Fraction(18, 27).toString());
    }

    @Test
    public void testSerial() {
        Fraction[] fractions = {
            new Fraction(3, 4), Fraction.ONE, Fraction.ZERO,
            new Fraction(17), new Fraction(FastMath.PI, 1000),
            new Fraction(-5, 2)
        };
        for (Fraction fraction : fractions) {
            Assert.assertEquals(fraction, UnitTestUtils.serializeAndRecover(fraction));
        }
    }

    @Test
    public void testConvergents() {
        // OEIS A002485, Numerators of convergents to Pi (https://oeis.org/A002485)
        // 0, 1, 3, 22, 333, 355, 103993, 104348, 208341, 312689, 833719, 1146408, 4272943, 5419351, 80143857, 165707065, 245850922
        // OEIS A002486, Apart from two leading terms (which are present by convention), denominators of convergents to Pi (https://oeis.org/A002486)
        // 1, 0, 1,  7, 106, 113,  33102,  33215,  66317,  99532, 265381,  364913, 1360120, 1725033, 25510582,  52746197, 78256779
        List<Fraction> convergents = Fraction.convergents(FastMath.PI, 20).collect(Collectors.toList());
        Assert.assertEquals(13, convergents.size());
        Assert.assertEquals(new Fraction(       3,        1), convergents.get( 0));
        Assert.assertEquals(new Fraction(      22,        7), convergents.get( 1));
        Assert.assertEquals(new Fraction(     333,      106), convergents.get( 2));
        Assert.assertEquals(new Fraction(     355,      113), convergents.get( 3));
        Assert.assertEquals(new Fraction(  103993,    33102), convergents.get( 4));
        Assert.assertEquals(new Fraction(  104348,    33215), convergents.get( 5));
        Assert.assertEquals(new Fraction(  208341,    66317), convergents.get( 6));
        Assert.assertEquals(new Fraction(  312689,    99532), convergents.get( 7));
        Assert.assertEquals(new Fraction(  833719,   265381), convergents.get( 8));
        Assert.assertEquals(new Fraction( 1146408,   364913), convergents.get( 9));
        Assert.assertEquals(new Fraction( 4272943,  1360120), convergents.get(10));
        Assert.assertEquals(new Fraction( 5419351,  1725033), convergents.get(11));
        Assert.assertEquals(new Fraction(80143857, 25510582), convergents.get(12));
    }

    @Test
    public void testLimitedConvergents() {
        double value = FastMath.PI;
        Assert.assertEquals(new Fraction(  208341,    66317),
                Fraction.convergent(value, 7, (p, q) -> Precision.equals(p / (double) q, value, 1)).getKey());
    }

    @Test
    public void testTruncatedConvergents() {
        final double value = FastMath.PI;
        Assert.assertEquals(new Fraction(   355,   113),
                Fraction.convergent(value, 20, (p, q) -> FastMath.abs(p / (double) q - value) < 1.0e-6).getKey());
        Assert.assertEquals(new Fraction(312689, 99532),
                Fraction.convergent(value, 20, (p, q) -> FastMath.abs(p / (double) q - value) < 1.0e-10).getKey());
    }

}
