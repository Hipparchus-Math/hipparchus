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
package org.hipparchus.util;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathRuntimeException;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;

import static org.junit.jupiter.api.Assertions.assertThrows;

public class BigRealTest {

    @Test
    public void testConstructor() {
        Assertions.assertEquals(1.625,
                            new BigReal(new BigDecimal("1.625")).doubleValue(),
                            1.0e-15);
        Assertions.assertEquals(-5.0,
                            new BigReal(new BigInteger("-5")).doubleValue(),
                            1.0e-15);
        Assertions.assertEquals(-5.0, new BigReal(new BigInteger("-5"),
                                              MathContext.DECIMAL64)
            .doubleValue(), 1.0e-15);
        Assertions
            .assertEquals(0.125,
                          new BigReal(new BigInteger("125"), 3).doubleValue(),
                          1.0e-15);
        Assertions.assertEquals(0.125, new BigReal(new BigInteger("125"), 3,
                                               MathContext.DECIMAL64)
            .doubleValue(), 1.0e-15);
        Assertions.assertEquals(1.625, new BigReal(new char[] {
            '1', '.', '6', '2', '5'
        }).doubleValue(), 1.0e-15);
        Assertions.assertEquals(1.625, new BigReal(new char[] {
            'A', 'A', '1', '.', '6', '2', '5', '9'
        }, 2, 5).doubleValue(), 1.0e-15);
        Assertions.assertEquals(1.625, new BigReal(new char[] {
            'A', 'A', '1', '.', '6', '2', '5', '9'
        }, 2, 5, MathContext.DECIMAL64).doubleValue(), 1.0e-15);
        Assertions.assertEquals(1.625, new BigReal(new char[] {
            '1', '.', '6', '2', '5'
        }, MathContext.DECIMAL64).doubleValue(), 1.0e-15);
        Assertions.assertEquals(1.625, new BigReal(1.625).doubleValue(), 1.0e-15);
        Assertions.assertEquals(1.625, new BigReal(1.625, MathContext.DECIMAL64)
            .doubleValue(), 1.0e-15);
        Assertions.assertEquals(-5.0, new BigReal(-5).doubleValue(), 1.0e-15);
        Assertions.assertEquals(-5.0, new BigReal(-5, MathContext.DECIMAL64)
            .doubleValue(), 1.0e-15);
        Assertions.assertEquals(-5.0, new BigReal(-5l).doubleValue(), 1.0e-15);
        Assertions.assertEquals(-5.0, new BigReal(-5l, MathContext.DECIMAL64)
            .doubleValue(), 1.0e-15);
        Assertions.assertEquals(1.625, new BigReal("1.625").doubleValue(), 1.0e-15);
        Assertions.assertEquals(1.625, new BigReal("1.625", MathContext.DECIMAL64)
            .doubleValue(), 1.0e-15);
    }

    @Test
    public void testCompareTo() {
        BigReal first = new BigReal(1.0 / 2.0);
        BigReal second = new BigReal(1.0 / 3.0);
        BigReal third = new BigReal(1.0 / 2.0);

        Assertions.assertEquals(0, first.compareTo(first));
        Assertions.assertEquals(0, first.compareTo(third));
        Assertions.assertEquals(1, first.compareTo(second));
        Assertions.assertEquals(-1, second.compareTo(first));

    }

    @Test
    public void testAdd() {
        BigReal a = new BigReal("1.2345678");
        BigReal b = new BigReal("8.7654321");
        Assertions.assertEquals(9.9999999, a.add(b).doubleValue(), 1.0e-15);
    }

    @Test
    public void testSubtract() {
        BigReal a = new BigReal("1.2345678");
        BigReal b = new BigReal("8.7654321");
        Assertions.assertEquals(-7.5308643, a.subtract(b).doubleValue(), 1.0e-15);
    }

    @Test
    public void testNegate() {
        BigReal a = new BigReal("1.2345678");
        BigReal zero = new BigReal("0.0000000");
        Assertions.assertEquals(a.negate().add(a), zero);
        Assertions.assertEquals(a.add(a.negate()), zero);
        Assertions.assertEquals(zero, zero.negate());
    }

    @Test
    public void testDivide() {
        BigReal a = new BigReal("1.0000000000");
        BigReal b = new BigReal("0.0009765625");
        Assertions.assertEquals(1024.0, a.divide(b).doubleValue(), 1.0e-15);
    }

    @Test
    public void testDivisionByZero() {
        assertThrows(MathRuntimeException.class, () -> {
            final BigReal a = BigReal.ONE;
            final BigReal b = BigReal.ZERO;
            a.divide(b);
        });
    }

    @Test
    public void testReciprocal() {
        BigReal a = new BigReal("1.2345678");
        double eps = FastMath.pow(10., -a.getScale());
        BigReal one = new BigReal("1.0000000");
        BigReal b = a.reciprocal();
        BigReal r = one.subtract(a.multiply(b));
        Assertions.assertTrue(FastMath.abs(r.doubleValue()) <= eps);
        r = one.subtract(b.multiply(a));
        Assertions.assertTrue(FastMath.abs(r.doubleValue()) <= eps);
    }

    @Test
    public void testReciprocalOfZero() {
        assertThrows(MathRuntimeException.class, () -> {
            BigReal.ZERO.reciprocal();
        });
    }

    @Test
    public void testMultiply() {
        BigReal a = new BigReal("1024.0");
        BigReal b = new BigReal("0.0009765625");
        Assertions.assertEquals(1.0, a.multiply(b).doubleValue(), 1.0e-15);
        int n = 1024;
        Assertions.assertEquals(1.0, b.multiply(n).doubleValue(), 1.0e-15);
    }

    @Test
    public void testDoubleValue() {
        Assertions.assertEquals(0.5, new BigReal(0.5).doubleValue(), 1.0e-15);
    }

    @Test
    public void testBigDecimalValue() {
        BigDecimal pi = new BigDecimal(
                                       "3.1415926535897932384626433832795028841971693993751");
        Assertions.assertEquals(pi, new BigReal(pi).bigDecimalValue());
        Assertions.assertEquals(new BigDecimal(0.5),
                            new BigReal(1.0 / 2.0).bigDecimalValue());
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    public void testEqualsAndHashCode() {
        BigReal zero = new BigReal(0.0);
        BigReal nullReal = null;
        Assertions.assertEquals(zero, zero);
        Assertions.assertNotEquals(zero, nullReal);
        Assertions.assertNotEquals(zero, Double.valueOf(0));
        BigReal zero2 = new BigReal(0.0);
        Assertions.assertEquals(zero, zero2);
        Assertions.assertEquals(zero.hashCode(), zero2.hashCode());
        BigReal one = new BigReal(1.0);
        Assertions.assertFalse((one.equals(zero) || zero.equals(one)));
        Assertions.assertEquals(BigReal.ONE, one);
    }

    @Test
    public void testSerial() {
        BigReal[] Reals = {
            new BigReal(3.0), BigReal.ONE, BigReal.ZERO, new BigReal(17),
            new BigReal(FastMath.PI), new BigReal(-2.5)
        };
        for (BigReal Real : Reals) {
            Assertions.assertEquals(Real, UnitTestUtils.serializeAndRecover(Real));
        }
    }
}
