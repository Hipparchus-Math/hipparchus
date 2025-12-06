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
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class BigRealTest {

    @Test
    void testConstructor() {
        assertEquals(1.625,
                            new BigReal(new BigDecimal("1.625")).doubleValue(),
                            1.0e-15);
        assertEquals(-5.0,
                            new BigReal(new BigInteger("-5")).doubleValue(),
                            1.0e-15);
        assertEquals(-5.0, new BigReal(new BigInteger("-5"),
                                              MathContext.DECIMAL64)
            .doubleValue(), 1.0e-15);
        assertEquals(0.125,
                          new BigReal(new BigInteger("125"), 3).doubleValue(),
                          1.0e-15);
        assertEquals(0.125, new BigReal(new BigInteger("125"), 3,
                                               MathContext.DECIMAL64)
            .doubleValue(), 1.0e-15);
        assertEquals(1.625, new BigReal(new char[] {
            '1', '.', '6', '2', '5'
        }).doubleValue(), 1.0e-15);
        assertEquals(1.625, new BigReal(new char[] {
            'A', 'A', '1', '.', '6', '2', '5', '9'
        }, 2, 5).doubleValue(), 1.0e-15);
        assertEquals(1.625, new BigReal(new char[] {
            'A', 'A', '1', '.', '6', '2', '5', '9'
        }, 2, 5, MathContext.DECIMAL64).doubleValue(), 1.0e-15);
        assertEquals(1.625, new BigReal(new char[] {
            '1', '.', '6', '2', '5'
        }, MathContext.DECIMAL64).doubleValue(), 1.0e-15);
        assertEquals(1.625, new BigReal(1.625).doubleValue(), 1.0e-15);
        assertEquals(1.625, new BigReal(1.625, MathContext.DECIMAL64)
            .doubleValue(), 1.0e-15);
        assertEquals(-5.0, new BigReal(-5).doubleValue(), 1.0e-15);
        assertEquals(-5.0, new BigReal(-5, MathContext.DECIMAL64)
            .doubleValue(), 1.0e-15);
        assertEquals(-5.0, new BigReal(-5l).doubleValue(), 1.0e-15);
        assertEquals(-5.0, new BigReal(-5l, MathContext.DECIMAL64)
            .doubleValue(), 1.0e-15);
        assertEquals(1.625, new BigReal("1.625").doubleValue(), 1.0e-15);
        assertEquals(1.625, new BigReal("1.625", MathContext.DECIMAL64)
            .doubleValue(), 1.0e-15);
    }

    @Test
    void testCompareTo() {
        BigReal first = new BigReal(1.0 / 2.0);
        BigReal second = new BigReal(1.0 / 3.0);
        BigReal third = new BigReal(1.0 / 2.0);

        assertEquals(0, first.compareTo(first));
        assertEquals(0, first.compareTo(third));
        assertEquals(1, first.compareTo(second));
        assertEquals(-1, second.compareTo(first));

    }

    @Test
    void testAdd() {
        BigReal a = new BigReal("1.2345678");
        BigReal b = new BigReal("8.7654321");
        assertEquals(9.9999999, a.add(b).doubleValue(), 1.0e-15);
    }

    @Test
    void testSubtract() {
        BigReal a = new BigReal("1.2345678");
        BigReal b = new BigReal("8.7654321");
        assertEquals(-7.5308643, a.subtract(b).doubleValue(), 1.0e-15);
    }

    @Test
    void testNegate() {
        BigReal a = new BigReal("1.2345678");
        BigReal zero = new BigReal("0.0000000");
        assertEquals(a.negate().add(a), zero);
        assertEquals(a.add(a.negate()), zero);
        assertEquals(zero, zero.negate());
    }

    @Test
    void testDivide() {
        BigReal a = new BigReal("1.0000000000");
        BigReal b = new BigReal("0.0009765625");
        assertEquals(1024.0, a.divide(b).doubleValue(), 1.0e-15);
    }

    @Test
    void testDivisionByZero() {
        assertThrows(MathRuntimeException.class, () -> {
            final BigReal a = BigReal.ONE;
            final BigReal b = BigReal.ZERO;
            a.divide(b);
        });
    }

    @Test
    void testReciprocal() {
        BigReal a = new BigReal("1.2345678");
        double eps = FastMath.pow(10., -a.getScale());
        BigReal one = new BigReal("1.0000000");
        BigReal b = a.reciprocal();
        BigReal r = one.subtract(a.multiply(b));
        assertTrue(FastMath.abs(r.doubleValue()) <= eps);
        r = one.subtract(b.multiply(a));
        assertTrue(FastMath.abs(r.doubleValue()) <= eps);
    }

    @Test
    void testReciprocalOfZero() {
        assertThrows(MathRuntimeException.class, BigReal.ZERO::reciprocal);
    }

    @Test
    void testMultiply() {
        BigReal a = new BigReal("1024.0");
        BigReal b = new BigReal("0.0009765625");
        assertEquals(1.0, a.multiply(b).doubleValue(), 1.0e-15);
        int n = 1024;
        assertEquals(1.0, b.multiply(n).doubleValue(), 1.0e-15);
    }

    @Test
    void testDoubleValue() {
        assertEquals(0.5, new BigReal(0.5).doubleValue(), 1.0e-15);
    }

    @Test
    void testBigDecimalValue() {
        BigDecimal pi = new BigDecimal(
                                       "3.1415926535897932384626433832795028841971693993751");
        assertEquals(pi, new BigReal(pi).bigDecimalValue());
        assertEquals(new BigDecimal(0.5),
                            new BigReal(1.0 / 2.0).bigDecimalValue());
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testEqualsAndHashCode() {
        BigReal zero = new BigReal(0.0);
        BigReal nullReal = null;
        assertEquals(zero, zero);
        assertNotEquals(zero, nullReal);
        assertNotEquals(zero, Double.valueOf(0));
        BigReal zero2 = new BigReal(0.0);
        assertEquals(zero, zero2);
        assertEquals(zero.hashCode(), zero2.hashCode());
        BigReal one = new BigReal(1.0);
        assertFalse((one.equals(zero) || zero.equals(one)));
        assertEquals(BigReal.ONE, one);
    }

    @Test
    void testSerial() {
        BigReal[] Reals = {
            new BigReal(3.0), BigReal.ONE, BigReal.ZERO, new BigReal(17),
            new BigReal(FastMath.PI), new BigReal(-2.5)
        };
        for (BigReal Real : Reals) {
            assertEquals(Real, UnitTestUtils.serializeAndRecover(Real));
        }
    }
}
