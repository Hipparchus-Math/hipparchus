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

import org.hipparchus.CalculusFieldElementAbstractTest;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class Binary64Test extends CalculusFieldElementAbstractTest<Binary64> {
    public static final double X = 1.2345;

    public static final Binary64 PLUS_X = new Binary64(X);

    public static final Binary64 MINUS_X = new Binary64(-X);

    public static final double Y = 6.789;

    public static final Binary64 PLUS_Y = new Binary64(Y);

    public static final Binary64 MINUS_Y = new Binary64(-Y);

    public static final Binary64 PLUS_ZERO = new Binary64(0.0);

    public static final Binary64 MINUS_ZERO = new Binary64(-0.0);

    @Override
    protected Binary64 build(final double x) {
        return new Binary64(x);
    }

    @Test
    public void testAdd() {
        Binary64 expected, actual;

        expected = new Binary64(X + Y);
        actual = PLUS_X.add(PLUS_Y);
        Assertions.assertEquals(expected, actual);
        actual = PLUS_Y.add(PLUS_X);
        Assertions.assertEquals(expected, actual);

        expected = new Binary64(X + (-Y));
        actual = PLUS_X.add(MINUS_Y);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_Y.add(PLUS_X);
        Assertions.assertEquals(expected, actual);

        expected = new Binary64((-X) + (-Y));
        actual = MINUS_X.add(MINUS_Y);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_Y.add(MINUS_X);
        Assertions.assertEquals(expected, actual);

        expected = Binary64.POSITIVE_INFINITY;
        actual = PLUS_X.add(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.POSITIVE_INFINITY.add(PLUS_X);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_X.add(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.POSITIVE_INFINITY.add(MINUS_X);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.POSITIVE_INFINITY.add(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);

        expected = Binary64.NEGATIVE_INFINITY;
        actual = PLUS_X.add(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY.add(PLUS_X);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY.add(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_X.add(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY.add(MINUS_X);
        Assertions.assertEquals(expected, actual);

        expected = Binary64.NAN;
        actual = Binary64.POSITIVE_INFINITY.add(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY.add(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = PLUS_X.add(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.add(PLUS_X);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_X.add(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.add(MINUS_X);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.POSITIVE_INFINITY.add(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.add(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY.add(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.add(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.add(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
    }

    @Test
    public void testSubtract() {
        Binary64 expected, actual;

        expected = new Binary64(X - Y);
        actual = PLUS_X.subtract(PLUS_Y);
        Assertions.assertEquals(expected, actual);

        expected = new Binary64(X - (-Y));
        actual = PLUS_X.subtract(MINUS_Y);
        Assertions.assertEquals(expected, actual);

        expected = new Binary64((-X) - Y);
        actual = MINUS_X.subtract(PLUS_Y);
        Assertions.assertEquals(expected, actual);

        expected = new Binary64((-X) - (-Y));
        actual = MINUS_X.subtract(MINUS_Y);
        Assertions.assertEquals(expected, actual);

        expected = Binary64.NEGATIVE_INFINITY;
        actual = PLUS_X.subtract(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_X.subtract(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY
                .subtract(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);

        expected = Binary64.POSITIVE_INFINITY;
        actual = PLUS_X.subtract(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_X.subtract(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.POSITIVE_INFINITY
                .subtract(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);

        expected = Binary64.NAN;
        actual = Binary64.POSITIVE_INFINITY
                .subtract(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY
                .subtract(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = PLUS_X.subtract(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.subtract(PLUS_X);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_X.subtract(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.subtract(MINUS_X);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.POSITIVE_INFINITY.subtract(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.subtract(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY.subtract(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.subtract(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.subtract(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
    }

    @Test
    public void testNegate() {
        Binary64 expected, actual;

        expected = MINUS_X;
        actual = PLUS_X.negate();
        Assertions.assertEquals(expected, actual);

        expected = PLUS_X;
        actual = MINUS_X.negate();
        Assertions.assertEquals(expected, actual);

        expected = MINUS_ZERO;
        actual = PLUS_ZERO.negate();
        Assertions.assertEquals(expected, actual);

        expected = PLUS_ZERO;
        actual = MINUS_ZERO.negate();
        Assertions.assertEquals(expected, actual);

        expected = Binary64.POSITIVE_INFINITY;
        actual = Binary64.NEGATIVE_INFINITY.negate();
        Assertions.assertEquals(expected, actual);

        expected = Binary64.NEGATIVE_INFINITY;
        actual = Binary64.POSITIVE_INFINITY.negate();
        Assertions.assertEquals(expected, actual);

        expected = Binary64.NAN;
        actual = Binary64.NAN.negate();
        Assertions.assertEquals(expected, actual);
    }

    @Test
    public void testMultiply() {
        Binary64 expected, actual;

        expected = new Binary64(X * Y);
        actual = PLUS_X.multiply(PLUS_Y);
        Assertions.assertEquals(expected, actual);
        actual = PLUS_Y.multiply(PLUS_X);
        Assertions.assertEquals(expected, actual);

        expected = new Binary64(X * (-Y));
        actual = PLUS_X.multiply(MINUS_Y);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_Y.multiply(PLUS_X);
        Assertions.assertEquals(expected, actual);

        expected = new Binary64((-X) * (-Y));
        actual = MINUS_X.multiply(MINUS_Y);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_Y.multiply(MINUS_X);
        Assertions.assertEquals(expected, actual);

        expected = Binary64.POSITIVE_INFINITY;
        actual = PLUS_X.multiply(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.POSITIVE_INFINITY.multiply(PLUS_X);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_X.multiply(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY.multiply(MINUS_X);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.POSITIVE_INFINITY
                .multiply(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY
                .multiply(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);

        expected = Binary64.NEGATIVE_INFINITY;
        actual = PLUS_X.multiply(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY.multiply(PLUS_X);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_X.multiply(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.POSITIVE_INFINITY.multiply(MINUS_X);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.POSITIVE_INFINITY
                .multiply(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY
                .multiply(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);

        expected = Binary64.NAN;
        actual = PLUS_X.multiply(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.multiply(PLUS_X);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_X.multiply(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.multiply(MINUS_X);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.POSITIVE_INFINITY.multiply(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.multiply(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY.multiply(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.multiply(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.multiply(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
    }

    @Test
    public void testDivide() {
        Binary64 expected, actual;

        expected = new Binary64(X / Y);
        actual = PLUS_X.divide(PLUS_Y);
        Assertions.assertEquals(expected, actual);

        expected = new Binary64(X / (-Y));
        actual = PLUS_X.divide(MINUS_Y);
        Assertions.assertEquals(expected, actual);

        expected = new Binary64((-X) / Y);
        actual = MINUS_X.divide(PLUS_Y);
        Assertions.assertEquals(expected, actual);

        expected = new Binary64((-X) / (-Y));
        actual = MINUS_X.divide(MINUS_Y);
        Assertions.assertEquals(expected, actual);

        expected = PLUS_ZERO;
        actual = PLUS_X.divide(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_X.divide(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);

        expected = MINUS_ZERO;
        actual = MINUS_X.divide(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = PLUS_X.divide(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);

        expected = Binary64.POSITIVE_INFINITY;
        actual = Binary64.POSITIVE_INFINITY.divide(PLUS_X);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY.divide(MINUS_X);
        Assertions.assertEquals(expected, actual);
        actual = PLUS_X.divide(PLUS_ZERO);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_X.divide(MINUS_ZERO);
        Assertions.assertEquals(expected, actual);

        expected = Binary64.NEGATIVE_INFINITY;
        actual = Binary64.POSITIVE_INFINITY.divide(MINUS_X);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY.divide(PLUS_X);
        Assertions.assertEquals(expected, actual);
        actual = PLUS_X.divide(MINUS_ZERO);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_X.divide(PLUS_ZERO);
        Assertions.assertEquals(expected, actual);

        expected = Binary64.NAN;
        actual = Binary64.POSITIVE_INFINITY
                .divide(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.POSITIVE_INFINITY
                .divide(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY
                .divide(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY
                .divide(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = PLUS_X.divide(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.divide(PLUS_X);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_X.divide(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.divide(MINUS_X);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.POSITIVE_INFINITY.divide(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.divide(Binary64.POSITIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NEGATIVE_INFINITY.divide(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.divide(Binary64.NEGATIVE_INFINITY);
        Assertions.assertEquals(expected, actual);
        actual = Binary64.NAN.divide(Binary64.NAN);
        Assertions.assertEquals(expected, actual);
        actual = PLUS_ZERO.divide(PLUS_ZERO);
        Assertions.assertEquals(expected, actual);
        actual = PLUS_ZERO.divide(MINUS_ZERO);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_ZERO.divide(PLUS_ZERO);
        Assertions.assertEquals(expected, actual);
        actual = MINUS_ZERO.divide(MINUS_ZERO);
        Assertions.assertEquals(expected, actual);
    }

    @Test
    public void testReciprocal() {
        Binary64 expected, actual;

        expected = new Binary64(1.0 / X);
        actual = PLUS_X.reciprocal();
        Assertions.assertEquals(expected, actual);

        expected = new Binary64(1.0 / (-X));
        actual = MINUS_X.reciprocal();
        Assertions.assertEquals(expected, actual);

        expected = PLUS_ZERO;
        actual = Binary64.POSITIVE_INFINITY.reciprocal();
        Assertions.assertEquals(expected, actual);

        expected = MINUS_ZERO;
        actual = Binary64.NEGATIVE_INFINITY.reciprocal();
        Assertions.assertEquals(expected, actual);
    }

    @Test
    public void testIsInfinite() {
        Assertions.assertFalse(MINUS_X.isInfinite());
        Assertions.assertFalse(PLUS_X.isInfinite());
        Assertions.assertFalse(MINUS_Y.isInfinite());
        Assertions.assertFalse(PLUS_Y.isInfinite());
        Assertions.assertFalse(Binary64.NAN.isInfinite());

        Assertions.assertTrue(Binary64.NEGATIVE_INFINITY.isInfinite());
        Assertions.assertTrue(Binary64.POSITIVE_INFINITY.isInfinite());
    }

    @Test
    public void testIsNaN() {
        Assertions.assertFalse(MINUS_X.isNaN());
        Assertions.assertFalse(PLUS_X.isNaN());
        Assertions.assertFalse(MINUS_Y.isNaN());
        Assertions.assertFalse(PLUS_Y.isNaN());
        Assertions.assertFalse(Binary64.NEGATIVE_INFINITY.isNaN());
        Assertions.assertFalse(Binary64.POSITIVE_INFINITY.isNaN());

        Assertions.assertTrue(Binary64.NAN.isNaN());
    }

    @Test
    public void testToDegreesDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            for (double x = 0.1; x < 1.2; x += 0.001) {
                Binary64 value = new Binary64(x);
                Assertions.assertEquals(FastMath.toDegrees(x), value.toDegrees().getReal(), epsilon);
            }
        }
    }

    @Test
    public void testToRadiansDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            for (double x = 0.1; x < 1.2; x += 0.001) {
                Binary64 value = new Binary64(x);
                Assertions.assertEquals(FastMath.toRadians(x), value.toRadians().getReal(), epsilon);
            }
        }
    }

    @Test
    public void testDegRad() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            Binary64 value = new Binary64(x);
            Binary64 rebuilt = value.toDegrees().toRadians();
            Binary64 zero = rebuilt.subtract(value);
            Assertions.assertEquals(0, zero.getReal(), 3.0e-16);
        }
    }

    @Test
    public void testRootNegative() {
        final Binary64 neg64      = new Binary64(-64);
        final Binary64 root3Neg64 = neg64.rootN(3);
        final Binary64 root2Neg64 = neg64.rootN(2);
        Assertions.assertEquals(-4.0, root3Neg64.getReal(), 1.0e-15);
        Assertions.assertTrue(root2Neg64.isNaN());
    }

    @Test
    public void testSignedZeroEquality() {
        Assertions.assertFalse(new Binary64(1.0).isZero());
        Assertions.assertTrue(new Binary64(-0.0).isZero());
        Assertions.assertTrue(new Binary64(+0.0).isZero());
        Assertions.assertNotEquals(new Binary64(+0.0), new Binary64(-0.0));
    }

    @Test
    public void testValues() {
        Assertions.assertEquals(1,    new Binary64(1.2).byteValue());
        Assertions.assertEquals(1,    new Binary64(1.2).shortValue());
        Assertions.assertEquals(1,    new Binary64(1.2).intValue());
        Assertions.assertEquals(1l,   new Binary64(1.2).longValue());
        Assertions.assertEquals(1.2f, new Binary64(1.2).floatValue(),  0.00001f);
        Assertions.assertEquals(1.2 , new Binary64(1.2).doubleValue(), 1.0e-15);
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    public void testEquals() {
        Assertions.assertEquals(new Binary64(1.25), new Binary64(1.0).add(new Binary64(0.25)));
        Assertions.assertNotEquals(new Binary64(1.25), new Binary64(1.0).add(new Binary64(1.25)));
        Assertions.assertNotEquals("1.25", new Binary64(1.25));
    }

}
