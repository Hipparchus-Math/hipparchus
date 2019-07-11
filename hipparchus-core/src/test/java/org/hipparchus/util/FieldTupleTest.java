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
package org.hipparchus.util;

import org.hipparchus.ExtendedFieldElementAbstractTest;
import org.junit.Assert;
import org.junit.Test;

public class FieldTupleTest extends ExtendedFieldElementAbstractTest<FieldTuple<Decimal64>> {
    public static final Decimal64 X = new Decimal64(1.2345);

    public static final FieldTuple<Decimal64> PLUS_X = new FieldTuple<Decimal64>(X, X);

    public static final FieldTuple<Decimal64> MINUS_X = new FieldTuple<Decimal64>(X.negate(), X.negate());

    public static final Decimal64 Y = new Decimal64(6.789);

    public static final FieldTuple<Decimal64> PLUS_Y = new FieldTuple<Decimal64>(Y, Y);

    public static final FieldTuple<Decimal64> MINUS_Y = new FieldTuple<Decimal64>(Y.negate(), Y.negate());

    public static final FieldTuple<Decimal64> PLUS_ZERO = new FieldTuple<Decimal64>(new Decimal64(0.0), new Decimal64(0.0));

    public static final FieldTuple<Decimal64> MINUS_ZERO = new FieldTuple<Decimal64>(new Decimal64(-0.0), new Decimal64(-0.0));

    @Override
    protected FieldTuple<Decimal64> build(final double x) {
        return new FieldTuple<Decimal64>(new Decimal64(x), new Decimal64(x));
    }

    @Test
    public void testComponents() {
        Assert.assertEquals(2, PLUS_ZERO.getDimension());
        final FieldTuple<Decimal64> oneToFive = new FieldTuple<Decimal64>(new Decimal64(1),
                                                                          new Decimal64(2),
                                                                          new Decimal64(3),
                                                                          new Decimal64(4),
                                                                          new Decimal64(5));
        Assert.assertEquals(5, oneToFive.getDimension());
        Assert.assertArrayEquals(new Decimal64[] {
            new Decimal64(1), new Decimal64(2), new Decimal64(3), new Decimal64(4), new Decimal64(5) },
                                 oneToFive.getComponents());
        Assert.assertEquals(1, oneToFive.getComponent(0).getReal(), 1.0e-10);
        Assert.assertEquals(2, oneToFive.getComponent(1).getReal(), 1.0e-10);
        Assert.assertEquals(3, oneToFive.getComponent(2).getReal(), 1.0e-10);
        Assert.assertEquals(4, oneToFive.getComponent(3).getReal(), 1.0e-10);
        Assert.assertEquals(5, oneToFive.getComponent(4).getReal(), 1.0e-10);
    }

    @Test
    public void testEquals() {
        Assert.assertNotEquals(PLUS_ZERO, null);
        Assert.assertEquals(PLUS_ZERO, PLUS_ZERO);
        Assert.assertEquals(PLUS_X,    PLUS_X);
        Assert.assertEquals(PLUS_Y,    PLUS_Y);
        Assert.assertEquals(MINUS_X,   MINUS_X);
        Assert.assertEquals(MINUS_Y,   MINUS_Y);
        Assert.assertNotEquals(PLUS_X,  new FieldTuple<>(new Decimal64(1),
                                                         new Decimal64(2),
                                                         new Decimal64(3),
                                                         new Decimal64(4),
                                                         new Decimal64(5)));
        Assert.assertNotEquals(PLUS_X,  new FieldTuple<>(PLUS_X.getComponent(0), new Decimal64(999.999)));
        Assert.assertNotEquals(PLUS_ZERO.getField(), null);
        Assert.assertNotEquals(PLUS_X.getField(),  new FieldTuple<>(new Decimal64(1),
                        new Decimal64(2),
                        new Decimal64(3),
                        new Decimal64(4),
                        new Decimal64(5)).getField());
        Assert.assertEquals(PLUS_ZERO.getField(), MINUS_Y.getField());
    }

    @Test
    public void testHashcode() {
        Assert.assertEquals(-1264241693, PLUS_ZERO.getField().hashCode());
        Assert.assertEquals(1492525478, PLUS_ZERO.hashCode());
        Assert.assertEquals(1492525478, MINUS_ZERO.hashCode());
        Assert.assertEquals(211878758, PLUS_X.hashCode());
        Assert.assertEquals(211878758, MINUS_X.hashCode());
        Assert.assertEquals(1460019622, new FieldTuple<>(new Decimal64(1), new Decimal64(2)).hashCode());
        Assert.assertEquals(1491476902, new FieldTuple<>(new Decimal64(2), new Decimal64(1)).hashCode());
    }

    @Test
    public void testAdd() {
        FieldTuple<Decimal64> expected, actual;

        expected = new FieldTuple<Decimal64>(X.add(Y), X.add(Y));
        actual = PLUS_X.add(PLUS_Y);
        Assert.assertEquals(expected, actual);
        actual = PLUS_Y.add(PLUS_X);
        Assert.assertEquals(expected, actual);

        expected = new FieldTuple<Decimal64>(X.add(Y.negate()), X.add(Y.negate()));
        actual = PLUS_X.add(MINUS_Y);
        Assert.assertEquals(expected, actual);
        actual = MINUS_Y.add(PLUS_X);
        Assert.assertEquals(expected, actual);

        expected = new FieldTuple<Decimal64>(X.negate().add(Y.negate()), X.negate().add(Y.negate()));
        actual = MINUS_X.add(MINUS_Y);
        Assert.assertEquals(expected, actual);
        actual = MINUS_Y.add(MINUS_X);
        Assert.assertEquals(expected, actual);

    }

    @Test
    public void testSubtract() {
        FieldTuple<Decimal64> expected, actual;

        expected = new FieldTuple<Decimal64>(X.subtract(Y), X.subtract(Y));
        actual = PLUS_X.subtract(PLUS_Y);
        Assert.assertEquals(expected, actual);

        expected = new FieldTuple<Decimal64>(X.subtract(Y.negate()), X.subtract(Y.negate()));
        actual = PLUS_X.subtract(MINUS_Y);
        Assert.assertEquals(expected, actual);

        expected = new FieldTuple<Decimal64>(X.negate().subtract(Y), X.negate().subtract(Y));
        actual = MINUS_X.subtract(PLUS_Y);
        Assert.assertEquals(expected, actual);

        expected = new FieldTuple<Decimal64>(X.negate().subtract(Y.negate()), X.negate().subtract(Y.negate()));
        actual = MINUS_X.subtract(MINUS_Y);
        Assert.assertEquals(expected, actual);

    }

    @Test
    public void testNegate() {
        FieldTuple<Decimal64> expected, actual;

        expected = MINUS_X;
        actual = PLUS_X.negate();
        Assert.assertEquals(expected, actual);

        expected = PLUS_X;
        actual = MINUS_X.negate();
        Assert.assertEquals(expected, actual);

        expected = MINUS_ZERO;
        actual = PLUS_ZERO.negate();
        Assert.assertEquals(expected, actual);

        expected = PLUS_ZERO;
        actual = MINUS_ZERO.negate();
        Assert.assertEquals(expected, actual);

    }

    @Test
    public void testMultiply() {
        FieldTuple<Decimal64> expected, actual;

        expected = new FieldTuple<Decimal64>(X.multiply(Y), X.multiply(Y));
        actual = PLUS_X.multiply(PLUS_Y);
        Assert.assertEquals(expected, actual);
        actual = PLUS_Y.multiply(PLUS_X);
        Assert.assertEquals(expected, actual);

        expected = new FieldTuple<Decimal64>(X.multiply(Y.negate()), X.multiply(Y.negate()));
        actual = PLUS_X.multiply(MINUS_Y);
        Assert.assertEquals(expected, actual);
        actual = MINUS_Y.multiply(PLUS_X);
        Assert.assertEquals(expected, actual);

        expected = new FieldTuple<Decimal64>(X.negate().multiply(Y.negate()), X.negate().multiply(Y.negate()));
        actual = MINUS_X.multiply(MINUS_Y);
        Assert.assertEquals(expected, actual);
        actual = MINUS_Y.multiply(MINUS_X);
        Assert.assertEquals(expected, actual);

    }

    @Test
    public void testDivide() {
        FieldTuple<Decimal64> expected, actual;

        expected = new FieldTuple<Decimal64>(X.divide(Y), X.divide(Y));
        actual = PLUS_X.divide(PLUS_Y);
        Assert.assertEquals(expected, actual);

        expected = new FieldTuple<Decimal64>(X.divide(Y.negate()), X.divide(Y.negate()));
        actual = PLUS_X.divide(MINUS_Y);
        Assert.assertEquals(expected, actual);

        expected = new FieldTuple<Decimal64>(X.negate().divide(Y), X.negate().divide(Y));
        actual = MINUS_X.divide(PLUS_Y);
        Assert.assertEquals(expected, actual);

        expected = new FieldTuple<Decimal64>(X.negate().divide(Y.negate()), X.negate().divide(Y.negate()));
        actual = MINUS_X.divide(MINUS_Y);
        Assert.assertEquals(expected, actual);

    }

    @Test
    public void testReciprocal() {
        FieldTuple<Decimal64> expected, actual;

        expected = new FieldTuple<Decimal64>(X.reciprocal(), X.reciprocal());
        actual = PLUS_X.reciprocal();
        Assert.assertEquals(expected, actual);

        expected = new FieldTuple<Decimal64>(X.negate().reciprocal(), X.negate().reciprocal());
        actual = MINUS_X.reciprocal();
        Assert.assertEquals(expected, actual);

        expected = PLUS_ZERO;
        actual = new FieldTuple<Decimal64>(new Decimal64(Double.POSITIVE_INFINITY), new Decimal64(Double.POSITIVE_INFINITY)).reciprocal();
        Assert.assertEquals(expected, actual);

        expected = MINUS_ZERO;
        actual = new FieldTuple<Decimal64>(new Decimal64(Double.NEGATIVE_INFINITY), new Decimal64(Double.NEGATIVE_INFINITY)).reciprocal();
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testToDegreesDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            for (double x = 0.1; x < 1.2; x += 0.001) {
                final Decimal64 dec64 = new Decimal64(x);
                FieldTuple<Decimal64> value = new FieldTuple<Decimal64>(dec64, dec64);
                Assert.assertEquals(FastMath.toDegrees(x), value.toDegrees().getReal(), epsilon);
            }
        }
    }

    @Test
    public void testToRadiansDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            for (double x = 0.1; x < 1.2; x += 0.001) {
                final Decimal64 dec64 = new Decimal64(x);
                FieldTuple<Decimal64> value = new FieldTuple<Decimal64>(dec64, dec64);
                Assert.assertEquals(FastMath.toRadians(x), value.toRadians().getReal(), epsilon);
            }
        }
    }

    @Test
    public void testDegRad() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            final Decimal64 dec64 = new Decimal64(x);
            FieldTuple<Decimal64> value = new FieldTuple<Decimal64>(dec64, dec64);
            FieldTuple<Decimal64> rebuilt = value.toDegrees().toRadians();
            FieldTuple<Decimal64> zero = rebuilt.subtract(value);
            Assert.assertEquals(zero.getReal(), 0, 3.0e-16);
        }
    }

}
