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
import org.hipparchus.util.Tuple;
import org.junit.Assert;
import org.junit.Test;

public class TupleTest extends ExtendedFieldElementAbstractTest<Tuple> {
    public static final double X = 1.2345;

    public static final Tuple PLUS_X = new Tuple(X, X);

    public static final Tuple MINUS_X = new Tuple(-X, -X);

    public static final double Y = 6.789;

    public static final Tuple PLUS_Y = new Tuple(Y, Y);

    public static final Tuple MINUS_Y = new Tuple(-Y, -Y);

    public static final Tuple PLUS_ZERO = new Tuple(0.0, 0.0);

    public static final Tuple MINUS_ZERO = new Tuple(-0.0, -0.0);

    @Override
    protected Tuple build(final double x) {
        return new Tuple(x, x);
    }

    @Test
    public void testComponents() {
        Assert.assertEquals(2, PLUS_ZERO.getDimension());
        final Tuple oneToFive = new Tuple(1, 2, 3, 4, 5);
        Assert.assertEquals(5, oneToFive.getDimension());
        Assert.assertArrayEquals(new double[] { 1,  2,  3, 4, 5 }, oneToFive.getComponents(), 1.0e-10);
        Assert.assertEquals(1, oneToFive.getComponent(0), 1.0e-10);
        Assert.assertEquals(2, oneToFive.getComponent(1), 1.0e-10);
        Assert.assertEquals(3, oneToFive.getComponent(2), 1.0e-10);
        Assert.assertEquals(4, oneToFive.getComponent(3), 1.0e-10);
        Assert.assertEquals(5, oneToFive.getComponent(4), 1.0e-10);
    }

    @Test
    public void testEquals() {
        Assert.assertNotEquals(PLUS_ZERO, null);
        Assert.assertEquals(PLUS_ZERO, PLUS_ZERO);
        Assert.assertEquals(PLUS_X,    PLUS_X);
        Assert.assertEquals(PLUS_Y,    PLUS_Y);
        Assert.assertEquals(MINUS_X,   MINUS_X);
        Assert.assertEquals(MINUS_Y,   MINUS_Y);
        Assert.assertNotEquals(PLUS_X,  new Tuple(1, 2, 3, 4, 5));
        Assert.assertNotEquals(PLUS_X,  new Tuple(PLUS_X.getComponent(0), 999.999));
        Assert.assertNotEquals(PLUS_ZERO.getField(), null);
        Assert.assertNotEquals(PLUS_X.getField(), new Tuple(1, 2, 3, 4, 5).getField());
        Assert.assertEquals(PLUS_ZERO.getField(), MINUS_Y.getField());
    }

    @Test
    public void testHashcode() {
        Assert.assertEquals(1718765887, PLUS_ZERO.getField().hashCode());
        Assert.assertEquals(884058117, PLUS_ZERO.hashCode());
        Assert.assertEquals(884058117, MINUS_ZERO.hashCode());
        Assert.assertEquals(-396588603, PLUS_X.hashCode());
        Assert.assertEquals(-396588603, MINUS_X.hashCode());
        Assert.assertEquals(851552261, new Tuple(1, 2).hashCode());
        Assert.assertEquals(883009541, new Tuple(2, 1).hashCode());
    }

    @Test
    public void testAdd() {
        Tuple expected, actual;

        expected = new Tuple(X + Y, X + Y);
        actual = PLUS_X.add(PLUS_Y);
        Assert.assertEquals(expected, actual);
        actual = PLUS_Y.add(PLUS_X);
        Assert.assertEquals(expected, actual);

        expected = new Tuple(X + (-Y), X + (-Y));
        actual = PLUS_X.add(MINUS_Y);
        Assert.assertEquals(expected, actual);
        actual = MINUS_Y.add(PLUS_X);
        Assert.assertEquals(expected, actual);

        expected = new Tuple((-X) + (-Y), (-X) + (-Y));
        actual = MINUS_X.add(MINUS_Y);
        Assert.assertEquals(expected, actual);
        actual = MINUS_Y.add(MINUS_X);
        Assert.assertEquals(expected, actual);

    }

    @Test
    public void testSubtract() {
        Tuple expected, actual;

        expected = new Tuple(X - Y, X - Y);
        actual = PLUS_X.subtract(PLUS_Y);
        Assert.assertEquals(expected, actual);

        expected = new Tuple(X - (-Y), X - (-Y));
        actual = PLUS_X.subtract(MINUS_Y);
        Assert.assertEquals(expected, actual);

        expected = new Tuple((-X) - Y, (-X) - Y);
        actual = MINUS_X.subtract(PLUS_Y);
        Assert.assertEquals(expected, actual);

        expected = new Tuple((-X) - (-Y), (-X) - (-Y));
        actual = MINUS_X.subtract(MINUS_Y);
        Assert.assertEquals(expected, actual);

    }

    @Test
    public void testNegate() {
        Tuple expected, actual;

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
        Tuple expected, actual;

        expected = new Tuple(X * Y, X * Y);
        actual = PLUS_X.multiply(PLUS_Y);
        Assert.assertEquals(expected, actual);
        actual = PLUS_Y.multiply(PLUS_X);
        Assert.assertEquals(expected, actual);

        expected = new Tuple(X * (-Y), X * (-Y));
        actual = PLUS_X.multiply(MINUS_Y);
        Assert.assertEquals(expected, actual);
        actual = MINUS_Y.multiply(PLUS_X);
        Assert.assertEquals(expected, actual);

        expected = new Tuple((-X) * (-Y), (-X) * (-Y));
        actual = MINUS_X.multiply(MINUS_Y);
        Assert.assertEquals(expected, actual);
        actual = MINUS_Y.multiply(MINUS_X);
        Assert.assertEquals(expected, actual);

    }

    @Test
    public void testDivide() {
        Tuple expected, actual;

        expected = new Tuple(X / Y, X / Y);
        actual = PLUS_X.divide(PLUS_Y);
        Assert.assertEquals(expected, actual);

        expected = new Tuple(X / (-Y), X / (-Y));
        actual = PLUS_X.divide(MINUS_Y);
        Assert.assertEquals(expected, actual);

        expected = new Tuple((-X) / Y, (-X) / Y);
        actual = MINUS_X.divide(PLUS_Y);
        Assert.assertEquals(expected, actual);

        expected = new Tuple((-X) / (-Y), (-X) / (-Y));
        actual = MINUS_X.divide(MINUS_Y);
        Assert.assertEquals(expected, actual);

    }

    @Test
    public void testReciprocal() {
        Tuple expected, actual;

        expected = new Tuple(1.0 / X, 1.0 / X);
        actual = PLUS_X.reciprocal();
        Assert.assertEquals(expected, actual);

        expected = new Tuple(1.0 / (-X), 1.0 / (-X));
        actual = MINUS_X.reciprocal();
        Assert.assertEquals(expected, actual);

        expected = PLUS_ZERO;
        actual = new Tuple(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY).reciprocal();
        Assert.assertEquals(expected, actual);

        expected = MINUS_ZERO;
        actual = new Tuple(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY).reciprocal();
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testToDegreesDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            for (double x = 0.1; x < 1.2; x += 0.001) {
                Tuple value = new Tuple(x, x);
                Assert.assertEquals(FastMath.toDegrees(x), value.toDegrees().getReal(), epsilon);
            }
        }
    }

    @Test
    public void testToRadiansDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            for (double x = 0.1; x < 1.2; x += 0.001) {
                Tuple value = new Tuple(x, x);
                Assert.assertEquals(FastMath.toRadians(x), value.toRadians().getReal(), epsilon);
            }
        }
    }

    @Test
    public void testDegRad() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            Tuple value = new Tuple(x, x);
            Tuple rebuilt = value.toDegrees().toRadians();
            Tuple zero = rebuilt.subtract(value);
            Assert.assertEquals(zero.getReal(), 0, 3.0e-16);
        }
    }

}
