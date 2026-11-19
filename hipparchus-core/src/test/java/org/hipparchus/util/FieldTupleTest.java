/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.util;

import org.hipparchus.CalculusFieldElementAbstractTest;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

public class FieldTupleTest extends CalculusFieldElementAbstractTest<FieldTuple<Binary64>> {
    public static final Binary64 X = new Binary64(1.2345);

    public static final FieldTuple<Binary64> PLUS_X = new FieldTuple<>(X, X);

    public static final FieldTuple<Binary64> MINUS_X = new FieldTuple<>(X.negate(), X.negate());

    public static final Binary64 Y = new Binary64(6.789);

    public static final FieldTuple<Binary64> PLUS_Y = new FieldTuple<>(Y, Y);

    public static final FieldTuple<Binary64> MINUS_Y = new FieldTuple<>(Y.negate(), Y.negate());

    public static final FieldTuple<Binary64> PLUS_ZERO = new FieldTuple<>(new Binary64(0.0), new Binary64(0.0));

    public static final FieldTuple<Binary64> MINUS_ZERO = new FieldTuple<>(new Binary64(-0.0), new Binary64(-0.0));

    @Override
    protected FieldTuple<Binary64> build(final double x) {
        return new FieldTuple<>(new Binary64(x), new Binary64(x));
    }

    @Test
    void testNorm() {
        final FieldTuple<Binary64> oneToFive = new FieldTuple<>(Binary64.ONE, new Binary64(2));
        assertEquals(3, oneToFive.norm());
    }

    @Test
    void testComponents() {
        assertEquals(2, PLUS_ZERO.getDimension());
        final FieldTuple<Binary64> oneToFive = new FieldTuple<>(new Binary64(1),
                                                                          new Binary64(2),
                                                                          new Binary64(3),
                                                                          new Binary64(4),
                                                                          new Binary64(5));
        assertEquals(5, oneToFive.getDimension());
        assertArrayEquals(new Binary64[] {
            new Binary64(1), new Binary64(2), new Binary64(3), new Binary64(4), new Binary64(5) },
                                 oneToFive.getComponents());
        assertEquals(1, oneToFive.getComponent(0).getReal(), 1.0e-10);
        assertEquals(2, oneToFive.getComponent(1).getReal(), 1.0e-10);
        assertEquals(3, oneToFive.getComponent(2).getReal(), 1.0e-10);
        assertEquals(4, oneToFive.getComponent(3).getReal(), 1.0e-10);
        assertEquals(5, oneToFive.getComponent(4).getReal(), 1.0e-10);
    }

    @Test
    void testEquals() {
        assertNotEquals(PLUS_ZERO, null);
        assertEquals(PLUS_ZERO, PLUS_ZERO);
        assertEquals(PLUS_X,    PLUS_X);
        assertEquals(PLUS_Y,    PLUS_Y);
        assertEquals(MINUS_X,   MINUS_X);
        assertEquals(MINUS_Y,   MINUS_Y);
        assertNotEquals(PLUS_X,  new FieldTuple<>(new Binary64(1),
                                                         new Binary64(2),
                                                         new Binary64(3),
                                                         new Binary64(4),
                                                         new Binary64(5)));
        assertNotEquals(PLUS_X,  new FieldTuple<>(PLUS_X.getComponent(0), new Binary64(999.999)));
        assertNotEquals(null, PLUS_ZERO.getField());
        assertNotEquals(PLUS_X.getField(),  new FieldTuple<>(new Binary64(1),
                        new Binary64(2),
                        new Binary64(3),
                        new Binary64(4),
                        new Binary64(5)).getField());
        assertEquals(PLUS_ZERO.getField(), MINUS_Y.getField());
    }

    @Test
    void testHashcode() {
        assertEquals(-1264241693, PLUS_ZERO.getField().hashCode());
        assertEquals(1492525478, PLUS_ZERO.hashCode());
        assertEquals(1492525478, MINUS_ZERO.hashCode());
        assertEquals(211878758, PLUS_X.hashCode());
        assertEquals(211878758, MINUS_X.hashCode());
        assertEquals(1460019622, new FieldTuple<>(new Binary64(1), new Binary64(2)).hashCode());
        assertEquals(1491476902, new FieldTuple<>(new Binary64(2), new Binary64(1)).hashCode());
    }

    @Test
    void testAdd() {
        FieldTuple<Binary64> expected, actual;

        expected = new FieldTuple<>(X.add(Y), X.add(Y));
        actual = PLUS_X.add(PLUS_Y);
        assertEquals(expected, actual);
        actual = PLUS_Y.add(PLUS_X);
        assertEquals(expected, actual);

        expected = new FieldTuple<>(X.add(Y.negate()), X.add(Y.negate()));
        actual = PLUS_X.add(MINUS_Y);
        assertEquals(expected, actual);
        actual = MINUS_Y.add(PLUS_X);
        assertEquals(expected, actual);

        expected = new FieldTuple<>(X.negate().add(Y.negate()), X.negate().add(Y.negate()));
        actual = MINUS_X.add(MINUS_Y);
        assertEquals(expected, actual);
        actual = MINUS_Y.add(MINUS_X);
        assertEquals(expected, actual);

    }

    @Test
    void testSubtract() {
        FieldTuple<Binary64> expected, actual;

        expected = new FieldTuple<>(X.subtract(Y), X.subtract(Y));
        actual = PLUS_X.subtract(PLUS_Y);
        assertEquals(expected, actual);

        expected = new FieldTuple<>(X.subtract(Y.negate()), X.subtract(Y.negate()));
        actual = PLUS_X.subtract(MINUS_Y);
        assertEquals(expected, actual);

        expected = new FieldTuple<>(X.negate().subtract(Y), X.negate().subtract(Y));
        actual = MINUS_X.subtract(PLUS_Y);
        assertEquals(expected, actual);

        expected = new FieldTuple<>(X.negate().subtract(Y.negate()), X.negate().subtract(Y.negate()));
        actual = MINUS_X.subtract(MINUS_Y);
        assertEquals(expected, actual);

    }

    @Test
    void testNegate() {
        FieldTuple<Binary64> expected, actual;

        expected = MINUS_X;
        actual = PLUS_X.negate();
        assertEquals(expected, actual);

        expected = PLUS_X;
        actual = MINUS_X.negate();
        assertEquals(expected, actual);

        expected = MINUS_ZERO;
        actual = PLUS_ZERO.negate();
        assertEquals(expected, actual);

        expected = PLUS_ZERO;
        actual = MINUS_ZERO.negate();
        assertEquals(expected, actual);

    }

    @Test
    void testMultiply() {
        FieldTuple<Binary64> expected, actual;

        expected = new FieldTuple<>(X.multiply(Y), X.multiply(Y));
        actual = PLUS_X.multiply(PLUS_Y);
        assertEquals(expected, actual);
        actual = PLUS_Y.multiply(PLUS_X);
        assertEquals(expected, actual);

        expected = new FieldTuple<>(X.multiply(Y.negate()), X.multiply(Y.negate()));
        actual = PLUS_X.multiply(MINUS_Y);
        assertEquals(expected, actual);
        actual = MINUS_Y.multiply(PLUS_X);
        assertEquals(expected, actual);

        expected = new FieldTuple<>(X.negate().multiply(Y.negate()), X.negate().multiply(Y.negate()));
        actual = MINUS_X.multiply(MINUS_Y);
        assertEquals(expected, actual);
        actual = MINUS_Y.multiply(MINUS_X);
        assertEquals(expected, actual);

    }

    @Test
    void testDivide() {
        FieldTuple<Binary64> expected, actual;

        expected = new FieldTuple<>(X.divide(Y), X.divide(Y));
        actual = PLUS_X.divide(PLUS_Y);
        assertEquals(expected, actual);

        expected = new FieldTuple<>(X.divide(Y.negate()), X.divide(Y.negate()));
        actual = PLUS_X.divide(MINUS_Y);
        assertEquals(expected, actual);

        expected = new FieldTuple<>(X.negate().divide(Y), X.negate().divide(Y));
        actual = MINUS_X.divide(PLUS_Y);
        assertEquals(expected, actual);

        expected = new FieldTuple<>(X.negate().divide(Y.negate()), X.negate().divide(Y.negate()));
        actual = MINUS_X.divide(MINUS_Y);
        assertEquals(expected, actual);

    }

    @Test
    void testReciprocal() {
        FieldTuple<Binary64> expected, actual;

        expected = new FieldTuple<>(X.reciprocal(), X.reciprocal());
        actual = PLUS_X.reciprocal();
        assertEquals(expected, actual);

        expected = new FieldTuple<>(X.negate().reciprocal(), X.negate().reciprocal());
        actual = MINUS_X.reciprocal();
        assertEquals(expected, actual);

        expected = PLUS_ZERO;
        actual = new FieldTuple<>(new Binary64(Double.POSITIVE_INFINITY), new Binary64(Double.POSITIVE_INFINITY)).reciprocal();
        assertEquals(expected, actual);

        expected = MINUS_ZERO;
        actual = new FieldTuple<>(new Binary64(Double.NEGATIVE_INFINITY), new Binary64(Double.NEGATIVE_INFINITY)).reciprocal();
        assertEquals(expected, actual);
    }

    @Test
    void testToDegreesDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            for (double x = 0.1; x < 1.2; x += 0.001) {
                final Binary64 dec64 = new Binary64(x);
                FieldTuple<Binary64> value = new FieldTuple<>(dec64, dec64);
                assertEquals(FastMath.toDegrees(x), value.toDegrees().getReal(), epsilon);
            }
        }
    }

    @Test
    void testToRadiansDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            for (double x = 0.1; x < 1.2; x += 0.001) {
                final Binary64 dec64 = new Binary64(x);
                FieldTuple<Binary64> value = new FieldTuple<>(dec64, dec64);
                assertEquals(FastMath.toRadians(x), value.toRadians().getReal(), epsilon);
            }
        }
    }

    @Test
    void testDegRad() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            final Binary64 dec64 = new Binary64(x);
            FieldTuple<Binary64> value = new FieldTuple<>(dec64, dec64);
            FieldTuple<Binary64> rebuilt = value.toDegrees().toRadians();
            FieldTuple<Binary64> zero = rebuilt.subtract(value);
            assertEquals(0, zero.getReal(), 3.0e-16);
        }
    }

    @Test
    public void testIsSmallDimensionMismatch() {
        try {
            new FieldTuple<>(new Binary64(0.0)).isSmall(new FieldTuple<>(new Binary64(0.0), new Binary64(1.0)), 1.0e-15);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, e.getSpecifier());
        }
    }

    @Test
    public void testIsSmall() {
        Assertions.assertTrue(new FieldTuple<>(new Binary64(1.0e-13), new Binary64(2.5e-13)).
                              isSmall(new FieldTuple<>(new Binary64(1.0), new Binary64(1.0)), 1.0e-12));
        Assertions.assertFalse(new FieldTuple<>(new Binary64(1.0e-11), new Binary64(2.5e-13)).
                               isSmall(new FieldTuple<>(new Binary64(1.0), new Binary64(1.0)), 1.0e-12));
        Assertions.assertFalse(new FieldTuple<>(new Binary64(1.0e-13), new Binary64(2.5e-11)).
                               isSmall(new FieldTuple<>(new Binary64(1.0), new Binary64(1.0)), 1.0e-12));
    }

}
