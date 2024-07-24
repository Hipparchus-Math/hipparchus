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

package org.hipparchus.geometry.euclidean.oned;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.geometry.Space;
import org.hipparchus.geometry.euclidean.twod.FieldVector2D;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Locale;

public class Vector1DTest {
    @Test
    public void testConstructors() throws MathIllegalArgumentException {
        checkVector(new Vector1D(3, new Vector1D(FastMath.PI / 3)),
                    FastMath.PI);
        checkVector(new Vector1D(2, Vector1D.ONE, -3, new Vector1D(2)),
                    -4);
        checkVector(new Vector1D(2, Vector1D.ONE,
                                 5, new Vector1D(2),
                                 -3, new Vector1D(3)),
                    3);
        checkVector(new Vector1D(2, Vector1D.ONE,
                                 5, new Vector1D(2),
                                 5, new Vector1D(-2),
                                 -3, new Vector1D(-3)),
                    11);
    }

    @Test
    public void testSpace() {
        Space space = new Vector1D(1).getSpace();
        Assertions.assertEquals(1, space.getDimension());
    }

    @Test
    public void testZero() {
        Assertions.assertEquals(0, new Vector1D(1).getZero().getNorm(), 1.0e-15);
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    public void testEquals() {
        Vector1D u1 = new Vector1D(1);
        Vector1D u2 = new Vector1D(1);
        Assertions.assertEquals(u1, u1);
        Assertions.assertEquals(u1, u2);
        Assertions.assertNotEquals(u1, FieldVector2D.getPlusI(Binary64Field.getInstance()));
        Assertions.assertNotEquals(u1, new Vector1D(1 + 10 * Precision.EPSILON));
        Assertions.assertEquals(new Vector1D(Double.NaN), new Vector1D(Double.NaN));
        Assertions.assertEquals(Vector1D.NaN, Vector1D.NaN);
    }

    @Test
    public void testEqualsIeee754() {
        Vector1D u1 = new Vector1D(1);
        Vector1D u2 = new Vector1D(1);
        Assertions.assertTrue(u1.equalsIeee754(u1));
        Assertions.assertTrue(u1.equalsIeee754(u2));
        Assertions.assertFalse(u1.equalsIeee754(FieldVector2D.getPlusI(Binary64Field.getInstance())));
        Assertions.assertFalse(u1.equalsIeee754(new Vector1D(1 + 10 * Precision.EPSILON)));
        Assertions.assertFalse(new Vector1D(Double.NaN).equalsIeee754(new Vector1D(Double.NaN)));
        Assertions.assertFalse(Vector1D.NaN.equalsIeee754(Vector1D.NaN));
        Assertions.assertFalse(Vector1D.NaN.equalsIeee754(Vector1D.NaN));
    }

    @Test
    public void testHash() {
        Assertions.assertEquals(new Vector1D(Double.NaN).hashCode(), new Vector1D(Double.NaN).hashCode());
        Vector1D u = new Vector1D(1);
        Vector1D v = new Vector1D(1 + 10 * Precision.EPSILON);
        Assertions.assertTrue(u.hashCode() != v.hashCode());
    }

    @Test
    public void testInfinite() {
        Assertions.assertTrue(new Vector1D(Double.NEGATIVE_INFINITY).isInfinite());
        Assertions.assertTrue(new Vector1D(Double.POSITIVE_INFINITY).isInfinite());
        Assertions.assertFalse(new Vector1D(1).isInfinite());
        Assertions.assertFalse(new Vector1D(Double.NaN).isInfinite());
    }

    @Test
    public void testNaN() {
        Assertions.assertTrue(new Vector1D(Double.NaN).isNaN());
        Assertions.assertFalse(new Vector1D(1).isNaN());
        Assertions.assertFalse(new Vector1D(Double.NEGATIVE_INFINITY).isNaN());
    }

    @Test
    public void testToString() {
        Assertions.assertEquals("{3}", new Vector1D(3).toString());
        NumberFormat format = new DecimalFormat("0.000", new DecimalFormatSymbols(Locale.US));
        Assertions.assertEquals("{3.000}", new Vector1D(3).toString(format));
    }

    @Test
    public void testCoordinates() {
        Vector1D v = new Vector1D(1);
        Assertions.assertTrue(FastMath.abs(v.getX() - 1) < 1.0e-12);
    }

    @Test
    public void testNorm1() {
        Assertions.assertEquals(0.0, Vector1D.ZERO.getNorm1(), 0);
        Assertions.assertEquals(6.0, new Vector1D(6).getNorm1(), 0);
    }

    @Test
    public void testNorm() {
        Assertions.assertEquals(0.0, Vector1D.ZERO.getNorm(), 0);
        Assertions.assertEquals(3.0, new Vector1D(-3).getNorm(), 1.0e-12);
    }

    @Test
    public void testNormSq() {
        Assertions.assertEquals(0.0, new Vector1D(0).getNormSq(), 0);
        Assertions.assertEquals(9.0, new Vector1D(-3).getNormSq(), 1.0e-12);
    }

    @Test
    public void testNormInf() {
        Assertions.assertEquals(0.0, Vector1D.ZERO.getNormInf(), 0);
        Assertions.assertEquals(3.0, new Vector1D(-3).getNormInf(), 0);
    }

    @Test
    public void testDistance1() {
        Vector1D v1 = new Vector1D(1);
        Vector1D v2 = new Vector1D(-4);
        Assertions.assertEquals(0.0, new Vector1D(-1).distance1(new Vector1D(-1)), 0);
        Assertions.assertEquals(5.0, v1.distance1(v2), 1.0e-12);
        Assertions.assertEquals(v1.subtract(v2).getNorm1(), v1.distance1(v2), 1.0e-12);
    }

    @Test
    public void testDistance() {
        Vector1D v1 = new Vector1D(1);
        Vector1D v2 = new Vector1D(-4);
        Assertions.assertEquals(0.0, Vector1D.distance(new Vector1D(-1), new Vector1D(-1)), 0);
        Assertions.assertEquals(5.0, Vector1D.distance(v1, v2), 1.0e-12);
        Assertions.assertEquals(v1.subtract(v2).getNorm(), Vector1D.distance(v1, v2), 1.0e-12);
    }

    @Test
    public void testDistanceSq() {
        Vector1D v1 = new Vector1D(1);
        Vector1D v2 = new Vector1D(-4);
        Assertions.assertEquals(0.0, Vector1D.distanceSq(new Vector1D(-1), new Vector1D(-1)), 0);
        Assertions.assertEquals(25.0, Vector1D.distanceSq(v1, v2), 1.0e-12);
        Assertions.assertEquals(Vector1D.distance(v1, v2) * Vector1D.distance(v1, v2),
                            Vector1D.distanceSq(v1, v2), 1.0e-12);
  }

    @Test
    public void testDistanceInf() {
        Vector1D v1 = new Vector1D(1);
        Vector1D v2 = new Vector1D(-4);
        Assertions.assertEquals(0.0, Vector1D.distanceInf(new Vector1D(-1), new Vector1D(-1)), 0);
        Assertions.assertEquals(5.0, Vector1D.distanceInf(v1, v2), 1.0e-12);
        Assertions.assertEquals(v1.subtract(v2).getNormInf(), Vector1D.distanceInf(v1, v2), 1.0e-12);
    }

    @Test
    public void testSubtract() {
        Vector1D v1 = new Vector1D(1);
        Vector1D v2 = new Vector1D(-3);
        v1 = v1.subtract(v2);
        checkVector(v1, 4);

        checkVector(v2.subtract(v1), -7);
        checkVector(v2.subtract(3, v1), -15);
    }

    @Test
    public void testAdd() {
        Vector1D v1 = new Vector1D(1);
        Vector1D v2 = new Vector1D(-3);
        v1 = v1.add(v2);
        checkVector(v1, -2);

        checkVector(v2.add(v1), -5);
        checkVector(v2.add(3, v1), -9);
    }

    @Test
    public void testScalarProduct() {
        Vector1D v = new Vector1D(1);
        v = v.scalarMultiply(3);
        checkVector(v, 3);

        checkVector(v.scalarMultiply(0.5), 1.5);
    }

    @Test
    public void testNormalize() throws MathRuntimeException {
        Assertions.assertEquals(1.0, new Vector1D(5).normalize().getNorm(), 1.0e-12);
        try {
            Vector1D.ZERO.normalize();
            Assertions.fail("an exception should have been thrown");
        } catch (MathRuntimeException ae) {
            // expected behavior
        }
    }

    @Test
    public void testNegate() {
        checkVector(new Vector1D(0.1).negate(), -0.1);
    }

    @Test
    public void testArithmeticBlending() {

        // Given
        final Vector1D v1 = new Vector1D(1);
        final Vector1D v2 = new Vector1D(2);

        final double blendingValue = 0.7;

        // When
        final Vector1D blendedVector = v1.blendArithmeticallyWith(v2, blendingValue);

        // Then
        checkVector(blendedVector, 1.7);
    }

    private void checkVector(Vector1D v, double x) {
        Assertions.assertEquals(x, v.getX(), 1.0e-12);
    }
}
