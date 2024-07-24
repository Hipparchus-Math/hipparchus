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
package org.hipparchus.geometry.euclidean.twod;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.geometry.LocalizedGeometryFormats;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.SinCos;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.text.NumberFormat;
import java.util.Locale;

public class Vector2DTest {

    @Test public void testConstructors() {
        final double p40 =  4.0;
        final double p20 =  2.0;
        final double p25 =  2.5;
        final double p10 =  1.0;
        final double m05 = -0.5;
        final double m30 = -3.0;
        check(new Vector2D(p25, m05), 2.5, -0.5, 1.0e-15);
        final double[] a = new double[2];
        a[0] = 1.0;
        a[1] = 0.0;
        check(new Vector2D(a), 1.0, 0.0, 1.0e-15);
        try {
            new Vector2D(new double[3]);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            Assertions.assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            Assertions.assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
        check(new Vector2D(p20, new Vector2D(p25, m05)), 5.0, -1.0, 1.0e-15);
        check(new Vector2D(p20, new Vector2D(2.5, -0.5)), 5.0, -1.0, 1.0e-15);
        check(new Vector2D(2, new Vector2D(p25, m05)), 5.0, -1.0, 1.0e-15);
        check(new Vector2D(p20, new Vector2D(p25, m05), m30, new Vector2D(m05, p40)),
              6.5, -13.0, 1.0e-15);
        check(new Vector2D(p20, new Vector2D(2.5, -0.5), m30, new Vector2D(-0.5, 4.0)),
              6.5, -13.0, 1.0e-15);
        check(new Vector2D(2.0, new Vector2D(p25, m05), -3.0, new Vector2D(m05, p40)),
              6.5, -13.0, 1.0e-15);
        check(new Vector2D(p20, new Vector2D(p25, m05), m30, new Vector2D(m05, p40),
                           p40, new Vector2D(p25, m30)),
              16.5, -25.0, 1.0e-15);
        check(new Vector2D(p20, new Vector2D(2.5, -0.5), m30, new Vector2D(-0.5, 4.0),
                           p40, new Vector2D(2.5, -3.0)),
              16.5, -25.0, 1.0e-15);
        check(new Vector2D(2.0, new Vector2D(p25, m05), -3.0, new Vector2D(m05, p40),
                           4.0, new Vector2D(p25, m30)),
              16.5, -25.0, 1.0e-15);
        check(new Vector2D(p20, new Vector2D(p25, m05), m30, new Vector2D(m05, p40),
                           p40, new Vector2D(p25, m30), p10, new Vector2D(p10, p10)),
              17.5, -24.0, 1.0e-15);
        check(new Vector2D(p20, new Vector2D(2.5, -0.5), m30, new Vector2D(-0.5, 4.0),
                           p40, new Vector2D(2.5, -3.0), p10, new Vector2D(1.0, 1.0)),
              17.5, -24.0, 1.0e-15);
        check(new Vector2D(2.0, new Vector2D(p25, m05), -3.0, new Vector2D(m05, p40),
                           4.0, new Vector2D(p25, m30),  1.0, new Vector2D(p10, p10)),
              17.5, -24.0, 1.0e-15);
    }

    @Test public void testConstants() {
        check(Vector2D.ZERO,    0.0,  0.0, 1.0e-15);
        check(Vector2D.PLUS_I,   1.0,  0.0, 1.0e-15);
        check(Vector2D.MINUS_I, -1.0,  0.0, 1.0e-15);
        check(Vector2D.PLUS_J,   0.0,  1.0, 1.0e-15);
        check(Vector2D.MINUS_J,  0.0, -1.0, 1.0e-15);
        check(Vector2D.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0e-15);
        check(Vector2D.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, 1.0e-15);
        Assertions.assertTrue(Double.isNaN(Vector2D.NaN.getX()));
        Assertions.assertTrue(Double.isNaN(Vector2D.NaN.getY()));
        Assertions.assertSame(Euclidean2D.getInstance(), Vector2D.NaN.getSpace());
        Assertions.assertSame(Vector2D.ZERO, new Vector2D(1.0, 2.0).getZero());
    }

    @Test public void testToMethods() {
        final Vector2D v = new Vector2D(2.5, -0.5);
        Assertions.assertEquals( 2,   v.toArray().length);
        Assertions.assertEquals( 2.5, v.toArray()[0], 1.0e-15);
        Assertions.assertEquals(-0.5, v.toArray()[1], 1.0e-15);
        Assertions.assertEquals("{2.5; -0.5}", v.toString().replaceAll(",", "."));
        Assertions.assertEquals("{2,5; -0,5}", v.toString(NumberFormat.getInstance(Locale.FRENCH)));
    }

    @Test public void testNorms() {
        final Vector2D v = new Vector2D(3.0, -4.0);
        Assertions.assertEquals( 7.0, v.getNorm1(),   1.0e-15);
        Assertions.assertEquals( 5.0, v.getNorm(),    1.0e-15);
        Assertions.assertEquals(25.0, v.getNormSq(),  1.0e-15);
        Assertions.assertEquals( 4.0, v.getNormInf(), 1.0e-15);
    }

    @Test public void testDistances() {
        final Vector2D u = new Vector2D( 2.0, -2.0);
        final Vector2D v = new Vector2D(-1.0,  2.0);
        Assertions.assertEquals( 7.0, Vector2D.distance1(u, v),   1.0e-15);
        Assertions.assertEquals( 5.0, Vector2D.distance(u, v),    1.0e-15);
        Assertions.assertEquals(25.0, Vector2D.distanceSq(u, v),  1.0e-15);
        Assertions.assertEquals( 4.0, Vector2D.distanceInf(u, v), 1.0e-15);
    }

    @Test public void testAdd() {
        final Vector2D u = new Vector2D( 2.0, -2.0);
        final Vector2D v = new Vector2D(-1.0,  2.0);
        check(u.add(v), 1.0, 0.0, 1.0e-15);
        check(u.add(5.0, v), -3.0, 8.0, 1.0e-15);
    }

    @Test public void testSubtract() {
        final Vector2D u = new Vector2D( 2.0, -2.0);
        final Vector2D v = new Vector2D( 1.0, -2.0);
        check(u.subtract(v), 1.0, 0.0, 1.0e-15);
        check(u.subtract(5, v), -3.0, 8.0, 1.0e-15);
    }

    @Test public void testNormalize() {
        try {
            Vector2D.ZERO.normalize();
            Assertions.fail("an exception should habe been thrown");
        } catch (MathRuntimeException mre) {
            Assertions.assertEquals(LocalizedGeometryFormats.CANNOT_NORMALIZE_A_ZERO_NORM_VECTOR, mre.getSpecifier());
        }
        check(new Vector2D(3, -4).normalize(), 0.6, -0.8, 1.0e-15);
    }

    @Test public void testAngle() {
        try {
            Vector2D.angle(Vector2D.ZERO, Vector2D.PLUS_I);
            Assertions.fail("an exception should habe been thrown");
        } catch (MathRuntimeException mre) {
            Assertions.assertEquals(LocalizedCoreFormats.ZERO_NORM, mre.getSpecifier());
        }
        final double alpha = 0.01;
        final SinCos sc = FastMath.sinCos(alpha);
        Assertions.assertEquals(alpha,
                            Vector2D.angle(new Vector2D(sc.cos(), sc.sin()),
                                           Vector2D.PLUS_I),
                            1.0e-15);
        Assertions.assertEquals(FastMath.PI - alpha,
                            Vector2D.angle(new Vector2D(-sc.cos(), sc.sin()),
                                           Vector2D.PLUS_I),
                            1.0e-15);
        Assertions.assertEquals(0.5 * FastMath.PI - alpha,
                            Vector2D.angle(new Vector2D(sc.sin(), sc.cos()),
                                           Vector2D.PLUS_I),
                            1.0e-15);
        Assertions.assertEquals(0.5 * FastMath.PI + alpha,
                            Vector2D.angle(new Vector2D(-sc.sin(), sc.cos()),
                                           Vector2D.PLUS_I),
                            1.0e-15);
    }


    @Test public void testNegate() {
        check(new Vector2D(3.0, -4.0).negate(), -3.0, 4.0, 1.0e-15);
    }

    @Test public void testScalarMultiply() {
        check(new Vector2D(3.0, -4.0).scalarMultiply(2.0), 6.0, -8.0, 1.0e-15);
    }

    @Test public void testIsNaN() {
        Assertions.assertTrue(new Vector2D(Double.NaN, 0.0).isNaN());
        Assertions.assertTrue(new Vector2D(0.0, Double.NaN).isNaN());
        Assertions.assertTrue(new Vector2D(Double.NaN, Double.NaN).isNaN());
        Assertions.assertTrue(Vector2D.NaN.isNaN());
        Assertions.assertFalse(new Vector2D(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY).isNaN());
        Assertions.assertFalse(Vector2D.MINUS_I.isNaN());
    }

    @Test public void testIsInfinite() {
        Assertions.assertFalse(new Vector2D(Double.NaN, 0.0).isInfinite());
        Assertions.assertTrue(new Vector2D(Double.POSITIVE_INFINITY, 0.0).isInfinite());
        Assertions.assertTrue(new Vector2D(0.0, Double.POSITIVE_INFINITY).isInfinite());
        Assertions.assertTrue(new Vector2D(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY).isInfinite());
        Assertions.assertTrue(new Vector2D(Double.NEGATIVE_INFINITY, 0.0).isInfinite());
        Assertions.assertTrue(new Vector2D(0.0, Double.NEGATIVE_INFINITY).isInfinite());
        Assertions.assertTrue(new Vector2D(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY).isInfinite());
        Assertions.assertFalse(Vector2D.NaN.isInfinite());
        Assertions.assertFalse(Vector2D.MINUS_I.isInfinite());
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test public void testEquals() {
        final Vector2D u1 = Vector2D.PLUS_I;
        final Vector2D u2 = Vector2D.MINUS_I.negate();
        final Vector2D v1 = new Vector2D(1.0, 0.001);
        final Vector2D v2 = new Vector2D(0.001, 1.0);
        Assertions.assertEquals(u1, u1);
        Assertions.assertEquals(u1, u2);
        Assertions.assertNotEquals(u1, v1);
        Assertions.assertNotEquals(u1, v2);
        Assertions.assertNotEquals(u1, FieldVector2D.getPlusI(Binary64Field.getInstance()));
        Assertions.assertEquals(Vector2D.NaN, new Vector2D(Double.NaN, u1));
        Assertions.assertNotEquals(Vector2D.NaN, u1);
        Assertions.assertNotEquals(Vector2D.NaN, v2);
    }

    @Test public void testEqualsIeee754() {
        final Vector2D u1 = Vector2D.PLUS_I;
        final Vector2D u2 = Vector2D.MINUS_I.negate();
        final Vector2D v1 = new Vector2D(1.0, 0.001);
        final Vector2D v2 = new Vector2D(0.001, 1.0);
        Assertions.assertTrue(u1.equalsIeee754(u1));
        Assertions.assertTrue(u1.equalsIeee754(u2));
        Assertions.assertFalse(u1.equalsIeee754(v1));
        Assertions.assertFalse(u1.equalsIeee754(v2));
        Assertions.assertFalse(u1.equalsIeee754(FieldVector2D.getPlusI(Binary64Field.getInstance())));
        Assertions.assertFalse(new Vector2D(Double.NaN, u1).equalsIeee754(Vector2D.NaN));
        Assertions.assertFalse(u1.equalsIeee754(Vector2D.NaN));
        Assertions.assertFalse(Vector2D.NaN.equalsIeee754(v2));
        Assertions.assertFalse(Vector2D.NaN.equalsIeee754(Vector2D.NaN));
    }

    @Test public void testHashCode() {
        Assertions.assertEquals(542, Vector2D.NaN.hashCode());
        Assertions.assertEquals(1325400064, new Vector2D(1.5, -0.5).hashCode());
    }

    @Test public void testCrossProduct() {
        final double epsilon = 1e-10;

        Vector2D p1 = new Vector2D(1, 1);
        Vector2D p2 = new Vector2D(2, 2);

        Vector2D p3 = new Vector2D(3, 3);
        Assertions.assertEquals(0.0, p3.crossProduct(p1, p2), epsilon);

        Vector2D p4 = new Vector2D(1, 2);
        Assertions.assertEquals(1.0, p4.crossProduct(p1, p2), epsilon);

        Vector2D p5 = new Vector2D(2, 1);
        Assertions.assertEquals(-1.0, p5.crossProduct(p1, p2), epsilon);

    }

    @Test public void testOrientation() {
        Assertions.assertTrue(Vector2D.orientation(new Vector2D(0, 0),
                                               new Vector2D(1, 0),
                                               new Vector2D(1, 1)) > 0);
        Assertions.assertTrue(Vector2D.orientation(new Vector2D(1, 0),
                                               new Vector2D(0, 0),
                                               new Vector2D(1, 1)) < 0);
        Assertions.assertEquals(0.0,
                            Vector2D.orientation(new Vector2D(0, 0),
                                                 new Vector2D(1, 0),
                                                 new Vector2D(1, 0)),
                            1.0e-15);
        Assertions.assertEquals(0.0,
                            Vector2D.orientation(new Vector2D(0, 0),
                                                 new Vector2D(1, 0),
                                                 new Vector2D(2, 0)),
                            1.0e-15);
    }

    @Test
    public void testArithmeticBlending() {

        // Given
        final Vector2D v1 = new Vector2D(1,2);
        final Vector2D v2 = new Vector2D(3,4);

        final double blendingValue = 0.7;

        // When
        final Vector2D blendedVector = v1.blendArithmeticallyWith(v2, blendingValue);

        // Then
        check(blendedVector, 2.4, 3.4, 1e-12);
    }

    private void check(final Vector2D v, final double x, final double y, final double tol) {
        Assertions.assertEquals(x, v.getX(), tol);
        Assertions.assertEquals(y, v.getY(), tol);
    }

}
