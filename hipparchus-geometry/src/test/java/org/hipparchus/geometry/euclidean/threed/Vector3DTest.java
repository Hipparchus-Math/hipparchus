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

package org.hipparchus.geometry.euclidean.threed;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.geometry.Space;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Test;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Locale;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class Vector3DTest {
    @Test
    void testConstructors() throws MathIllegalArgumentException {
        double r = FastMath.sqrt(2) /2;
        checkVector(new Vector3D(2, new Vector3D(FastMath.PI / 3, -FastMath.PI / 4)),
                    r, r * FastMath.sqrt(3), -2 * r);
        checkVector(new Vector3D(2, Vector3D.PLUS_I,
                                 -3, Vector3D.MINUS_K),
                    2, 0, 3);
        checkVector(new Vector3D(2, Vector3D.PLUS_I,
                                 5, Vector3D.PLUS_J,
                                 -3, Vector3D.MINUS_K),
                    2, 5, 3);
        checkVector(new Vector3D(2, Vector3D.PLUS_I,
                                 5, Vector3D.PLUS_J,
                                 5, Vector3D.MINUS_J,
                                 -3, Vector3D.MINUS_K),
                    2, 0, 3);
        checkVector(new Vector3D(new double[] { 2,  5,  -3 }),
                    2, 5, -3);
    }

    @Test
    void testSpace() {
        Space space = new Vector3D(1, 2, 2).getSpace();
        assertEquals(3, space.getDimension());
        assertEquals(2, space.getSubSpace().getDimension());
        Space deserialized = (Space) UnitTestUtils.serializeAndRecover(space);
        assertTrue(space == deserialized);
    }

    @Test
    void testZero() {
        assertEquals(0, new Vector3D(1, 2, 2).getZero().getNorm(), 1.0e-15);
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testEquals() {
        Vector3D u1 = new Vector3D(1, 2, 3);
        Vector3D u2 = new Vector3D(1, 2, 3);
        assertEquals(u1, u1);
        assertEquals(u1, u2);
        assertNotEquals(u1, new Rotation(1, 0, 0, 0, false));
        assertNotEquals(u1, new Vector3D(1, 2, 3 + 10 * Precision.EPSILON));
        assertNotEquals(u1, new Vector3D(1, 2 + 10 * Precision.EPSILON, 3));
        assertNotEquals(u1, new Vector3D(1 + 10 * Precision.EPSILON, 2, 3));
        assertEquals(new Vector3D(0, Double.NaN, 0), new Vector3D(0, 0, Double.NaN));
    }

    @Test
    void testEqualsIeee754() {
        Vector3D u1 = new Vector3D(1, 2, 3);
        Vector3D u2 = new Vector3D(1, 2, 3);
        assertTrue(u1.equalsIeee754(u1));
        assertTrue(u1.equalsIeee754(u2));
        assertFalse(u1.equalsIeee754(new Rotation(1, 0, 0, 0, false)));
        assertFalse(u1.equalsIeee754(new Vector3D(1, 2, 3 + 10 * Precision.EPSILON)));
        assertFalse(u1.equalsIeee754(new Vector3D(1, 2 + 10 * Precision.EPSILON, 3)));
        assertFalse(u1.equalsIeee754(new Vector3D(1 + 10 * Precision.EPSILON, 2, 3)));
        assertFalse(new Vector3D(0, Double.NaN, 0).equalsIeee754(new Vector3D(0, 0, Double.NaN)));
        assertFalse(Vector3D.NaN.equalsIeee754(Vector3D.NaN));
    }

    @Test
    void testHash() {
        assertEquals(new Vector3D(0, Double.NaN, 0).hashCode(), new Vector3D(0, 0, Double.NaN).hashCode());
        Vector3D u = new Vector3D(1, 2, 3);
        Vector3D v = new Vector3D(1, 2, 3 + 10 * Precision.EPSILON);
        assertTrue(u.hashCode() != v.hashCode());
    }

    @Test
    void testInfinite() {
        assertTrue(new Vector3D(1, 1, Double.NEGATIVE_INFINITY).isInfinite());
        assertTrue(new Vector3D(1, Double.NEGATIVE_INFINITY, 1).isInfinite());
        assertTrue(new Vector3D(Double.NEGATIVE_INFINITY, 1, 1).isInfinite());
        assertFalse(new Vector3D(1, 1, 2).isInfinite());
        assertFalse(new Vector3D(1, Double.NaN, Double.NEGATIVE_INFINITY).isInfinite());
    }

    @Test
    void testNaN() {
        assertTrue(new Vector3D(1, 1, Double.NaN).isNaN());
        assertTrue(new Vector3D(1, Double.NaN, 1).isNaN());
        assertTrue(new Vector3D(Double.NaN, 1, 1).isNaN());
        assertFalse(new Vector3D(1, 1, 2).isNaN());
        assertFalse(new Vector3D(1, 1, Double.NEGATIVE_INFINITY).isNaN());
    }

    @Test
    void testToString() {
        assertEquals("{3; 2; 1}", new Vector3D(3, 2, 1).toString());
        NumberFormat format = new DecimalFormat("0.000", new DecimalFormatSymbols(Locale.US));
        assertEquals("{3.000; 2.000; 1.000}", new Vector3D(3, 2, 1).toString(format));
    }

    @Test
    void testWrongDimension() throws MathIllegalArgumentException {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new Vector3D(new double[]{2, 5});
        });
    }

    @Test
    void testCoordinates() {
        Vector3D v = new Vector3D(1, 2, 3);
        assertTrue(FastMath.abs(v.getX() - 1) < 1.0e-12);
        assertTrue(FastMath.abs(v.getY() - 2) < 1.0e-12);
        assertTrue(FastMath.abs(v.getZ() - 3) < 1.0e-12);
        double[] coordinates = v.toArray();
        assertTrue(FastMath.abs(coordinates[0] - 1) < 1.0e-12);
        assertTrue(FastMath.abs(coordinates[1] - 2) < 1.0e-12);
        assertTrue(FastMath.abs(coordinates[2] - 3) < 1.0e-12);
    }

    @Test
    void testNorm1() {
        assertEquals(0.0, Vector3D.ZERO.getNorm1(), 0);
        assertEquals(6.0, new Vector3D(1, -2, 3).getNorm1(), 0);
    }

    @Test
    void testNorm() {
        assertEquals(0.0, Vector3D.ZERO.getNorm(), 0);
        assertEquals(FastMath.sqrt(14), new Vector3D(1, 2, 3).getNorm(), 1.0e-12);
    }

    @Test
    void testNormSq() {
        assertEquals(0.0, new Vector3D(0, 0, 0).getNormSq(), 0);
        assertEquals(14, new Vector3D(1, 2, 3).getNormSq(), 1.0e-12);
    }

    @Test
    void testNormInf() {
        assertEquals(0.0, Vector3D.ZERO.getNormInf(), 0);
        assertEquals(3.0, new Vector3D(1, -2, 3).getNormInf(), 0);
    }

    @Test
    void testDistance1() {
        Vector3D v1 = new Vector3D(1, -2, 3);
        Vector3D v2 = new Vector3D(-4, 2, 0);
        assertEquals(0.0, Vector3D.distance1(Vector3D.MINUS_I, Vector3D.MINUS_I), 0);
        assertEquals(12.0, Vector3D.distance1(v1, v2), 1.0e-12);
        assertEquals(v1.subtract(v2).getNorm1(), Vector3D.distance1(v1, v2), 1.0e-12);
    }

    @Test
    void testDistance() {
        Vector3D v1 = new Vector3D(1, -2, 3);
        Vector3D v2 = new Vector3D(-4, 2, 0);
        assertEquals(0.0, Vector3D.distance(Vector3D.MINUS_I, Vector3D.MINUS_I), 0);
        assertEquals(FastMath.sqrt(50), Vector3D.distance(v1, v2), 1.0e-12);
        assertEquals(v1.subtract(v2).getNorm(), Vector3D.distance(v1, v2), 1.0e-12);
    }

    @Test
    void testDistanceSq() {
        Vector3D v1 = new Vector3D(1, -2, 3);
        Vector3D v2 = new Vector3D(-4, 2, 0);
        assertEquals(0.0, Vector3D.distanceSq(Vector3D.MINUS_I, Vector3D.MINUS_I), 0);
        assertEquals(50.0, Vector3D.distanceSq(v1, v2), 1.0e-12);
        assertEquals(Vector3D.distance(v1, v2) * Vector3D.distance(v1, v2),
                            Vector3D.distanceSq(v1, v2), 1.0e-12);
  }

    @Test
    void testDistanceInf() {
        Vector3D v1 = new Vector3D(1, -2, 3);
        Vector3D v2 = new Vector3D(-4, 2, 0);
        assertEquals(0.0, Vector3D.distanceInf(Vector3D.MINUS_I, Vector3D.MINUS_I), 0);
        assertEquals(5.0, Vector3D.distanceInf(v1, v2), 1.0e-12);
        assertEquals(v1.subtract(v2).getNormInf(), Vector3D.distanceInf(v1, v2), 1.0e-12);
    }

    @Test
    void testSubtract() {
        Vector3D v1 = new Vector3D(1, 2, 3);
        Vector3D v2 = new Vector3D(-3, -2, -1);
        v1 = v1.subtract(v2);
        checkVector(v1, 4, 4, 4);

        checkVector(v2.subtract(v1), -7, -6, -5);
        checkVector(v2.subtract(3, v1), -15, -14, -13);
    }

    @Test
    void testAdd() {
        Vector3D v1 = new Vector3D(1, 2, 3);
        Vector3D v2 = new Vector3D(-3, -2, -1);
        v1 = v1.add(v2);
        checkVector(v1, -2, 0, 2);

        checkVector(v2.add(v1), -5, -2, 1);
        checkVector(v2.add(3, v1), -9, -2, 5);
    }

    @Test
    void testScalarProduct() {
        Vector3D v = new Vector3D(1, 2, 3);
        v = v.scalarMultiply(3);
        checkVector(v, 3, 6, 9);

        checkVector(v.scalarMultiply(0.5), 1.5, 3, 4.5);
    }

    @Test
    void testVectorialProducts() {
        Vector3D v1 = new Vector3D(2, 1, -4);
        Vector3D v2 = new Vector3D(3, 1, -1);

        assertTrue(FastMath.abs(Vector3D.dotProduct(v1, v2) - 11) < 1.0e-12);

        Vector3D v3 = Vector3D.crossProduct(v1, v2);
        checkVector(v3, 3, -10, -1);

        assertTrue(FastMath.abs(Vector3D.dotProduct(v1, v3)) < 1.0e-12);
        assertTrue(FastMath.abs(Vector3D.dotProduct(v2, v3)) < 1.0e-12);
    }

    @Test
    void testCrossProductCancellation() {
        Vector3D v1 = new Vector3D(9070467121.0, 4535233560.0, 1);
        Vector3D v2 = new Vector3D(9070467123.0, 4535233561.0, 1);
        checkVector(Vector3D.crossProduct(v1, v2), -1, 2, 1);

        double scale    = FastMath.scalb(1.0, 100);
        Vector3D big1   = new Vector3D(scale, v1);
        Vector3D small2 = new Vector3D(1 / scale, v2);
        checkVector(Vector3D.crossProduct(big1, small2), -1, 2, 1);

    }

    @Test
    void testAngular() {
        assertEquals(0,           Vector3D.PLUS_I.getAlpha(), 1.0e-10);
        assertEquals(0,           Vector3D.PLUS_I.getDelta(), 1.0e-10);
        assertEquals(FastMath.PI / 2, Vector3D.PLUS_J.getAlpha(), 1.0e-10);
        assertEquals(0,           Vector3D.PLUS_J.getDelta(), 1.0e-10);
        assertEquals(0,           Vector3D.PLUS_K.getAlpha(), 1.0e-10);
        assertEquals(FastMath.PI / 2, Vector3D.PLUS_K.getDelta(), 1.0e-10);

        Vector3D u = new Vector3D(-1, 1, -1);
        assertEquals(3 * FastMath.PI /4, u.getAlpha(), 1.0e-10);
        assertEquals(-1.0 / FastMath.sqrt(3), FastMath.sin(u.getDelta()), 1.0e-10);
    }

    @Test
    void testAngularSeparation() throws MathRuntimeException {
        Vector3D v1 = new Vector3D(2, -1, 4);

        Vector3D  k = v1.normalize();
        Vector3D  i = k.orthogonal();
        Vector3D v2 = k.scalarMultiply(FastMath.cos(1.2)).add(i.scalarMultiply(FastMath.sin(1.2)));

        assertTrue(FastMath.abs(Vector3D.angle(v1, v2) - 1.2) < 1.0e-12);
  }

    @Test
    void testNormalize() throws MathRuntimeException {
        assertEquals(1.0, new Vector3D(5, -4, 2).normalize().getNorm(), 1.0e-12);
        try {
            Vector3D.ZERO.normalize();
            fail("an exception should have been thrown");
        } catch (MathRuntimeException ae) {
            // expected behavior
        }
    }

    @Test
    void testNegate() {
        checkVector(new Vector3D(0.1, 2.5, 1.3).negate(), -0.1, -2.5, -1.3);
    }

    @Test
    void testOrthogonal() throws MathRuntimeException {
        Vector3D v1 = new Vector3D(0.1, 2.5, 1.3);
        assertEquals(0.0, Vector3D.dotProduct(v1, v1.orthogonal()), 1.0e-12);
        Vector3D v2 = new Vector3D(2.3, -0.003, 7.6);
        assertEquals(0.0, Vector3D.dotProduct(v2, v2.orthogonal()), 1.0e-12);
        Vector3D v3 = new Vector3D(-1.7, 1.4, 0.2);
        assertEquals(0.0, Vector3D.dotProduct(v3, v3.orthogonal()), 1.0e-12);
        Vector3D v4 = new Vector3D(4.2, 0.1, -1.8);
        assertEquals(0.0, Vector3D.dotProduct(v4, v4.orthogonal()), 1.0e-12);
        try {
            new Vector3D(0, 0, 0).orthogonal();
            fail("an exception should have been thrown");
        } catch (MathRuntimeException ae) {
            // expected behavior
        }
    }

    @Test
    void testAngle() throws MathRuntimeException {
        assertEquals(0.22572612855273393616,
                            Vector3D.angle(new Vector3D(1, 2, 3), new Vector3D(4, 5, 6)),
                            1.0e-12);
        assertEquals(7.98595620686106654517199e-8,
                            Vector3D.angle(new Vector3D(1, 2, 3), new Vector3D(2, 4, 6.000001)),
                            1.0e-12);
        assertEquals(3.14159257373023116985197793156,
                            Vector3D.angle(new Vector3D(1, 2, 3), new Vector3D(-2, -4, -6.000001)),
                            1.0e-12);
        try {
            Vector3D.angle(Vector3D.ZERO, Vector3D.PLUS_I);
            fail("an exception should have been thrown");
        } catch (MathRuntimeException ae) {
            // expected behavior
        }
    }

    @Test
    void testAccurateDotProduct() {
        // the following two vectors are nearly but not exactly orthogonal
        // naive dot product (i.e. computing u1.x * u2.x + u1.y * u2.y + u1.z * u2.z
        // leads to a result of 0.0, instead of the correct -1.855129...
        Vector3D u1 = new Vector3D(-1321008684645961.0 /  268435456.0,
                                   -5774608829631843.0 /  268435456.0,
                                   -7645843051051357.0 / 8589934592.0);
        Vector3D u2 = new Vector3D(-5712344449280879.0 /    2097152.0,
                                   -4550117129121957.0 /    2097152.0,
                                    8846951984510141.0 /     131072.0);
        double sNaive = u1.getX() * u2.getX() + u1.getY() * u2.getY() + u1.getZ() * u2.getZ();
        double sAccurate = u1.dotProduct(u2);
        assertEquals(0.0, sNaive, 1.0e-30);
        assertEquals(-2088690039198397.0 / 1125899906842624.0, sAccurate, 1.0e-15);
    }

    @Test
    void testDotProduct() {
        // we compare accurate versus naive dot product implementations
        // on regular vectors (i.e. not extreme cases like in the previous test)
        Well1024a random = new Well1024a(553267312521321234l);
        for (int i = 0; i < 10000; ++i) {
            double ux = 10000 * random.nextDouble();
            double uy = 10000 * random.nextDouble();
            double uz = 10000 * random.nextDouble();
            double vx = 10000 * random.nextDouble();
            double vy = 10000 * random.nextDouble();
            double vz = 10000 * random.nextDouble();
            double sNaive = ux * vx + uy * vy + uz * vz;
            double sAccurate = new Vector3D(ux, uy, uz).dotProduct(new Vector3D(vx, vy, vz));
            assertEquals(sNaive, sAccurate, 2.5e-16 * sAccurate);
        }
    }

    @Test
    void testAccurateCrossProduct() {
        // the vectors u1 and u2 are nearly but not exactly anti-parallel
        // (7.31e-16 degrees from 180 degrees) naive cross product (i.e.
        // computing u1.x * u2.x + u1.y * u2.y + u1.z * u2.z
        // leads to a result of   [0.0009765, -0.0001220, -0.0039062],
        // instead of the correct [0.0006913, -0.0001254, -0.0007909]
        final Vector3D u1 = new Vector3D(-1321008684645961.0 /   268435456.0,
                                         -5774608829631843.0 /   268435456.0,
                                         -7645843051051357.0 /  8589934592.0);
        final Vector3D u2 = new Vector3D( 1796571811118507.0 /  2147483648.0,
                                          7853468008299307.0 /  2147483648.0,
                                          2599586637357461.0 / 17179869184.0);
        final Vector3D u3 = new Vector3D(12753243807587107.0 / 18446744073709551616.0,
                                         -2313766922703915.0 / 18446744073709551616.0,
                                          -227970081415313.0 /   288230376151711744.0);
        Vector3D cNaive = new Vector3D(u1.getY() * u2.getZ() - u1.getZ() * u2.getY(),
                                       u1.getZ() * u2.getX() - u1.getX() * u2.getZ(),
                                       u1.getX() * u2.getY() - u1.getY() * u2.getX());
        Vector3D cAccurate = u1.crossProduct(u2);
        assertTrue(u3.distance(cNaive) > 2.9 * u3.getNorm());
        assertEquals(0.0, u3.distance(cAccurate), 1.0e-30 * cAccurate.getNorm());
    }

    @Test
    void testCrossProduct() {
        // we compare accurate versus naive cross product implementations
        // on regular vectors (i.e. not extreme cases like in the previous test)
        Well1024a random = new Well1024a(885362227452043214l);
        for (int i = 0; i < 10000; ++i) {
            double ux = 10000 * random.nextDouble();
            double uy = 10000 * random.nextDouble();
            double uz = 10000 * random.nextDouble();
            double vx = 10000 * random.nextDouble();
            double vy = 10000 * random.nextDouble();
            double vz = 10000 * random.nextDouble();
            Vector3D cNaive = new Vector3D(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx);
            Vector3D cAccurate = new Vector3D(ux, uy, uz).crossProduct(new Vector3D(vx, vy, vz));
            assertEquals(0.0, cAccurate.distance(cNaive), 6.0e-15 * cAccurate.getNorm());
        }
    }

    @Test
    void testArithmeticBlending() {

        // Given
        final Vector3D v1 = new Vector3D(1,2,3);
        final Vector3D v2 = new Vector3D(4,5,6);

        final double blendingValue = 0.7;

        // When
        final Vector3D blendedVector = v1.blendArithmeticallyWith(v2, blendingValue);

        // Then
        checkVector(blendedVector, 3.1, 4.1, 5.1);
    }

    private void checkVector(Vector3D v, double x, double y, double z) {
        assertEquals(x, v.getX(), 1.0e-12);
        assertEquals(y, v.getY(), 1.0e-12);
        assertEquals(z, v.getZ(), 1.0e-12);
    }
}
