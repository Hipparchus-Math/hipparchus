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
package org.hipparchus.geometry.spherical.twod;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.geometry.spherical.oned.Sphere1D;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class S2PointTest {


    @Test
    void testS2Point() {
        for (int k = -2; k < 3; ++k) {
            S2Point p = new S2Point(1.0 + k * MathUtils.TWO_PI, 1.4);
            assertEquals(1.0 + k * MathUtils.TWO_PI, p.getTheta(), 1.0e-10);
            assertEquals(1.4, p.getPhi(), 1.0e-10);
            assertEquals(FastMath.cos(1.0) * FastMath.sin(1.4), p.getVector().getX(), 1.0e-10);
            assertEquals(FastMath.sin(1.0) * FastMath.sin(1.4), p.getVector().getY(), 1.0e-10);
            assertEquals(FastMath.cos(1.4), p.getVector().getZ(), 1.0e-10);
            assertFalse(p.isNaN());
        }
    }

    @Test
    void testNegativePolarAngle() {
        assertThrows(MathIllegalArgumentException.class, () -> new S2Point(1.0, -1.0));
    }

    @Test
    void testTooLargePolarAngle() {
        assertThrows(MathIllegalArgumentException.class, () -> new S2Point(1.0, 3.5));
    }

    @Test
    void testNaN() {
        assertTrue(S2Point.NaN.isNaN());
        assertEquals(S2Point.NaN, new S2Point(Double.NaN, 1.0));
        assertNotEquals(S2Point.NaN, new S2Point(1.0, 1.3));
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testEquals() {
        S2Point a = new S2Point(1.0, 1.0);
        S2Point b = new S2Point(1.0, 1.0);
        assertEquals(a.hashCode(), b.hashCode());
        assertNotSame(a, b);
        assertEquals(a, b);
        assertEquals(a, a);
        assertNotEquals('a', a);
        assertEquals(S2Point.NaN, S2Point.NaN);
        assertEquals(S2Point.NaN, new S2Point(Double.NaN, 0.0));
        assertEquals(S2Point.NaN, new S2Point(0.0, Double.NaN));
    }

    @Test
    void testEqualsIeee754() {
        S2Point a = new S2Point(1.0, 1.0);
        S2Point b = new S2Point(1.0, 1.0);
        assertEquals(a.hashCode(), b.hashCode());
        assertNotSame(a, b);
        assertTrue(a.equalsIeee754(b));
        assertTrue(a.equalsIeee754(a));
        assertFalse(a.equalsIeee754('a'));
        assertFalse(S2Point.NaN.equalsIeee754(S2Point.NaN));
        assertFalse(S2Point.NaN.equalsIeee754(new S2Point(Double.NaN, 0.0)));
        assertFalse(S2Point.NaN.equalsIeee754(new S2Point(0.0, Double.NaN)));
    }

    @Test
    void testDistance() {
        S2Point a = new S2Point(1.0, MathUtils.SEMI_PI);
        S2Point b = new S2Point(a.getTheta() + MathUtils.SEMI_PI, a.getPhi());
        assertEquals(MathUtils.SEMI_PI, a.distance(b), 1.0e-10);
        assertEquals(FastMath.PI, a.distance(a.negate()), 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, S2Point.MINUS_I.distance(S2Point.MINUS_K), 1.0e-10);
        assertEquals(0.0, new S2Point(1.0, 0).distance(new S2Point(2.0, 0)), 1.0e-10);
    }

    @Test
    void testNegate() {
        RandomGenerator generator = new Well1024a(0x79d1bc2e0999d238L);
        for (int i = 0; i < 100000; ++i) {
            S2Point p = new S2Point(MathUtils.TWO_PI * generator.nextDouble(),
                                    FastMath.PI * generator.nextDouble());
            S2Point np = new S2Point(p.negate().getTheta(), p.negate().getPhi());
            assertEquals(FastMath.PI, p.distance(np), 1.4e-15);
        }
    }

    @Test
    void testSpace() {
        S2Point a = new S2Point(1.0, 1.0);
        assertInstanceOf(Sphere2D.class, a.getSpace());
        assertEquals(2, a.getSpace().getDimension());
        assertInstanceOf(Sphere1D.class, a.getSpace().getSubSpace());
    }

    @Test
    void testMoveTowards() {
        final S2Point s1 = new S2Point(2.0, 0.5);
        final S2Point s2 = new S2Point(4.0, 2.5);
        for (double r = 0.0; r <= 1.0; r += FastMath.scalb(1.0, -10)) {
            // motion should be linear according to angles
            assertEquals(r * s1.distance(s2), s1.distance(s1.moveTowards(s2, r)), 1.0e-14);
        }
    }

    @Test
    void testMoveTowardsSpecialCase() {
        final S2Point s = new S2Point(2.0, 0.5);
        for (double r = 0.0; r <= 1.0; r += FastMath.scalb(1.0, -10)) {
            // motion should be linear according to angles
            assertEquals(0.0, s.distance(s.moveTowards(s, r)), 1.0e-20);
        }
    }

}
