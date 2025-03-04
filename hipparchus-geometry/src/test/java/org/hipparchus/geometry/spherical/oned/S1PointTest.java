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
package org.hipparchus.geometry.spherical.oned;

import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class S1PointTest {

    @Test
    void testS1Point() {
        for (int k = -2; k < 3; ++k) {
            S1Point p = new S1Point(1.0 + k * MathUtils.TWO_PI);
            assertEquals(FastMath.cos(1.0), p.getVector().getX(), 1.0e-10);
            assertEquals(FastMath.sin(1.0), p.getVector().getY(), 1.0e-10);
            assertFalse(p.isNaN());
        }
    }

    @Test
    void testNaN() {
        assertTrue(S1Point.NaN.isNaN());
        assertEquals(S1Point.NaN, new S1Point(Double.NaN));
        assertNotEquals(S1Point.NaN, new S1Point(1.0));
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testEquals() {
        S1Point a = new S1Point(1.0);
        S1Point b = new S1Point(1.0);
        assertEquals(a.hashCode(), b.hashCode());
        assertNotSame(a, b);
        assertEquals(a, b);
        assertEquals(a, a);
        assertNotEquals('a', a);
        assertEquals(S1Point.NaN, S1Point.NaN);
        assertEquals(S1Point.NaN, new S1Point(Double.NaN));
    }

    @Test
    void testEqualsIeee754() {
        S1Point a = new S1Point(1.0);
        S1Point b = new S1Point(1.0);
        assertEquals(a.hashCode(), b.hashCode());
        assertNotSame(a, b);
        assertTrue(a.equalsIeee754(b));
        assertTrue(a.equalsIeee754(a));
        assertFalse(a.equalsIeee754('a'));
        assertFalse(S1Point.NaN.equalsIeee754(S1Point.NaN));
        assertFalse(S1Point.NaN.equalsIeee754(new S1Point(Double.NaN)));
    }

    @Test
    void testDistance() {
        S1Point a = new S1Point(1.0);
        S1Point b = new S1Point(a.getAlpha() + MathUtils.SEMI_PI);
        assertEquals(MathUtils.SEMI_PI, a.distance(b), 1.0e-10);
    }

    @Test
    void testSpace() {
        S1Point a = new S1Point(1.0);
        assertInstanceOf(Sphere1D.class, a.getSpace());
        assertEquals(1, a.getSpace().getDimension());
        try {
            a.getSpace().getSubSpace();
            fail("an exception should have been thrown");
        } catch (MathRuntimeException muoe) {
            // expected
        }
    }

    @Test
    void testMoveTowards() {
        final S1Point s1 = new S1Point(2.0);
        final S1Point s2 = new S1Point(4.0);
        for (double r = 0.0; r <= 1.0; r += FastMath.scalb(1.0, -10)) {
            // motion should be linear according to angles
            assertEquals(r * s1.distance(s2), s1.distance(s1.moveTowards(s2, r)), 1.0e-14);
        }
    }

}
