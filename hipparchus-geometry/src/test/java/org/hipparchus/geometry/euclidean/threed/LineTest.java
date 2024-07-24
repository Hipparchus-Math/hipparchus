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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

class LineTest {

    @Test
    void testContains() throws MathIllegalArgumentException, MathRuntimeException {
        Vector3D p1 = new Vector3D(0, 0, 1);
        Line l = new Line(p1, new Vector3D(0, 0, 2), 1.0e-10);
        assertTrue(l.contains(p1));
        assertTrue(l.contains(new Vector3D(1.0, p1, 0.3, l.getDirection())));
        Vector3D u = l.getDirection().orthogonal();
        Vector3D v = Vector3D.crossProduct(l.getDirection(), u);
        for (double alpha = 0; alpha < 2 * FastMath.PI; alpha += 0.3) {
            assertFalse(l.contains(p1.add(new Vector3D(FastMath.cos(alpha), u,
                FastMath.sin(alpha), v))));
        }
    }

    @Test
    void testSimilar() throws MathIllegalArgumentException, MathRuntimeException {
        Vector3D p1  = new Vector3D (1.2, 3.4, -5.8);
        Vector3D p2  = new Vector3D (3.4, -5.8, 1.2);
        Line     lA  = new Line(p1, p2, 1.0e-10);
        Line     lB  = new Line(p2, p1, 1.0e-10);
        assertTrue(lA.isSimilarTo(lB));
        assertFalse(lA.isSimilarTo(new Line(p1, p1.add(lA.getDirection().orthogonal()), 1.0e-10)));
    }

    @Test
    void testPointDistance() throws MathIllegalArgumentException {
        Line l = new Line(new Vector3D(0, 1, 1), new Vector3D(0, 2, 2), 1.0e-10);
        assertEquals(FastMath.sqrt(3.0 / 2.0), l.distance(new Vector3D(1, 0, 1)), 1.0e-10);
        assertEquals(0, l.distance(new Vector3D(0, -4, -4)), 1.0e-10);
    }

    @Test
    void testLineDistance() throws MathIllegalArgumentException {
        Line l = new Line(new Vector3D(0, 1, 1), new Vector3D(0, 2, 2), 1.0e-10);
        assertEquals(1.0,
                            l.distance(new Line(new Vector3D(1, 0, 1), new Vector3D(1, 0, 2), 1.0e-10)),
                            1.0e-10);
        assertEquals(0.5,
                            l.distance(new Line(new Vector3D(-0.5, 0, 0), new Vector3D(-0.5, -1, -1), 1.0e-10)),
                            1.0e-10);
        assertEquals(0.0,
                            l.distance(l),
                            1.0e-10);
        assertEquals(0.0,
                            l.distance(new Line(new Vector3D(0, -4, -4), new Vector3D(0, -5, -5), 1.0e-10)),
                            1.0e-10);
        assertEquals(0.0,
                            l.distance(new Line(new Vector3D(0, -4, -4), new Vector3D(0, -3, -4), 1.0e-10)),
                            1.0e-10);
        assertEquals(0.0,
                            l.distance(new Line(new Vector3D(0, -4, -4), new Vector3D(1, -4, -4), 1.0e-10)),
                            1.0e-10);
        assertEquals(FastMath.sqrt(8),
                            l.distance(new Line(new Vector3D(0, -4, 0), new Vector3D(1, -4, 0), 1.0e-10)),
                            1.0e-10);
    }

    @Test
    void testClosest() throws MathIllegalArgumentException {
        Line l = new Line(new Vector3D(0, 1, 1), new Vector3D(0, 2, 2), 1.0e-10);
        assertEquals(0.0,
                            l.closestPoint(new Line(new Vector3D(1, 0, 1), new Vector3D(1, 0, 2), 1.0e-10)).distance(new Vector3D(0, 0, 0)),
                            1.0e-10);
        assertEquals(0.5,
                            l.closestPoint(new Line(new Vector3D(-0.5, 0, 0), new Vector3D(-0.5, -1, -1), 1.0e-10)).distance(new Vector3D(-0.5, 0, 0)),
                            1.0e-10);
        assertEquals(0.0,
                            l.closestPoint(l).distance(new Vector3D(0, 0, 0)),
                            1.0e-10);
        assertEquals(0.0,
                            l.closestPoint(new Line(new Vector3D(0, -4, -4), new Vector3D(0, -5, -5), 1.0e-10)).distance(new Vector3D(0, 0, 0)),
                            1.0e-10);
        assertEquals(0.0,
                            l.closestPoint(new Line(new Vector3D(0, -4, -4), new Vector3D(0, -3, -4), 1.0e-10)).distance(new Vector3D(0, -4, -4)),
                            1.0e-10);
        assertEquals(0.0,
                            l.closestPoint(new Line(new Vector3D(0, -4, -4), new Vector3D(1, -4, -4), 1.0e-10)).distance(new Vector3D(0, -4, -4)),
                            1.0e-10);
        assertEquals(0.0,
                            l.closestPoint(new Line(new Vector3D(0, -4, 0), new Vector3D(1, -4, 0), 1.0e-10)).distance(new Vector3D(0, -2, -2)),
                            1.0e-10);
    }

    @Test
    void testIntersection() throws MathIllegalArgumentException {
        Line l = new Line(new Vector3D(0, 1, 1), new Vector3D(0, 2, 2), 1.0e-10);
        assertNull(l.intersection(new Line(new Vector3D(1, 0, 1), new Vector3D(1, 0, 2), 1.0e-10)));
        assertNull(l.intersection(new Line(new Vector3D(-0.5, 0, 0), new Vector3D(-0.5, -1, -1), 1.0e-10)));
        assertEquals(0.0,
                            l.intersection(l).distance(new Vector3D(0, 0, 0)),
                            1.0e-10);
        assertEquals(0.0,
                            l.intersection(new Line(new Vector3D(0, -4, -4), new Vector3D(0, -5, -5), 1.0e-10)).distance(new Vector3D(0, 0, 0)),
                            1.0e-10);
        assertEquals(0.0,
                            l.intersection(new Line(new Vector3D(0, -4, -4), new Vector3D(0, -3, -4), 1.0e-10)).distance(new Vector3D(0, -4, -4)),
                            1.0e-10);
        assertEquals(0.0,
                            l.intersection(new Line(new Vector3D(0, -4, -4), new Vector3D(1, -4, -4), 1.0e-10)).distance(new Vector3D(0, -4, -4)),
                            1.0e-10);
        assertNull(l.intersection(new Line(new Vector3D(0, -4, 0), new Vector3D(1, -4, 0), 1.0e-10)));
    }

    @Test
    void testRevert() {

        // setup
        Line line = new Line(new Vector3D(1653345.6696423641, 6170370.041579291, 90000),
                             new Vector3D(1650757.5050732433, 6160710.879908984, 0.9),
                             1.0e-10);
        Vector3D expected = line.getDirection().negate();

        // action
        Line reverted = line.revert();

        // verify
        assertArrayEquals(expected.toArray(), reverted.getDirection().toArray(), 0);

    }

    /** Test for issue #78. */
    @Test
    void testCancellation() {
         // setup
        // point
        Vector3D p = new Vector3D(1e16, 1e16, 1e16);
        // unit vector
        Vector3D u = new Vector3D(-1, -1, -1);

        // action
        Line line = Line.fromDirection(p, u, 0);

        // verify
        double ulp = FastMath.ulp(p.getNorm());
        assertArrayEquals(line.getOrigin().toArray(), Vector3D.ZERO.toArray(), ulp);
        assertArrayEquals(line.getDirection().toArray(), u.normalize().toArray(), 0);
        assertEquals(0, line.getTolerance(), 0);
    }

}
