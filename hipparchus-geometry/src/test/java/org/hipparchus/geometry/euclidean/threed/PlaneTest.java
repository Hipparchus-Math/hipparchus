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

import org.hipparchus.exception.MathRuntimeException;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

class PlaneTest {

    @Test
    void testContains() throws MathRuntimeException {
        Plane p = new Plane(new Vector3D(0, 0, 1), new Vector3D(0, 0, 1), 1.0e-10);
        assertTrue(p.contains(new Vector3D(0, 0, 1)));
        assertTrue(p.contains(new Vector3D(17, -32, 1)));
        assertFalse(p.contains(new Vector3D(17, -32, 1.001)));
    }

    @Test
    void testOffset() throws MathRuntimeException {
        Vector3D p1 = new Vector3D(1, 1, 1);
        Plane p = new Plane(p1, new Vector3D(0.2, 0, 0), 1.0e-10);
        assertEquals(-5.0, p.getOffset(new Vector3D(-4, 0, 0)), 1.0e-10);
        assertEquals(+5.0, p.getOffset(new Vector3D(6, 10, -12)), 1.0e-10);
        assertEquals(0.3,
                            p.getOffset(new Vector3D(1.0, p1, 0.3, p.getNormal())),
                            1.0e-10);
        assertEquals(-0.3,
                            p.getOffset(new Vector3D(1.0, p1, -0.3, p.getNormal())),
                            1.0e-10);
    }

    @Test
    void testPoint() throws MathRuntimeException {
        Plane p = new Plane(new Vector3D(2, -3, 1), new Vector3D(1, 4, 9), 1.0e-10);
        assertTrue(p.contains(p.getOrigin()));
    }

    @Test
    void testThreePoints() throws MathRuntimeException {
        Vector3D p1 = new Vector3D(1.2, 3.4, -5.8);
        Vector3D p2 = new Vector3D(3.4, -5.8, 1.2);
        Vector3D p3 = new Vector3D(-2.0, 4.3, 0.7);
        Plane    p  = new Plane(p1, p2, p3, 1.0e-10);
        assertTrue(p.contains(p1));
        assertTrue(p.contains(p2));
        assertTrue(p.contains(p3));
    }

    @Test
    void testRotate() throws MathRuntimeException {
        Vector3D p1 = new Vector3D(1.2, 3.4, -5.8);
        Vector3D p2 = new Vector3D(3.4, -5.8, 1.2);
        Vector3D p3 = new Vector3D(-2.0, 4.3, 0.7);
        Plane    p  = new Plane(p1, p2, p3, 1.0e-10);
        Vector3D oldNormal = p.getNormal();

        p = p.rotate(p2, new Rotation(p2.subtract(p1), 1.7, RotationConvention.VECTOR_OPERATOR));
        assertTrue(p.contains(p1));
        assertTrue(p.contains(p2));
        assertFalse(p.contains(p3));

        p = p.rotate(p2, new Rotation(oldNormal, 0.1, RotationConvention.VECTOR_OPERATOR));
        assertFalse(p.contains(p1));
        assertTrue(p.contains(p2));
        assertFalse(p.contains(p3));

        p = p.rotate(p1, new Rotation(oldNormal, 0.1, RotationConvention.VECTOR_OPERATOR));
        assertFalse(p.contains(p1));
        assertFalse(p.contains(p2));
        assertFalse(p.contains(p3));

    }

    @Test
    void testTranslate() throws MathRuntimeException {
        Vector3D p1 = new Vector3D(1.2, 3.4, -5.8);
        Vector3D p2 = new Vector3D(3.4, -5.8, 1.2);
        Vector3D p3 = new Vector3D(-2.0, 4.3, 0.7);
        Plane    p  = new Plane(p1, p2, p3, 1.0e-10);

        p = p.translate(new Vector3D(2.0, p.getU(), -1.5, p.getV()));
        assertTrue(p.contains(p1));
        assertTrue(p.contains(p2));
        assertTrue(p.contains(p3));

        p = p.translate(new Vector3D(-1.2, p.getNormal()));
        assertFalse(p.contains(p1));
        assertFalse(p.contains(p2));
        assertFalse(p.contains(p3));

        p = p.translate(new Vector3D(+1.2, p.getNormal()));
        assertTrue(p.contains(p1));
        assertTrue(p.contains(p2));
        assertTrue(p.contains(p3));

    }

    @Test
    void testIntersection() throws MathRuntimeException {
        Plane p = new Plane(new Vector3D(1, 2, 3), new Vector3D(-4, 1, -5), 1.0e-10);
        Line  l = new Line(new Vector3D(0.2, -3.5, 0.7), new Vector3D(1.2, -2.5, -0.3), 1.0e-10);
        Vector3D point = p.intersection(l);
        assertTrue(p.contains(point));
        assertTrue(l.contains(point));
        assertNull(p.intersection(new Line(new Vector3D(10, 10, 10),
                                                  new Vector3D(10, 10, 10).add(p.getNormal().orthogonal()),
                                                  1.0e-10)));
    }

    @Test
    void testIntersection2() throws MathRuntimeException {
        Vector3D p1  = new Vector3D (1.2, 3.4, -5.8);
        Vector3D p2  = new Vector3D (3.4, -5.8, 1.2);
        Plane    pA  = new Plane(p1, p2, new Vector3D (-2.0, 4.3, 0.7), 1.0e-10);
        Plane    pB  = new Plane(p1, new Vector3D (11.4, -3.8, 5.1), p2, 1.0e-10);
        Line     l   = pA.intersection(pB);
        assertTrue(l.contains(p1));
        assertTrue(l.contains(p2));
        assertNull(pA.intersection(pA));
    }

    @Test
    void testIntersection3() throws MathRuntimeException {
        Vector3D reference = new Vector3D (1.2, 3.4, -5.8);
        Plane p1 = new Plane(reference, new Vector3D(1, 3, 3), 1.0e-10);
        Plane p2 = new Plane(reference, new Vector3D(-2, 4, 0), 1.0e-10);
        Plane p3 = new Plane(reference, new Vector3D(7, 0, -4), 1.0e-10);
        Vector3D p = Plane.intersection(p1, p2, p3);
        assertEquals(reference.getX(), p.getX(), 1.0e-10);
        assertEquals(reference.getY(), p.getY(), 1.0e-10);
        assertEquals(reference.getZ(), p.getZ(), 1.0e-10);
    }

    @Test
    void testSimilar() throws MathRuntimeException {
        Vector3D p1  = new Vector3D (1.2, 3.4, -5.8);
        Vector3D p2  = new Vector3D (3.4, -5.8, 1.2);
        Vector3D p3  = new Vector3D (-2.0, 4.3, 0.7);
        Plane    pA  = new Plane(p1, p2, p3, 1.0e-10);
        Plane    pB  = new Plane(p1, new Vector3D (11.4, -3.8, 5.1), p2, 1.0e-10);
        assertFalse(pA.isSimilarTo(pB));
        assertTrue(pA.isSimilarTo(pA));
        assertTrue(pA.isSimilarTo(new Plane(p1, p3, p2, 1.0e-10)));
        Vector3D shift = new Vector3D(0.3, pA.getNormal());
        assertFalse(pA.isSimilarTo(new Plane(p1.add(shift),
            p3.add(shift),
            p2.add(shift),
            1.0e-10)));
        assertTrue(pA.sameOrientationAs(pB));
        Plane pC = pB.copySelf();
        pC.revertSelf();
        assertFalse(pA.sameOrientationAs(pC));
    }

    @Test
    public void testEmptyHyperplane() {
        Vector3D p1    = new Vector3D (1.2, 3.4, -5.8);
        Vector3D p2    = new Vector3D (3.4, -5.8, 1.2);
        Vector3D p3    = new Vector3D (-2.0, 4.3, 0.7);
        Plane    plane = new Plane(p1, p2, p3, 1.0e-10);
        SubPlane sp = plane.emptyHyperplane();
        assertTrue(sp.isEmpty());
    }

    @Test
    public void testMove() {
        Vector3D p1    = new Vector3D (1.2, 3.4, -5.8);
        Vector3D p2    = new Vector3D (3.4, -5.8, 1.2);
        Vector3D p3    = new Vector3D (-2.0, 4.3, 0.7);
        Plane    plane = new Plane(p1, p2, p3, 1.0e-10);
        assertEquals(0.0, plane.getOffset(plane.arbitraryPoint()), 1.0e-10);
        assertEquals(0.0, plane.getOffset(p1), 1.0e-10);
        assertEquals(0.0, plane.getOffset(p2), 1.0e-10);
        assertEquals(0.0, plane.getOffset(p3), 1.0e-10);

        assertEquals( 1.5, plane.getOffset(plane.moveToOffset(p1,  1.5)), 1.0e-10);
        assertEquals(-1.5, plane.getOffset(plane.moveToOffset(p2, -1.5)), 1.0e-10);
        assertEquals(17.0, plane.getOffset(plane.moveToOffset(p3, 17.0)), 1.0e-10);

    }

}
