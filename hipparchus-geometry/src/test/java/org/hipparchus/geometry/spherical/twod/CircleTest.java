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
import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.RotationConvention;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.geometry.partitioning.Transform;
import org.hipparchus.geometry.spherical.oned.Arc;
import org.hipparchus.geometry.spherical.oned.LimitAngle;
import org.hipparchus.geometry.spherical.oned.S1Point;
import org.hipparchus.geometry.spherical.oned.Sphere1D;
import org.hipparchus.geometry.spherical.oned.SubLimitAngle;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.UnitSphereRandomVectorGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class CircleTest {

    @Test
    void testEquator() {
        Circle circle = new Circle(new Vector3D(0, 0, 1000), 1.0e-10).copySelf();
        assertEquals(Vector3D.PLUS_K, circle.getPole());
        assertEquals(1.0e-10, circle.getTolerance(), 1.0e-20);
        circle.revertSelf();
        assertEquals(Vector3D.MINUS_K, circle.getPole());
        assertEquals(Vector3D.PLUS_K, circle.getReverse().getPole());
        assertEquals(Vector3D.MINUS_K, circle.getPole());
    }

    @Test
    void testXY() {
        Circle circle = new Circle(new S2Point(1.2, 2.5), new S2Point(-4.3, 0), 1.0e-10);
        assertEquals(0.0, circle.getPointAt(0).distance(circle.getXAxis()), 1.0e-10);
        assertEquals(0.0, circle.getPointAt(MathUtils.SEMI_PI).distance(circle.getYAxis()), 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, Vector3D.angle(circle.getXAxis(), circle.getYAxis()), 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, Vector3D.angle(circle.getXAxis(), circle.getPole()), 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, Vector3D.angle(circle.getPole(), circle.getYAxis()), 1.0e-10);
        assertEquals(0.0,
                            circle.getPole().distance(Vector3D.crossProduct(circle.getXAxis(), circle.getYAxis())),
                            1.0e-10);
    }

    @Test
    void testReverse() {
        Circle circle = new Circle(new S2Point(1.2, 2.5), new S2Point(-4.3, 0), 1.0e-10);
        Circle reversed = circle.getReverse();
        assertEquals(0.0, reversed.getPointAt(0).distance(reversed.getXAxis()), 1.0e-10);
        assertEquals(0.0, reversed.getPointAt(MathUtils.SEMI_PI).distance(reversed.getYAxis()), 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, Vector3D.angle(reversed.getXAxis(), reversed.getYAxis()), 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, Vector3D.angle(reversed.getXAxis(), reversed.getPole()), 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, Vector3D.angle(reversed.getPole(), reversed.getYAxis()), 1.0e-10);
        assertEquals(0.0,
                            reversed.getPole().distance(Vector3D.crossProduct(reversed.getXAxis(), reversed.getYAxis())),
                            1.0e-10);

        assertEquals(0, Vector3D.angle(circle.getXAxis(), reversed.getXAxis()), 1.0e-10);
        assertEquals(FastMath.PI, Vector3D.angle(circle.getYAxis(), reversed.getYAxis()), 1.0e-10);
        assertEquals(FastMath.PI, Vector3D.angle(circle.getPole(), reversed.getPole()), 1.0e-10);

        assertTrue(circle.sameOrientationAs(circle));
        assertFalse(circle.sameOrientationAs(reversed));

    }

    @Test
    void testPhase() {
        Circle circle = new Circle(new S2Point(1.2, 2.5), new S2Point(-4.3, 0), 1.0e-10);
        Vector3D p = new Vector3D(1, 2, -4);
        Vector3D samePhase = circle.getPointAt(circle.getPhase(p));
        assertEquals(0.0,
                            Vector3D.angle(Vector3D.crossProduct(circle.getPole(), p),
                                           Vector3D.crossProduct(circle.getPole(), samePhase)),
                            1.0e-10);
        assertEquals(MathUtils.SEMI_PI, Vector3D.angle(circle.getPole(), samePhase), 1.0e-10);
        assertEquals(circle.getPhase(p), circle.getPhase(samePhase), 1.0e-10);
        assertEquals(0.0, circle.getPhase(circle.getXAxis()), 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, circle.getPhase(circle.getYAxis()), 1.0e-10);

    }

    @Test
    void testSubSpace() {
        Circle circle = new Circle(new S2Point(1.2, 2.5), new S2Point(-4.3, 0), 1.0e-10);
        assertEquals(0.0, circle.toSubSpace(new S2Point(circle.getXAxis())).getAlpha(), 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, circle.toSubSpace(new S2Point(circle.getYAxis())).getAlpha(), 1.0e-10);
        Vector3D p = new Vector3D(1, 2, -4);
        assertEquals(circle.getPhase(p), circle.toSubSpace(new S2Point(p)).getAlpha(), 1.0e-10);
    }

    @Test
    void testSpace() {
        Circle circle = new Circle(new S2Point(1.2, 2.5), new S2Point(-4.3, 0), 1.0e-10);
        for (double alpha = 0; alpha < MathUtils.TWO_PI; alpha += 0.1) {
            Vector3D p = new Vector3D(FastMath.cos(alpha), circle.getXAxis(),
                                      FastMath.sin(alpha), circle.getYAxis());
            Vector3D q = circle.toSpace(new S1Point(alpha)).getVector();
            assertEquals(0.0, p.distance(q), 1.0e-10);
            assertEquals(MathUtils.SEMI_PI, Vector3D.angle(circle.getPole(), q), 1.0e-10);
        }
    }

    @Test
    void testOffset() {
        Circle circle = new Circle(Vector3D.PLUS_K, 1.0e-10);
        assertEquals(0.0,                circle.getOffset(new S2Point(Vector3D.PLUS_I)),  1.0e-10);
        assertEquals(0.0,                circle.getOffset(new S2Point(Vector3D.MINUS_I)), 1.0e-10);
        assertEquals(0.0,                circle.getOffset(new S2Point(Vector3D.PLUS_J)),  1.0e-10);
        assertEquals(0.0,                circle.getOffset(new S2Point(Vector3D.MINUS_J)), 1.0e-10);
        assertEquals(-MathUtils.SEMI_PI, circle.getOffset(new S2Point(Vector3D.PLUS_K)),  1.0e-10);
        assertEquals( MathUtils.SEMI_PI, circle.getOffset(new S2Point(Vector3D.MINUS_K)), 1.0e-10);

    }

    @Test
    void testInsideArc() {
        RandomGenerator random = new Well1024a(0xbfd34e92231bbcfeL);
        UnitSphereRandomVectorGenerator sphRandom = new UnitSphereRandomVectorGenerator(3, random);
        for (int i = 0; i < 100; ++i) {
            Circle c1 = new Circle(new Vector3D(sphRandom.nextVector()), 1.0e-10);
            Circle c2 = new Circle(new Vector3D(sphRandom.nextVector()), 1.0e-10);
            checkArcIsInside(c1, c2);
            checkArcIsInside(c2, c1);
        }
    }

    private void checkArcIsInside(final Circle arcCircle, final Circle otherCircle) {
        Arc arc = arcCircle.getInsideArc(otherCircle);
        assertEquals(FastMath.PI, arc.getSize(), 1.0e-10);
        for (double alpha = arc.getInf(); alpha < arc.getSup(); alpha += 0.1) {
            assertTrue(otherCircle.getOffset(arcCircle.getPointAt(alpha)) <= 2.0e-15);
        }
        for (double alpha = arc.getSup(); alpha < arc.getInf() + MathUtils.TWO_PI; alpha += 0.1) {
            assertTrue(otherCircle.getOffset(arcCircle.getPointAt(alpha)) >= -2.0e-15);
        }
    }

    @Test
    void testTransform() {
        RandomGenerator random = new Well1024a(0x16992fc4294bf2f1L);
        UnitSphereRandomVectorGenerator sphRandom = new UnitSphereRandomVectorGenerator(3, random);
        for (int i = 0; i < 100; ++i) {

            Rotation r = new Rotation(new Vector3D(sphRandom.nextVector()),
                                      FastMath.PI * random.nextDouble(),
                                      RotationConvention.VECTOR_OPERATOR);
            Transform<Sphere2D, S2Point, Circle, SubCircle, Sphere1D, S1Point, LimitAngle, SubLimitAngle> t = Circle.getTransform(r);

            S2Point  p = new S2Point(new Vector3D(sphRandom.nextVector()));
            S2Point tp = t.apply(p);
            assertEquals(0.0, r.applyTo(p.getVector()).distance(tp.getVector()), 1.0e-10);

            Circle  c = new Circle(new Vector3D(sphRandom.nextVector()), 1.0e-10);
            Circle tc = t.apply(c);
            assertEquals(0.0, r.applyTo(c.getPole()).distance(tc.getPole()),   1.0e-10);
            assertEquals(0.0, r.applyTo(c.getXAxis()).distance(tc.getXAxis()), 1.0e-10);
            assertEquals(0.0, r.applyTo(c.getYAxis()).distance(tc.getYAxis()), 1.0e-10);
            assertEquals(c.getTolerance(), t.apply(c).getTolerance(), 1.0e-10);

            SubLimitAngle  sub = new LimitAngle(new S1Point(MathUtils.TWO_PI * random.nextDouble()),
                                                random.nextBoolean(), 1.0e-10).wholeHyperplane();
            Vector3D psub = c.getPointAt(sub.getHyperplane().getLocation().getAlpha());
            SubLimitAngle tsub = t.apply(sub, c, tc);
            Vector3D ptsub = tc.getPointAt(tsub.getHyperplane().getLocation().getAlpha());
            assertEquals(0.0, r.applyTo(psub).distance(ptsub), 1.0e-10);

        }
    }

    @Test
    void testTooSmallTolerance() {
        assertThrows(MathIllegalArgumentException.class, () -> new Circle(Vector3D.PLUS_K, 0.9 * Sphere2D.SMALLEST_TOLERANCE));
    }

    /** Check {@link Circle#getArc(S2Point, S2Point)}. */
    @Test
    void testGetArc() {
        // setup
        double tol = 1e-6;
        Circle circle = new Circle(Vector3D.PLUS_K, tol);

        // action
        Arc arc = circle.getArc(new S2Point(circle.getPointAt(0)), new S2Point(circle.getPointAt(1)));
        // verify
        assertEquals(0.5, arc.getBarycenter(), tol);
        assertEquals(0, arc.getInf(), tol);
        assertEquals(1, arc.getSup(), tol);
        assertEquals(arc.getTolerance(), tol, 0);

        // action, crossing discontinuity
        arc = circle.getArc(new S2Point(circle.getPointAt(3)), new S2Point(circle.getPointAt(-3)));
        // verify
        assertEquals(FastMath.PI, arc.getBarycenter(), tol);
        assertEquals(3, arc.getInf(), tol);
        assertEquals(arc.getSup(), 2 * FastMath.PI - 3, tol);
        assertEquals(arc.getTolerance(), tol, 0);
    }

}
