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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.geometry.LocalizedGeometryFormats;
import org.hipparchus.geometry.enclosing.EnclosingBall;
import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.geometry.partitioning.Region;
import org.hipparchus.geometry.partitioning.Region.Location;
import org.hipparchus.geometry.partitioning.RegionFactory;
import org.hipparchus.geometry.spherical.oned.ArcsSet;
import org.hipparchus.geometry.spherical.oned.LimitAngle;
import org.hipparchus.geometry.spherical.oned.S1Point;
import org.hipparchus.geometry.spherical.oned.Sphere1D;
import org.hipparchus.geometry.spherical.oned.SubLimitAngle;
import org.hipparchus.random.UnitSphereRandomVectorGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;
import java.util.function.IntPredicate;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class SphericalPolygonsSetTest {

    @Test
    void testFullSphere() {
        SphericalPolygonsSet full = new SphericalPolygonsSet(1.0e-10);
        UnitSphereRandomVectorGenerator random =
                new UnitSphereRandomVectorGenerator(3, new Well1024a(0x852fd2a0ed8d2f6dL));
        for (int i = 0; i < 1000; ++i) {
            Vector3D v = new Vector3D(random.nextVector());
            assertEquals(Location.INSIDE, full.checkPoint(new S2Point(v)));
        }
        assertEquals(4 * FastMath.PI, new SphericalPolygonsSet(0.01, new S2Point[0]).getSize(), 1.0e-10);
        assertEquals(0, new SphericalPolygonsSet(0.01, new S2Point[0]).getBoundarySize(), 1.0e-10);
        assertEquals(0, full.getBoundaryLoops().size());
        assertTrue(full.getEnclosingCap().getRadius() > 0);
        assertTrue(Double.isInfinite(full.getEnclosingCap().getRadius()));
        assertEquals(Location.INSIDE, full.checkPoint(full.getInteriorPoint()));
    }

    @Test
    void testEmpty() {
        SphericalPolygonsSet empty =
            (SphericalPolygonsSet) new RegionFactory<Sphere2D, S2Point, Circle, SubCircle>().
                    getComplement(new SphericalPolygonsSet(1.0e-10));
        UnitSphereRandomVectorGenerator random =
                new UnitSphereRandomVectorGenerator(3, new Well1024a(0x76d9205d6167b6ddL));
        for (int i = 0; i < 1000; ++i) {
            Vector3D v = new Vector3D(random.nextVector());
            assertEquals(Location.OUTSIDE, empty.checkPoint(new S2Point(v)));
        }
        assertEquals(0, empty.getSize(), 1.0e-10);
        assertEquals(0, empty.getBoundarySize(), 1.0e-10);
        assertEquals(0, empty.getBoundaryLoops().size());
        assertTrue(empty.getEnclosingCap().getRadius() < 0);
        assertTrue(Double.isInfinite(empty.getEnclosingCap().getRadius()));
        assertNull(empty.getInteriorPoint());
    }

    @Test
    void testSouthHemisphere() {
        double tol = 0.01;
        double sinTol = FastMath.sin(tol);
        SphericalPolygonsSet south = new SphericalPolygonsSet(Vector3D.MINUS_K, tol);
        UnitSphereRandomVectorGenerator random =
                new UnitSphereRandomVectorGenerator(3, new Well1024a(0x6b9d4a6ad90d7b0bL));
        for (int i = 0; i < 1000; ++i) {
            Vector3D v = new Vector3D(random.nextVector());
            if (v.getZ() < -sinTol) {
                assertEquals(Location.INSIDE, south.checkPoint(new S2Point(v)));
            } else if (v.getZ() > sinTol) {
                assertEquals(Location.OUTSIDE, south.checkPoint(new S2Point(v)));
            } else {
                assertEquals(Location.BOUNDARY, south.checkPoint(new S2Point(v)));
            }
        }
        assertEquals(1, south.getBoundaryLoops().size());

        EnclosingBall<Sphere2D, S2Point> southCap = south.getEnclosingCap();
        assertEquals(0.0, S2Point.MINUS_K.distance(southCap.getCenter()), 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, southCap.getRadius(), 1.0e-10);

        EnclosingBall<Sphere2D, S2Point> northCap =
                ((SphericalPolygonsSet) new RegionFactory<Sphere2D, S2Point, Circle, SubCircle>().
                        getComplement(south)).getEnclosingCap();
        assertEquals(0.0, S2Point.PLUS_K.distance(northCap.getCenter()), 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, northCap.getRadius(), 1.0e-10);

        assertEquals(0.0, S2Point.distance(S2Point.MINUS_K, south.getInteriorPoint()), 1.0e-15);
        assertEquals(Location.INSIDE, south.checkPoint(south.getInteriorPoint()));

        final RegionFactory<Sphere2D, S2Point, Circle, SubCircle> factory = new RegionFactory<>();
        final SphericalPolygonsSet reversed = (SphericalPolygonsSet) factory.getComplement(south);
        assertEquals(0.0, S2Point.distance(S2Point.PLUS_K, reversed.getInteriorPoint()), 1.0e-15);
        assertEquals(Location.INSIDE, reversed.checkPoint(reversed.getInteriorPoint()));

    }

    @Test
    void testPositiveOctantByIntersection() {
        double tol = 0.01;
        double sinTol = FastMath.sin(tol);
        RegionFactory<Sphere2D, S2Point, Circle, SubCircle> factory = new RegionFactory<>();
        SphericalPolygonsSet plusX = new SphericalPolygonsSet(Vector3D.PLUS_I, tol);
        SphericalPolygonsSet plusY = new SphericalPolygonsSet(Vector3D.PLUS_J, tol);
        SphericalPolygonsSet plusZ = new SphericalPolygonsSet(Vector3D.PLUS_K, tol);
        SphericalPolygonsSet octant =
                (SphericalPolygonsSet) factory.intersection(factory.intersection(plusX, plusY), plusZ);
        UnitSphereRandomVectorGenerator random =
                new UnitSphereRandomVectorGenerator(3, new Well1024a(0x9c9802fde3cbcf25L));
        for (int i = 0; i < 1000; ++i) {
            Vector3D v = new Vector3D(random.nextVector());
            if ((v.getX() > sinTol) && (v.getY() > sinTol) && (v.getZ() > sinTol)) {
                assertEquals(Location.INSIDE, octant.checkPoint(new S2Point(v)));
            } else if ((v.getX() < -sinTol) || (v.getY() < -sinTol) || (v.getZ() < -sinTol)) {
                assertEquals(Location.OUTSIDE, octant.checkPoint(new S2Point(v)));
            } else {
                assertEquals(Location.BOUNDARY, octant.checkPoint(new S2Point(v)));
            }
        }

        List<Vertex> loops = octant.getBoundaryLoops();
        assertEquals(1, loops.size());
        boolean xPFound = false;
        boolean yPFound = false;
        boolean zPFound = false;
        boolean xVFound = false;
        boolean yVFound = false;
        boolean zVFound = false;
        Vertex first = loops.get(0);
        int count = 0;
        for (Vertex v = first; count == 0 || v != first; v = v.getOutgoing().getEnd()) {
            ++count;
            Edge e = v.getIncoming();
            assertSame(v, e.getStart().getOutgoing().getEnd());
            xPFound = xPFound || e.getCircle().getPole().distance(Vector3D.PLUS_I) < 1.0e-10;
            yPFound = yPFound || e.getCircle().getPole().distance(Vector3D.PLUS_J) < 1.0e-10;
            zPFound = zPFound || e.getCircle().getPole().distance(Vector3D.PLUS_K) < 1.0e-10;
            assertEquals(MathUtils.SEMI_PI, e.getLength(), 1.0e-10);
            xVFound = xVFound || v.getLocation().getVector().distance(Vector3D.PLUS_I) < 1.0e-10;
            yVFound = yVFound || v.getLocation().getVector().distance(Vector3D.PLUS_J) < 1.0e-10;
            zVFound = zVFound || v.getLocation().getVector().distance(Vector3D.PLUS_K) < 1.0e-10;
        }
        assertTrue(xPFound);
        assertTrue(yPFound);
        assertTrue(zPFound);
        assertTrue(xVFound);
        assertTrue(yVFound);
        assertTrue(zVFound);
        assertEquals(3, count);

        assertEquals(0.0, octant.getBarycenter().distance(new S2Point(new Vector3D(1, 1, 1))), 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, octant.getSize(), 1.0e-10);

        EnclosingBall<Sphere2D, S2Point> cap = octant.getEnclosingCap();
        assertEquals(0.0, octant.getBarycenter().distance(cap.getCenter()), 1.0e-10);
        assertEquals(FastMath.acos(1.0 / FastMath.sqrt(3)), cap.getRadius(), 1.0e-10);

        EnclosingBall<Sphere2D, S2Point> reversedCap =
                ((SphericalPolygonsSet) factory.getComplement(octant)).getEnclosingCap();
        assertEquals(0, reversedCap.getCenter().distance(new S2Point(new Vector3D(-1, -1, -1))), 1.0e-10);
        assertEquals(FastMath.PI - FastMath.asin(1.0 / FastMath.sqrt(3)), reversedCap.getRadius(), 1.0e-10);

        assertEquals(Location.INSIDE, octant.checkPoint(octant.getInteriorPoint()));
    }

    @Test
    void testPositiveOctantByVertices() {
        double tol = 0.01;
        double sinTol = FastMath.sin(tol);
        SphericalPolygonsSet octant = new SphericalPolygonsSet(tol, S2Point.PLUS_I, S2Point.PLUS_J, S2Point.PLUS_K);
        UnitSphereRandomVectorGenerator random =
                new UnitSphereRandomVectorGenerator(3, new Well1024a(0xb8fc5acc91044308L));
        for (int i = 0; i < 1000; ++i) {
            Vector3D v = new Vector3D(random.nextVector());
            if ((v.getX() > sinTol) && (v.getY() > sinTol) && (v.getZ() > sinTol)) {
                assertEquals(Location.INSIDE, octant.checkPoint(new S2Point(v)));
            } else if ((v.getX() < -sinTol) || (v.getY() < -sinTol) || (v.getZ() < -sinTol)) {
                assertEquals(Location.OUTSIDE, octant.checkPoint(new S2Point(v)));
            } else {
                assertEquals(Location.BOUNDARY, octant.checkPoint(new S2Point(v)));
            }
        }

        assertEquals(Location.INSIDE, octant.checkPoint(octant.getInteriorPoint()));

    }

    @Test
    void testNonConvex() {
        double tol = 0.01;
        double sinTol = FastMath.sin(tol);
        RegionFactory<Sphere2D, S2Point, Circle, SubCircle> factory = new RegionFactory<>();
        SphericalPolygonsSet plusX = new SphericalPolygonsSet(Vector3D.PLUS_I, tol);
        SphericalPolygonsSet plusY = new SphericalPolygonsSet(Vector3D.PLUS_J, tol);
        SphericalPolygonsSet plusZ = new SphericalPolygonsSet(Vector3D.PLUS_K, tol);
        SphericalPolygonsSet threeOctants =
                (SphericalPolygonsSet) factory.difference(plusZ, factory.intersection(plusX, plusY));

        UnitSphereRandomVectorGenerator random =
                new UnitSphereRandomVectorGenerator(3, new Well1024a(0x9c9802fde3cbcf25L));
        for (int i = 0; i < 1000; ++i) {
            Vector3D v = new Vector3D(random.nextVector());
            if (((v.getX() < -sinTol) || (v.getY() < -sinTol)) && (v.getZ() > sinTol)) {
                assertEquals(Location.INSIDE, threeOctants.checkPoint(new S2Point(v)));
            } else if (((v.getX() > sinTol) && (v.getY() > sinTol)) || (v.getZ() < -sinTol)) {
                assertEquals(Location.OUTSIDE, threeOctants.checkPoint(new S2Point(v)));
            } else {
                assertEquals(Location.BOUNDARY, threeOctants.checkPoint(new S2Point(v)));
            }
        }

        List<Vertex> loops = threeOctants.getBoundaryLoops();
        assertEquals(1, loops.size());
        boolean xPFound = false;
        boolean yPFound = false;
        boolean zPFound = false;
        boolean xVFound = false;
        boolean yVFound = false;
        boolean zVFound = false;
        Vertex first = loops.get(0);
        int count = 0;
        double sumPoleX = 0;
        double sumPoleY = 0;
        double sumPoleZ = 0;
        for (Vertex v = first; count == 0 || v != first; v = v.getOutgoing().getEnd()) {
            ++count;
            Edge e = v.getIncoming();
            assertSame(v, e.getStart().getOutgoing().getEnd());
            if (e.getCircle().getPole().distance(Vector3D.MINUS_I) < 1.0e-10) {
                xPFound = true;
                sumPoleX += e.getLength();
            } else if (e.getCircle().getPole().distance(Vector3D.MINUS_J) < 1.0e-10) {
                yPFound = true;
                sumPoleY += e.getLength();
            } else {
                assertEquals(0.0, e.getCircle().getPole().distance(Vector3D.PLUS_K), 1.0e-10);
                zPFound = true;
                sumPoleZ += e.getLength();
            }
            xVFound = xVFound || v.getLocation().getVector().distance(Vector3D.PLUS_I) < 1.0e-10;
            yVFound = yVFound || v.getLocation().getVector().distance(Vector3D.PLUS_J) < 1.0e-10;
            zVFound = zVFound || v.getLocation().getVector().distance(Vector3D.PLUS_K) < 1.0e-10;
        }
        assertTrue(xPFound);
        assertTrue(yPFound);
        assertTrue(zPFound);
        assertTrue(xVFound);
        assertTrue(yVFound);
        assertTrue(zVFound);
        assertEquals(MathUtils.SEMI_PI, sumPoleX, 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, sumPoleY, 1.0e-10);
        assertEquals(1.5 * FastMath.PI, sumPoleZ, 1.0e-10);

        assertEquals(1.5 * FastMath.PI, threeOctants.getSize(), 1.0e-10);

        assertEquals(Location.INSIDE, threeOctants.checkPoint(threeOctants.getInteriorPoint()));
    }

    @Test
    void testModeratlyComplexShape() {
        double tol = 0.01;
        List<SubCircle> boundary = new ArrayList<>();
        boundary.add(create(Vector3D.MINUS_J, Vector3D.PLUS_I,  Vector3D.PLUS_K,  tol, 0.0, MathUtils.SEMI_PI));
        boundary.add(create(Vector3D.MINUS_I, Vector3D.PLUS_K,  Vector3D.PLUS_J,  tol, 0.0, MathUtils.SEMI_PI));
        boundary.add(create(Vector3D.PLUS_K,  Vector3D.PLUS_J,  Vector3D.MINUS_I, tol, 0.0, MathUtils.SEMI_PI));
        boundary.add(create(Vector3D.MINUS_J, Vector3D.MINUS_I, Vector3D.MINUS_K, tol, 0.0, MathUtils.SEMI_PI));
        boundary.add(create(Vector3D.MINUS_I, Vector3D.MINUS_K, Vector3D.MINUS_J, tol, 0.0, MathUtils.SEMI_PI));
        boundary.add(create(Vector3D.PLUS_K,  Vector3D.MINUS_J, Vector3D.PLUS_I,  tol, 0.0, MathUtils.SEMI_PI));
        SphericalPolygonsSet polygon = new SphericalPolygonsSet(boundary, tol);

        assertEquals(Location.OUTSIDE, polygon.checkPoint(new S2Point(new Vector3D( 1,  1,  1).normalize())));
        assertEquals(Location.INSIDE,  polygon.checkPoint(new S2Point(new Vector3D(-1,  1,  1).normalize())));
        assertEquals(Location.INSIDE,  polygon.checkPoint(new S2Point(new Vector3D(-1, -1,  1).normalize())));
        assertEquals(Location.INSIDE,  polygon.checkPoint(new S2Point(new Vector3D( 1, -1,  1).normalize())));
        assertEquals(Location.OUTSIDE, polygon.checkPoint(new S2Point(new Vector3D( 1,  1, -1).normalize())));
        assertEquals(Location.OUTSIDE, polygon.checkPoint(new S2Point(new Vector3D(-1,  1, -1).normalize())));
        assertEquals(Location.INSIDE,  polygon.checkPoint(new S2Point(new Vector3D(-1, -1, -1).normalize())));
        assertEquals(Location.OUTSIDE, polygon.checkPoint(new S2Point(new Vector3D( 1, -1, -1).normalize())));

        assertEquals(MathUtils.TWO_PI, polygon.getSize(), 1.0e-10);
        assertEquals(3 * FastMath.PI, polygon.getBoundarySize(), 1.0e-10);

        List<Vertex> loops = polygon.getBoundaryLoops();
        assertEquals(1, loops.size());
        boolean pXFound = false;
        boolean mXFound = false;
        boolean pYFound = false;
        boolean mYFound = false;
        boolean pZFound = false;
        boolean mZFound = false;
        Vertex first = loops.get(0);
        int count = 0;
        for (Vertex v = first; count == 0 || v != first; v = v.getOutgoing().getEnd()) {
            ++count;
            Edge e = v.getIncoming();
            assertSame(v, e.getStart().getOutgoing().getEnd());
            pXFound = pXFound || v.getLocation().getVector().distance(Vector3D.PLUS_I)  < 1.0e-10;
            mXFound = mXFound || v.getLocation().getVector().distance(Vector3D.MINUS_I) < 1.0e-10;
            pYFound = pYFound || v.getLocation().getVector().distance(Vector3D.PLUS_J)  < 1.0e-10;
            mYFound = mYFound || v.getLocation().getVector().distance(Vector3D.MINUS_J) < 1.0e-10;
            pZFound = pZFound || v.getLocation().getVector().distance(Vector3D.PLUS_K)  < 1.0e-10;
            mZFound = mZFound || v.getLocation().getVector().distance(Vector3D.MINUS_K) < 1.0e-10;
            assertEquals(MathUtils.SEMI_PI, e.getLength(), 1.0e-10);
        }
        assertTrue(pXFound);
        assertTrue(mXFound);
        assertTrue(pYFound);
        assertTrue(mYFound);
        assertTrue(pZFound);
        assertTrue(mZFound);
        assertEquals(6, count);

        assertEquals(Location.INSIDE, polygon.checkPoint(polygon.getInteriorPoint()));

    }

    @Test
    void testSeveralParts() {
        double tol = 0.01;
        double sinTol = FastMath.sin(tol);
        List<SubCircle> boundary = new ArrayList<>();

        // first part: +X, +Y, +Z octant
        boundary.add(create(Vector3D.PLUS_J,  Vector3D.PLUS_K,  Vector3D.PLUS_I,  tol, 0.0, MathUtils.SEMI_PI));
        boundary.add(create(Vector3D.PLUS_K,  Vector3D.PLUS_I,  Vector3D.PLUS_J,  tol, 0.0, MathUtils.SEMI_PI));
        boundary.add(create(Vector3D.PLUS_I,  Vector3D.PLUS_J,  Vector3D.PLUS_K,  tol, 0.0, MathUtils.SEMI_PI));

        // first part: -X, -Y, -Z octant
        boundary.add(create(Vector3D.MINUS_J, Vector3D.MINUS_I, Vector3D.MINUS_K, tol, 0.0, MathUtils.SEMI_PI));
        boundary.add(create(Vector3D.MINUS_I, Vector3D.MINUS_K, Vector3D.MINUS_J, tol, 0.0, MathUtils.SEMI_PI));
        boundary.add(create(Vector3D.MINUS_K, Vector3D.MINUS_J, Vector3D.MINUS_I,  tol, 0.0, MathUtils.SEMI_PI));

        SphericalPolygonsSet polygon = new SphericalPolygonsSet(boundary, tol);

        UnitSphereRandomVectorGenerator random =
                new UnitSphereRandomVectorGenerator(3, new Well1024a(0xcc5ce49949e0d3ecL));
        for (int i = 0; i < 1000; ++i) {
            Vector3D v = new Vector3D(random.nextVector());
            if ((v.getX() < -sinTol) && (v.getY() < -sinTol) && (v.getZ() < -sinTol)) {
                assertEquals(Location.INSIDE, polygon.checkPoint(new S2Point(v)));
            } else if ((v.getX() < sinTol) && (v.getY() < sinTol) && (v.getZ() < sinTol)) {
                assertEquals(Location.BOUNDARY, polygon.checkPoint(new S2Point(v)));
            } else if ((v.getX() > sinTol) && (v.getY() > sinTol) && (v.getZ() > sinTol)) {
                assertEquals(Location.INSIDE, polygon.checkPoint(new S2Point(v)));
            } else if ((v.getX() > -sinTol) && (v.getY() > -sinTol) && (v.getZ() > -sinTol)) {
                assertEquals(Location.BOUNDARY, polygon.checkPoint(new S2Point(v)));
            } else {
                assertEquals(Location.OUTSIDE, polygon.checkPoint(new S2Point(v)));
            }
        }

        assertEquals(FastMath.PI, polygon.getSize(), 1.0e-10);
        assertEquals(3 * FastMath.PI, polygon.getBoundarySize(), 1.0e-10);

        // there should be two separate boundary loops
        assertEquals(2, polygon.getBoundaryLoops().size());

        assertEquals(Location.INSIDE, polygon.checkPoint(polygon.getInteriorPoint()));

    }

    @Test
    void testPartWithHole() {
        double tol = 0.01;
        double alpha = 0.7;
        S2Point center = new S2Point(new Vector3D(1, 1, 1));
        SphericalPolygonsSet hexa = new SphericalPolygonsSet(center.getVector(), Vector3D.PLUS_K, alpha, 6, tol);
        SphericalPolygonsSet hole  = new SphericalPolygonsSet(tol,
                                                              new S2Point(FastMath.PI / 6, FastMath.PI / 3),
                                                              new S2Point(FastMath.PI / 3, FastMath.PI / 3),
                                                              new S2Point(FastMath.PI / 4, FastMath.PI / 6));
        SphericalPolygonsSet hexaWithHole =
                (SphericalPolygonsSet) new RegionFactory<Sphere2D, S2Point, Circle, SubCircle>().
                        difference(hexa, hole);

        for (double phi = center.getPhi() - alpha + 0.1; phi < center.getPhi() + alpha - 0.1; phi += 0.07) {
            Location l = hexaWithHole.checkPoint(new S2Point(FastMath.PI / 4, phi));
            if (phi < FastMath.PI / 6 || phi > FastMath.PI / 3) {
                assertEquals(Location.INSIDE,  l);
            } else {
                assertEquals(Location.OUTSIDE, l);
            }
        }

        // there should be two separate boundary loops
        assertEquals(2, hexaWithHole.getBoundaryLoops().size());

        assertEquals(hexa.getBoundarySize() + hole.getBoundarySize(), hexaWithHole.getBoundarySize(), 1.0e-10);
        assertEquals(hexa.getSize() - hole.getSize(), hexaWithHole.getSize(), 1.0e-10);

        assertEquals(Location.INSIDE, hexaWithHole.checkPoint(hexaWithHole.getInteriorPoint()));

    }

    @Test
    void testConcentricSubParts() {
        double tol = 0.001;
        Vector3D center = new Vector3D(1, 1, 1);
        SphericalPolygonsSet hexaOut   = new SphericalPolygonsSet(center, Vector3D.PLUS_K, 0.9,  6, tol);
        SphericalPolygonsSet hexaIn    = new SphericalPolygonsSet(center, Vector3D.PLUS_K, 0.8,  6, tol);
        SphericalPolygonsSet pentaOut  = new SphericalPolygonsSet(center, Vector3D.PLUS_K, 0.7,  5, tol);
        SphericalPolygonsSet pentaIn   = new SphericalPolygonsSet(center, Vector3D.PLUS_K, 0.6,  5, tol);
        SphericalPolygonsSet quadriOut = new SphericalPolygonsSet(center, Vector3D.PLUS_K, 0.5,  4, tol);
        SphericalPolygonsSet quadriIn  = new SphericalPolygonsSet(center, Vector3D.PLUS_K, 0.4,  4, tol);
        SphericalPolygonsSet triOut    = new SphericalPolygonsSet(center, Vector3D.PLUS_K, 0.25, 3, tol);
        SphericalPolygonsSet triIn     = new SphericalPolygonsSet(center, Vector3D.PLUS_K, 0.15, 3, tol);

        RegionFactory<Sphere2D, S2Point, Circle, SubCircle> factory = new RegionFactory<>();
        SphericalPolygonsSet hexa   = (SphericalPolygonsSet) factory.difference(hexaOut,   hexaIn);
        SphericalPolygonsSet penta  = (SphericalPolygonsSet) factory.difference(pentaOut,  pentaIn);
        SphericalPolygonsSet quadri = (SphericalPolygonsSet) factory.difference(quadriOut, quadriIn);
        SphericalPolygonsSet tri    = (SphericalPolygonsSet) factory.difference(triOut,    triIn);
        SphericalPolygonsSet concentric =
                (SphericalPolygonsSet) factory.union(factory.union(hexa, penta), factory.union(quadri, tri));

        // there should be two separate boundary loops
        assertEquals(8, concentric.getBoundaryLoops().size());

        assertEquals(hexaOut.getBoundarySize()   + hexaIn.getBoundarySize()   +
                            pentaOut.getBoundarySize()  + pentaIn.getBoundarySize()  +
                            quadriOut.getBoundarySize() + quadriIn.getBoundarySize() +
                            triOut.getBoundarySize()    + triIn.getBoundarySize(),
                            concentric.getBoundarySize(), 1.0e-10);
        assertEquals(hexaOut.getSize()   - hexaIn.getSize()   +
                            pentaOut.getSize()  - pentaIn.getSize()  +
                            quadriOut.getSize() - quadriIn.getSize() +
                            triOut.getSize()    - triIn.getSize(),
                            concentric.getSize(), 1.0e-10);

        // we expect lots of sign changes as we traverse all concentric rings
        double phi = new S2Point(center).getPhi();
        assertEquals(+0.207, concentric.projectToBoundary(new S2Point(-0.60,  phi)).getOffset(), 0.01);
        assertEquals(-0.048, concentric.projectToBoundary(new S2Point(-0.21,  phi)).getOffset(), 0.01);
        assertEquals(+0.027, concentric.projectToBoundary(new S2Point(-0.10,  phi)).getOffset(), 0.01);
        assertEquals(-0.041, concentric.projectToBoundary(new S2Point( 0.01,  phi)).getOffset(), 0.01);
        assertEquals(+0.049, concentric.projectToBoundary(new S2Point( 0.16,  phi)).getOffset(), 0.01);
        assertEquals(-0.038, concentric.projectToBoundary(new S2Point( 0.29,  phi)).getOffset(), 0.01);
        assertEquals(+0.097, concentric.projectToBoundary(new S2Point( 0.48,  phi)).getOffset(), 0.01);
        assertEquals(-0.022, concentric.projectToBoundary(new S2Point( 0.64,  phi)).getOffset(), 0.01);
        assertEquals(+0.072, concentric.projectToBoundary(new S2Point( 0.79,  phi)).getOffset(), 0.01);
        assertEquals(-0.022, concentric.projectToBoundary(new S2Point( 0.93,  phi)).getOffset(), 0.01);
        assertEquals(+0.091, concentric.projectToBoundary(new S2Point( 1.08,  phi)).getOffset(), 0.01);
        assertEquals(-0.037, concentric.projectToBoundary(new S2Point( 1.28,  phi)).getOffset(), 0.01);
        assertEquals(+0.051, concentric.projectToBoundary(new S2Point( 1.40,  phi)).getOffset(), 0.01);
        assertEquals(-0.041, concentric.projectToBoundary(new S2Point( 1.55,  phi)).getOffset(), 0.01);
        assertEquals(+0.027, concentric.projectToBoundary(new S2Point( 1.67,  phi)).getOffset(), 0.01);
        assertEquals(-0.044, concentric.projectToBoundary(new S2Point( 1.79,  phi)).getOffset(), 0.01);
        assertEquals(+0.201, concentric.projectToBoundary(new S2Point( 2.16,  phi)).getOffset(), 0.01);

        assertEquals(Location.INSIDE, concentric.checkPoint(concentric.getInteriorPoint()));

    }

    @Test
    void testGeographicalMap() {

        SphericalPolygonsSet continental = buildSimpleZone(new double[][] {
          { 51.14850,  2.51357 }, { 50.94660,  1.63900 }, { 50.12717,  1.33876 }, { 49.34737, -0.98946 },
          { 49.77634, -1.93349 }, { 48.64442, -1.61651 }, { 48.90169, -3.29581 }, { 48.68416, -4.59234 },
          { 47.95495, -4.49155 }, { 47.57032, -2.96327 }, { 46.01491, -1.19379 }, { 44.02261, -1.38422 },
          { 43.42280, -1.90135 }, { 43.03401, -1.50277 }, { 42.34338,  1.82679 }, { 42.47301,  2.98599 },
          { 43.07520,  3.10041 }, { 43.39965,  4.55696 }, { 43.12889,  6.52924 }, { 43.69384,  7.43518 },
          { 44.12790,  7.54959 }, { 45.02851,  6.74995 }, { 45.33309,  7.09665 }, { 46.42967,  6.50009 },
          { 46.27298,  6.02260 }, { 46.72577,  6.03738 }, { 47.62058,  7.46675 }, { 49.01778,  8.09927 },
          { 49.20195,  6.65822 }, { 49.44266,  5.89775 }, { 49.98537,  4.79922 }
        });
        SphericalPolygonsSet corsica = buildSimpleZone(new double[][] {
          { 42.15249,  9.56001 }, { 43.00998,  9.39000 }, { 42.62812,  8.74600 }, { 42.25651,  8.54421 },
          { 41.58361,  8.77572 }, { 41.38000,  9.22975 }
        });
        RegionFactory<Sphere2D, S2Point, Circle, SubCircle> factory = new RegionFactory<>();
        SphericalPolygonsSet zone = (SphericalPolygonsSet) factory.union(continental, corsica);
        EnclosingBall<Sphere2D, S2Point> enclosing = zone.getEnclosingCap();
        Vector3D enclosingCenter = enclosing.getCenter().getVector();

        double step = FastMath.toRadians(0.1);
        for (Vertex loopStart : zone.getBoundaryLoops()) {
            int count = 0;
            for (Vertex v = loopStart; count == 0 || v != loopStart; v = v.getOutgoing().getEnd()) {
                ++count;
                for (int i = 0; i < FastMath.ceil(v.getOutgoing().getLength() / step); ++i) {
                    Vector3D p = v.getOutgoing().getPointAt(i * step);
                    assertTrue(Vector3D.angle(p, enclosingCenter) <= enclosing.getRadius());
                }
            }
        }

        S2Point supportPointA = s2Point(48.68416, -4.59234);
        S2Point supportPointB = s2Point(41.38000,  9.22975);
        assertEquals(enclosing.getRadius(), supportPointA.distance(enclosing.getCenter()), 1.0e-10);
        assertEquals(enclosing.getRadius(), supportPointB.distance(enclosing.getCenter()), 1.0e-10);
        assertEquals(0.5 * supportPointA.distance(supportPointB), enclosing.getRadius(), 1.0e-10);
        assertEquals(2, enclosing.getSupportSize());

        EnclosingBall<Sphere2D, S2Point> continentalInscribed =
                ((SphericalPolygonsSet) factory.getComplement(continental)).getEnclosingCap();
        Vector3D continentalCenter = continentalInscribed.getCenter().getVector();
        assertEquals(2.2, FastMath.toDegrees(FastMath.PI - continentalInscribed.getRadius()), 0.1);
        for (Vertex loopStart : continental.getBoundaryLoops()) {
            int count = 0;
            for (Vertex v = loopStart; count == 0 || v != loopStart; v = v.getOutgoing().getEnd()) {
                ++count;
                for (int i = 0; i < FastMath.ceil(v.getOutgoing().getLength() / step); ++i) {
                    Vector3D p = v.getOutgoing().getPointAt(i * step);
                    assertTrue(Vector3D.angle(p, continentalCenter) <= continentalInscribed.getRadius());
                }
            }
        }

        EnclosingBall<Sphere2D, S2Point> corsicaInscribed =
                ((SphericalPolygonsSet) factory.getComplement(corsica)).getEnclosingCap();
        Vector3D corsicaCenter = corsicaInscribed.getCenter().getVector();
        assertEquals(0.34, FastMath.toDegrees(FastMath.PI - corsicaInscribed.getRadius()), 0.01);
        for (Vertex loopStart : corsica.getBoundaryLoops()) {
            int count = 0;
            for (Vertex v = loopStart; count == 0 || v != loopStart; v = v.getOutgoing().getEnd()) {
                ++count;
                for (int i = 0; i < FastMath.ceil(v.getOutgoing().getLength() / step); ++i) {
                    Vector3D p = v.getOutgoing().getPointAt(i * step);
                    assertTrue(Vector3D.angle(p, corsicaCenter) <= corsicaInscribed.getRadius());
                }
            }
        }

        assertEquals(Location.INSIDE, zone.checkPoint(zone.getInteriorPoint()));

        // TODO: this triggers an exception as we cannot find an interior point
//        checkInterior(FastMath.toRadians(-2), FastMath.toRadians(5), FastMath.toRadians(0.4),
//                      FastMath.toRadians(90 - 55), FastMath.toRadians(90 - 40), FastMath.toRadians(0.4),
//                      zone, 1.0e-8);

    }

    @Test
    void testZigZagBoundary() {
        SphericalPolygonsSet zone = new SphericalPolygonsSet(1.0e-6,
                                                             new S2Point(-0.12630940610562444, 0.8998192093789258),
                                                             new S2Point(-0.12731320182988207, 0.8963735568774486),
                                                             new S2Point(-0.1351107624622557,  0.8978258663483273),
                                                             new S2Point(-0.13545331405131725, 0.8966781238246179),
                                                             new S2Point(-0.14324883017454967, 0.8981309629283796),
                                                             new S2Point(-0.14359875625524995, 0.896983965573036),
                                                             new S2Point(-0.14749650541159384, 0.8977109994666864),
                                                             new S2Point(-0.14785037758231825, 0.8965644005442432),
                                                             new S2Point(-0.15369807257448784, 0.8976550608135502),
                                                             new S2Point(-0.1526225554339386,  0.9010934265410458),
                                                             new S2Point(-0.14679028466684121, 0.9000043396997698),
                                                             new S2Point(-0.14643807494172612, 0.9011511073761742),
                                                             new S2Point(-0.1386609051963748,  0.8996991539048602),
                                                             new S2Point(-0.13831601655974668, 0.9008466623902937),
                                                             new S2Point(-0.1305365419828323,  0.8993961857946309),
                                                             new S2Point(-0.1301989630405964,  0.9005444294061787));
        assertEquals(Region.Location.INSIDE, zone.checkPoint(new S2Point(-0.145, 0.898)));
        assertEquals(6.463e-5, zone.getSize(),         1.0e-7);
        assertEquals(5.487e-2, zone.getBoundarySize(), 1.0e-4);
    }

    @Test
    void testGitHubIssue41() {
        RegionFactory<Sphere2D, S2Point, Circle, SubCircle> regionFactory = new RegionFactory<>();
        S2Point[] s2pA = new S2Point[]{
                new S2Point(new Vector3D(0.2122954606, -0.629606302,  0.7473463333)),
                new S2Point(new Vector3D(0.2120220248, -0.6296445493, 0.747391733)),
                new S2Point(new Vector3D(0.2119838016, -0.6298173178, 0.7472569934)),
                new S2Point(new Vector3D(0.2122571927, -0.6297790738, 0.7472116182))};

        S2Point[] s2pB = new S2Point[]{
                new S2Point(new Vector3D(0.2120291561, -0.629952069,  0.7471305292)),
                new S2Point(new Vector3D(0.2123026002, -0.6299138005, 0.7470851423)),
                new S2Point(new Vector3D(0.2123408927, -0.6297410403, 0.7472198923)),
                new S2Point(new Vector3D(0.2120674039, -0.6297793122, 0.7472653037))};

        final double tol = 0.0001;
        final SphericalPolygonsSet spsA = new SphericalPolygonsSet(tol, s2pA);
        final SphericalPolygonsSet spsB = new SphericalPolygonsSet(tol, s2pB);
        assertEquals(0.61254e-7, spsA.getSize(),          1.0e-12);
        assertEquals(1.00437e-3, spsA.getBoundarySize(),  1.0e-08);
        assertEquals(0.61269e-7, spsB.getSize(),          1.0e-12);
        assertEquals(1.00452e-3, spsB.getBoundarySize(),  1.0e-08);
        SphericalPolygonsSet union = (SphericalPolygonsSet) regionFactory.union(spsA, spsB);

        // as the tolerance is very large with respect to polygons,
        // the union is not computed properly, which is EXPECTED
        // so the thresholds for the tests are large.
        // the reference values have been computed with a much lower tolerance
        assertEquals(1.15628e-7, union.getSize(),         4.0e-9);
        assertEquals(1.53824e-3, union.getBoundarySize(), 3.0e-4);

    }

    @Test
    void testGitHubIssue42A() {
        // if building it was allowed (i.e. if the check for tolerance was removed)
        // the BSP tree would be wrong, it would include a large extra chunk that contains
        // a point that should really be outside
        try {
            doTestGitHubIssue42(1.0e-100);
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedGeometryFormats.TOO_SMALL_TOLERANCE, miae.getSpecifier());
            assertEquals(1.0e-100, (Double) miae.getParts()[0], 1.0e-110);
            assertEquals("Sphere2D.SMALLEST_TOLERANCE", miae.getParts()[1]);
            assertEquals(Sphere2D.SMALLEST_TOLERANCE, (Double) miae.getParts()[2], 1.0e-20);
        }
    }

    @Test
    void testGitHubIssue42B() {
        // the BSP tree is right, but size cannot be computed
        try {
            // success of this call is dependent on numerical noise.
            // If it fails it should fail predictably.
            doTestGitHubIssue42(9.0e-16);
        } catch (MathIllegalStateException e) {
            assertEquals(
                LocalizedGeometryFormats.OUTLINE_BOUNDARY_LOOP_OPEN, e.getSpecifier());
        }
        // Computations in Edge.split use angles up to 4 PI
        double tol = FastMath.ulp(4 * FastMath.PI);
        // works when tol >= ulp(largest angle used)
        doTestGitHubIssue42(tol);
    }

    private void doTestGitHubIssue42(double tolerance) throws MathIllegalArgumentException {
        S2Point[] s2pA = new S2Point[]{
            new S2Point(new Vector3D(0.1504230736114679,  -0.6603084987333554, 0.7357754993377947)),
            new S2Point(new Vector3D(0.15011191112224423, -0.6603400871954631, 0.7358106980616113)),
            new S2Point(new Vector3D(0.15008035620222715, -0.6605195692153062, 0.7356560238085725)),
            new S2Point(new Vector3D(0.1503914563063968,  -0.6604879854490165, 0.7356208472763267))
        };
        S2Point outsidePoint = new S2Point(new Vector3D( 2, s2pA[0].getVector(),
                                                        -1, s2pA[1].getVector(),
                                                        -1, s2pA[2].getVector(),
                                                         2, s2pA[3].getVector()).normalize());
        // test all permutations because order matters when tolerance is small
        for (int i = 0; i < s2pA.length; i++) {
            S2Point[] points = new S2Point[s2pA.length];
            for (int j = 0; j < s2pA.length; j++) {
                points[j] = s2pA[(i + j) % s2pA.length];
            }
            final SphericalPolygonsSet spsA = new SphericalPolygonsSet(tolerance, points);
            assertEquals(Location.OUTSIDE, spsA.checkPoint(outsidePoint));
            assertEquals(7.4547e-8, spsA.getSize(), 1.0e-12);
        }
    }

    @Test
    void testZigZagBoundaryOversampledIssue46() {
        final double tol = 1.0e-4;
        // sample region, non-convex, not too big, not too small
        final S2Point[] vertices = {
                new S2Point(-0.12630940610562444e1, (0.8998192093789258 - 0.89) * 100),
                new S2Point(-0.12731320182988207e1, (0.8963735568774486 - 0.89) * 100),
                new S2Point(-0.1351107624622557e1, (0.8978258663483273 - 0.89) * 100),
                new S2Point(-0.13545331405131725e1, (0.8966781238246179 - 0.89) * 100),
                new S2Point(-0.14324883017454967e1, (0.8981309629283796 - 0.89) * 100),
                new S2Point(-0.14359875625524995e1, (0.896983965573036 - 0.89) * 100),
                new S2Point(-0.14749650541159384e1, (0.8977109994666864 - 0.89) * 100),
                new S2Point(-0.14785037758231825e1, (0.8965644005442432 - 0.89) * 100),
                new S2Point(-0.15369807257448784e1, (0.8976550608135502 - 0.89) * 100),
                new S2Point(-0.1526225554339386e1, (0.9010934265410458 - 0.89) * 100),
                new S2Point(-0.14679028466684121e1, (0.9000043396997698 - 0.89) * 100),
                new S2Point(-0.14643807494172612e1, (0.9011511073761742 - 0.89) * 100),
                new S2Point(-0.1386609051963748e1, (0.8996991539048602 - 0.89) * 100),
                new S2Point(-0.13831601655974668e1, (0.9008466623902937 - 0.89) * 100),
                new S2Point(-0.1305365419828323e1, (0.8993961857946309 - 0.89) * 100),
                new S2Point(-0.1301989630405964e1, (0.9005444294061787 - 0.89) * 100)};
        SphericalPolygonsSet zone = new SphericalPolygonsSet(tol, vertices);
        // sample high resolution boundary
        List<S2Point> points = new ArrayList<>();
        final Vertex start = zone.getBoundaryLoops().get(0);
        Vertex v = start;
        double step = tol / 10;
        do {
            Edge outgoing = v.getOutgoing();
            final double length = outgoing.getLength();
            int n = (int) (length / step);
            for (int i = 0; i < n; i++) {
                points.add(new S2Point(outgoing.getPointAt(i * step)));
            }
            v = outgoing.getEnd();
        } while (v != start);
        // create zone from high resolution boundary
        zone = new SphericalPolygonsSet(tol, points.toArray(new S2Point[0]));
        EnclosingBall<Sphere2D, S2Point> cap = zone.getEnclosingCap();
        // check cap size is reasonable. The region is ~0.5 accross, could be < 0.25
        assertTrue(cap.getRadius() < 0.5);
        assertEquals(Location.INSIDE, zone.checkPoint(zone.getBarycenter()));
        assertEquals(Location.INSIDE, zone.checkPoint(cap.getCenter()));
        // extra tolerance at corners due to SPS tolerance being a hyperplaneThickness
        // a factor of 3.1 corresponds to a edge intersection angle of ~19 degrees
        final double cornerTol = 3.1 * tol;
        for (S2Point vertex : vertices) {
            // check original points are on the boundary
            assertEquals(Location.BOUNDARY, zone.checkPoint(vertex), "" + vertex);
            double offset = FastMath.abs(zone.projectToBoundary(vertex).getOffset());
            assertEquals(0, offset, cornerTol, vertex + " offset: " + offset);
            // check original points are within the cap
            assertTrue(
                    cap.contains(vertex, tol),
                    "vertex: " + vertex + " distance: " + (vertex.distance(cap.getCenter()) - cap.getRadius()));
        }
    }

    @Test
    void testPositiveOctantByVerticesDetailIssue46() {
        double tol = 0.01;
        double sinTol = FastMath.sin(tol);
        Circle x = new Circle(Vector3D.PLUS_I, tol);
        Circle z = new Circle(Vector3D.PLUS_K, tol);
        double length = FastMath.PI / 2;
        double step = tol / 10;
        // sample high resolution boundary
        int n = (int) (length / step);
        List<S2Point> points = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            double t = i * step;
            points.add(new S2Point(z.getPointAt(z.getPhase(Vector3D.PLUS_I) + t)));
        }
        for (int i = 0; i < n; i++) {
            double t = i * step;
            points.add(new S2Point(x.getPointAt(x.getPhase(Vector3D.PLUS_J) + t)));
        }
        points.add(S2Point.PLUS_K);

        SphericalPolygonsSet octant = new SphericalPolygonsSet(tol, points.toArray(new S2Point[0]));
        UnitSphereRandomVectorGenerator random =
                new UnitSphereRandomVectorGenerator(3, new Well1024a(0xb8fc5acc91044308L));
        /* Where exactly the boundaries fall depends on which points are kept from
         * decimation, which can vary by up to tol. So a point up to 2*tol away from a
         * input point may be on the boundary. All input points are guaranteed to be on
         * the boundary, just not the center of the boundary.
         */
        for (int i = 0; i < 1000; ++i) {
            Vector3D v = new Vector3D(random.nextVector());
            final Location actual = octant.checkPoint(new S2Point(v));
            if ((v.getX() > sinTol) && (v.getY() > sinTol) && (v.getZ() > sinTol)) {
                if ((v.getX() > 2*sinTol) && (v.getY() > 2*sinTol) && (v.getZ() > 2*sinTol)) {
                    // certainly inside
                    assertEquals(Location.INSIDE, actual, "" + v);
                } else {
                    // may be inside or boundary
                    assertNotEquals(Location.OUTSIDE, actual, "" + v);
                }
            } else if ((v.getX() < 0) || (v.getY() < 0) || (v.getZ() < 0)) {
                if ((v.getX() < -sinTol) || (v.getY() < -sinTol) || (v.getZ() < -sinTol)) {
                    // certainly outside
                    assertEquals(Location.OUTSIDE, actual);
                } else {
                    // may be outside or boundary
                    assertNotEquals(Location.INSIDE, actual);
                }
            } else {
                // certainly boundary
                assertEquals(Location.BOUNDARY, actual);
            }
        }
        // all input points are on the boundary
        for (S2Point point : points) {
            assertEquals(Location.BOUNDARY, octant.checkPoint(point), "" + point);
        }
    }

    /** Check that constructing a region from 0, 1, 2, 3 points along a circle works. */
    @Test
    void testConstructingFromFewPointsIssue46() {
        double tol  = 1e-9;
        List<S2Point> points = new ArrayList<>();
        Circle circle = new Circle(Vector3D.PLUS_K, tol);

        // no points
        SphericalPolygonsSet sps = new SphericalPolygonsSet(tol, new S2Point[0]);
        assertEquals(4*FastMath.PI, sps.getSize(), tol);

        // one point, does not define a valid boundary
        points.add(new S2Point(circle.getPointAt(0)));
        try {
            new SphericalPolygonsSet(tol, points.toArray(new S2Point[0]));
            fail("expcected exception");
        } catch (MathRuntimeException e) {
            // expected
        }

        // two points, defines hemisphere but not orientation
        points.add(new S2Point(circle.getPointAt(FastMath.PI / 2)));
        sps = new SphericalPolygonsSet(tol, points.toArray(new S2Point[0]));
        assertEquals(2*FastMath.PI, sps.getSize(), tol);

        // three points
        points.add(0, new S2Point(circle.getPointAt(FastMath.PI)));
        sps = new SphericalPolygonsSet(tol, points.toArray(new S2Point[0]));
        assertEquals(2*FastMath.PI, sps.getSize(), tol);
        assertEquals(0, sps.getBarycenter().distance(new S2Point(Vector3D.PLUS_K)), tol);

        // four points
        points.add(1, new S2Point(circle.getPointAt(3 * FastMath.PI / 2)));
        sps = new SphericalPolygonsSet(tol, points.toArray(new S2Point[0]));
        assertEquals(2 * FastMath.PI, sps.getSize(), tol);
        assertEquals(0, sps.getBarycenter().distance(new S2Point(Vector3D.PLUS_K)), tol);

        // many points in semi-circle
        sps = new SphericalPolygonsSet(tol,
                new S2Point(circle.getPointAt(0)),
                new S2Point(circle.getPointAt(-0.3)),
                new S2Point(circle.getPointAt(-0.2)),
                new S2Point(circle.getPointAt(-0.1)));
        assertEquals(2 * FastMath.PI, sps.getSize(), tol);
        assertEquals(0, sps.getBarycenter().distance(new S2Point(Vector3D.PLUS_K)), tol);
    }

    @Test
    void testDefensiveProgrammingCheck() {
        // this tests defensive programming code that seems almost unreachable otherwise
        try {
            Method searchHelper = SphericalPolygonsSet.class.getDeclaredMethod("searchHelper",
                                                                               IntPredicate.class,
                                                                               Integer.TYPE, Integer.TYPE);
            searchHelper.setAccessible(true);
            searchHelper.invoke(null, (IntPredicate) (n -> true), 1, 0);
            fail("an exception should have been thrown");
        } catch (InvocationTargetException ite) {
            MathIllegalArgumentException miae = (MathIllegalArgumentException) ite.getCause();
            assertEquals(LocalizedCoreFormats.LOWER_ENDPOINT_ABOVE_UPPER_ENDPOINT, miae.getSpecifier());
        } catch (NoSuchMethodException | IllegalAccessException | IllegalArgumentException e) {
            fail(e.getLocalizedMessage());
        }
    }

    /**
     * Tests the Hipparchus {@link RegionFactory#intersection(Region, Region)}
     * method.
     */
    @Test
    void TestIntersectionOrder() {

        final S2Point[] vertices1 = {
            new S2Point(0.0193428339344826, 1.5537209444301618),
            new S2Point(0.0178197572212936, 1.553415699912148),
            new S2Point(0.01628496406053076, 1.5531081515279537),
            new S2Point(0.016284670226196844, 1.5531096373947835),
            new S2Point(0.019342540199680208, 1.5537224293848613)
        };

        final S2Point[] vertices2 = {
            new S2Point(0.016, 1.555),
            new S2Point(0.017453292519943295, 1.555),
            new S2Point(0.017453292519943295, 1.5533430342749532),
            new S2Point(0.016, 1.5533430342749532)
        };

        final RegionFactory<Sphere2D, S2Point, Circle, SubCircle> regionFactory = new RegionFactory<>();

        // thickness is small enough for proper computation of very small intersection
        double thickness1 = 4.96740426e-11;
        final SphericalPolygonsSet sps1 = new SphericalPolygonsSet(thickness1, vertices1);
        final SphericalPolygonsSet sps2 = new SphericalPolygonsSet(thickness1, vertices2);
        assertEquals(1.4886e-12, regionFactory.intersection(sps1, sps2).getSize(), 1.0e-15);
        assertEquals(1.4881e-12, regionFactory.intersection(sps2, sps1).getSize(), 1.0e-15);

        // thickness is too large, very small intersection is not computed properly in one case
        double thickness2 = 4.96740427e-11;
        final SphericalPolygonsSet sps3 = new SphericalPolygonsSet(thickness2, vertices1);
        final SphericalPolygonsSet sps4 = new SphericalPolygonsSet(thickness2, vertices2);
        assertEquals(1.4886e-12, regionFactory.intersection(sps3, sps4).getSize(), 1.0e-15);
        assertEquals(0.0,        regionFactory.intersection(sps4, sps3).getSize(), 1.0e-15);

    }

    @Test
    public void testIssueOrekit1388A() {
        doTestIssueOrekit1388(true);
    }

    @Test
    public void testIssueOrekit1388B() {
        doTestIssueOrekit1388(false);
    }

    private void doTestIssueOrekit1388(final boolean order) {
        final double[][] coordinates1 = new double[][] {
            { 18.52684751402596,  -76.97880893719434 },
            { 18.451108862175584, -76.99484778988442 },
            { 18.375369256045143, -77.01087679277504 },
            { 18.299628701801268, -77.02689599675749 },
            { 18.223887203723567, -77.04290545732495 },
            { 18.148144771769385, -77.05890521550272 },
            { 18.072401410187016, -77.0748953260324  },
            { 17.99665712514784,  -77.09087584010511 },
            { 17.92091192468999,  -77.10684680260262 },
            { 17.84516581113023,  -77.12280827309398 },
            { 17.769418792522664, -77.13876029654094 },
            { 17.747659863099422, -77.14334087084347 },
            { 17.67571798336192,  -77.15846791369165 },
            { 17.624293265977183, -77.16928381433733 },
            { 17.5485398681768,   -77.18520934447962 },
            { 17.526779103104783, -77.1897823275402  },
            { 17.49650619905315,  -77.0342031192472  },
            { 17.588661518962343, -77.01473903648854 },
            { 17.728574326965138, -76.98517769352242 },
            { 17.80416324015021,  -76.96919708557023 },
            { 17.87969526622326,  -76.95321858415689 },
            { 17.955280973332677, -76.93721874766547 },
            { 18.030855567607098, -76.92121123297645 },
            { 18.106414929680927, -76.90519686376611 },
            { 18.182031502555215, -76.88916022728444 },
            { 18.257597934434987, -76.87312403715188 },
            { 18.3331742522667,   -76.85707550881591 },
            { 18.408750874895002, -76.84101662269072 },
            { 18.57249082100609,  -76.80620195239239 },
            { 18.602585205425896, -76.96276018365842 }
        };

        final double[][] coordinates2 = new double[][] {
            { 18.338614038907608, -78.37885677406668 },
            { 18.195574802144037, -78.24425107003432 },
            { 18.20775293886321,  -78.0711865934217  },
            { 18.07679345301507,  -77.95901517339438 },
            { 18.006705181057598, -77.85325354879791 },
            { 17.857293838883137, -77.73787723105598 },
            { 17.854243316622103, -77.57442744758828 },
            { 17.875595873376014, -77.38213358468467 },
            { 17.72607423578937,  -77.23470828979222 },
            { 17.71386286451302,  -77.12253686976543 },
            { 17.790170276013725, -77.14817605148616 },
            { 17.869495404611797, -77.14497115377101 },
            { 17.854243309397717, -76.9302429967729  },
            { 17.954882874700132, -76.84371075846688 },
            { 17.94268718313505,  -76.6898756681441  },
            { 17.869495397388064, -76.54886016868198 },
            { 17.863394719203555, -76.35015651034861 },
            { 17.93049065091843,  -76.23478019260665 },
            { 18.155989976776553, -76.32451732862788 },
            { 18.22601854027039,  -76.63218750927341 },
            { 18.33861403170316,  -76.85653034932697 },
            { 18.405527980074993, -76.97831646249921 },
            { 18.4541763474828,   -77.28598664314421 },
            { 18.496732365966466, -77.705828243816   },
            { 18.451136227912485, -78.00708862903122 },
            { 18.405527980074993, -78.25707065080552 }
        };

        final double[][] expectedIn = new double[][] {
                { 18.408, -77.003 },
                { 18.338, -76.857 },
                { 17.869, -77.117 },
                { 17.857, -76.959 },
                { 17.761, -77.139 },
                { 17.715, -77.125 },
                { 17.935, -77.055 }
        };

        final double[][] expectedOut = new double[][] {
                { 17.794, -77.145 },
                { 17.736, -76.981 },
                { 17.715, -77.138 },
                { 18.153, -77.059 },
                { 18.232, -76.877 },
                { 18.373, -76.917 },
                { 17.871, -77.261 } // this is the point that was wrongly inside the intersection
        };

        SphericalPolygonsSet shape1 = buildSimpleZone(coordinates1);
        SphericalPolygonsSet shape2 = buildSimpleZone(coordinates2);
        Region<Sphere2D, S2Point, Circle, SubCircle> intersection =
                order ?
                new RegionFactory<Sphere2D, S2Point, Circle, SubCircle>().intersection(shape1.copySelf(), shape2.copySelf()) :
                new RegionFactory<Sphere2D, S2Point, Circle, SubCircle>().intersection(shape2.copySelf(), shape1.copySelf());

        for (final double[] doubles : expectedIn) {
            Assertions.assertEquals(Location.INSIDE, intersection.checkPoint(s2Point(doubles[0], doubles[1])));
        }

        for (final double[] doubles : expectedOut) {
            Assertions.assertEquals(Location.OUTSIDE, intersection.checkPoint(s2Point(doubles[0], doubles[1])));
        }

    }

    private SubCircle create(Vector3D pole, Vector3D x, Vector3D y,
                             double tolerance, double ... limits) {
        RegionFactory<Sphere1D, S1Point, LimitAngle, SubLimitAngle> factory = new RegionFactory<>();
        Circle                           circle  = new Circle(pole, tolerance);
        Circle phased = Circle.getTransform(new Rotation(circle.getXAxis(), circle.getYAxis(), x, y)).apply(circle);
        ArcsSet set = (ArcsSet) factory.getComplement(new ArcsSet(tolerance));
        for (int i = 0; i < limits.length; i += 2) {
            set = (ArcsSet) factory.union(set, new ArcsSet(limits[i], limits[i + 1], tolerance));
        }
        return new SubCircle(phased, set);
    }

    private SphericalPolygonsSet buildSimpleZone(double[][] points) {
        final S2Point[] vertices = new S2Point[points.length];
        for (int i = 0; i < points.length; ++i) {
            vertices[i] = s2Point(points[i][0], points[i][1]);
        }
        return new SphericalPolygonsSet(1.0e-10, vertices);
    }

    private S2Point s2Point(double latitude, double longitude) {
        return new S2Point(FastMath.toRadians(longitude), FastMath.toRadians(90.0 - latitude));
    }

}
