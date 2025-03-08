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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.geometry.partitioning.BSPTree;
import org.hipparchus.geometry.partitioning.Region;
import org.hipparchus.geometry.partitioning.Region.Location;
import org.hipparchus.geometry.partitioning.RegionFactory;
import org.hipparchus.geometry.partitioning.Side;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class ArcsSetTest {

    @Test
    void testArc() {
        ArcsSet set = new ArcsSet(2.3, 5.7, 1.0e-10);
        assertEquals(3.4, set.getSize(), 1.0e-10);
        assertEquals(1.0e-10, set.getTolerance(), 1.0e-20);
        assertEquals(Region.Location.BOUNDARY, set.checkPoint(new S1Point(2.3)));
        assertEquals(Region.Location.BOUNDARY, set.checkPoint(new S1Point(5.7)));
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new S1Point(1.2)));
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new S1Point(8.5)));
        assertEquals(Region.Location.INSIDE,   set.checkPoint(new S1Point(8.7)));
        assertEquals(Region.Location.INSIDE,   set.checkPoint(new S1Point(3.0)));
        assertEquals(1, set.asList().size());
        assertEquals(2.3, set.asList().get(0).getInf(), 1.0e-10);
        assertEquals(5.7, set.asList().get(0).getSup(), 1.0e-10);
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));
    }

    @Test
    void testWrapAround2PiArc() {
        ArcsSet set = new ArcsSet(5.7 - MathUtils.TWO_PI, 2.3, 1.0e-10);
        assertEquals(MathUtils.TWO_PI - 3.4, set.getSize(), 1.0e-10);
        assertEquals(1.0e-10, set.getTolerance(), 1.0e-20);
        assertEquals(Region.Location.BOUNDARY, set.checkPoint(new S1Point(2.3)));
        assertEquals(Region.Location.BOUNDARY, set.checkPoint(new S1Point(5.7)));
        assertEquals(Region.Location.INSIDE,   set.checkPoint(new S1Point(1.2)));
        assertEquals(Region.Location.INSIDE,   set.checkPoint(new S1Point(8.5)));
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new S1Point(8.7)));
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new S1Point(3.0)));
        assertEquals(1, set.asList().size());
        assertEquals(5.7, set.asList().get(0).getInf(), 1.0e-10);
        assertEquals(2.3 + MathUtils.TWO_PI, set.asList().get(0).getSup(), 1.0e-10);
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));
    }

    @Test
    void testSplitOver2Pi() {
        ArcsSet set = new ArcsSet(1.0e-10);
        Arc     arc = new Arc(1.5 * FastMath.PI, 2.5 * FastMath.PI, 1.0e-10);
        ArcsSet.Split split = set.split(arc);
        for (double alpha = 0.0; alpha <= MathUtils.TWO_PI; alpha += 0.01) {
            S1Point p = new S1Point(alpha);
            if (alpha < MathUtils.SEMI_PI || alpha > 1.5 * FastMath.PI) {
                assertEquals(Location.OUTSIDE, split.getPlus().checkPoint(p));
                assertEquals(Location.INSIDE,  split.getMinus().checkPoint(p));
            } else {
                assertEquals(Location.INSIDE,  split.getPlus().checkPoint(p));
                assertEquals(Location.OUTSIDE, split.getMinus().checkPoint(p));
            }
        }
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));
    }

    @Test
    void testSplitAtEnd() {
        ArcsSet set = new ArcsSet(1.0e-10);
        Arc     arc = new Arc(FastMath.PI, MathUtils.TWO_PI, 1.0e-10);
        ArcsSet.Split split = set.split(arc);
        for (double alpha = 0.01; alpha < MathUtils.TWO_PI; alpha += 0.01) {
            S1Point p = new S1Point(alpha);
            if (alpha > FastMath.PI) {
                assertEquals(Location.OUTSIDE, split.getPlus().checkPoint(p));
                assertEquals(Location.INSIDE,  split.getMinus().checkPoint(p));
            } else {
                assertEquals(Location.INSIDE,  split.getPlus().checkPoint(p));
                assertEquals(Location.OUTSIDE, split.getMinus().checkPoint(p));
            }
        }

        S1Point zero = new S1Point(0.0);
        assertEquals(Location.BOUNDARY,  split.getPlus().checkPoint(zero));
        assertEquals(Location.BOUNDARY,  split.getMinus().checkPoint(zero));

        S1Point pi = new S1Point(FastMath.PI);
        assertEquals(Location.BOUNDARY,  split.getPlus().checkPoint(pi));
        assertEquals(Location.BOUNDARY,  split.getMinus().checkPoint(pi));

    }

    @Test
    void testWrongInterval() {
        assertThrows(MathIllegalArgumentException.class, () -> new ArcsSet(1.2, 0.0, 1.0e-10));
    }

    @Test
    void testTooSmallTolerance() {
        assertThrows(MathIllegalArgumentException.class, () -> new ArcsSet(0.0, 1.0, 0.9 * Sphere1D.SMALLEST_TOLERANCE));
    }

    @Test
    void testFullEqualEndPoints() {
        ArcsSet set = new ArcsSet(1.0, 1.0, 1.0e-10);
        assertEquals(1.0e-10, set.getTolerance(), 1.0e-20);
        assertEquals(Region.Location.INSIDE, set.checkPoint(new S1Point(9.0)));
        for (double alpha = -20.0; alpha <= 20.0; alpha += 0.1) {
            assertEquals(Region.Location.INSIDE, set.checkPoint(new S1Point(alpha)));
        }
        assertEquals(1, set.asList().size());
        assertEquals(0.0, set.asList().get(0).getInf(), 1.0e-10);
        assertEquals(2 * FastMath.PI, set.asList().get(0).getSup(), 1.0e-10);
        assertEquals(2 * FastMath.PI, set.getSize(), 1.0e-10);
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));
    }

    @Test
    void testFullCircle() {
        ArcsSet set = new ArcsSet(1.0e-10);
        assertEquals(1.0e-10, set.getTolerance(), 1.0e-20);
        assertEquals(Region.Location.INSIDE, set.checkPoint(new S1Point(9.0)));
        for (double alpha = -20.0; alpha <= 20.0; alpha += 0.1) {
            assertEquals(Region.Location.INSIDE, set.checkPoint(new S1Point(alpha)));
        }
        assertEquals(1, set.asList().size());
        assertEquals(0.0, set.asList().get(0).getInf(), 1.0e-10);
        assertEquals(2 * FastMath.PI, set.asList().get(0).getSup(), 1.0e-10);
        assertEquals(2 * FastMath.PI, set.getSize(), 1.0e-10);
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));
    }

    @Test
    void testEmpty() {
        ArcsSet empty = (ArcsSet) new RegionFactory<Sphere1D, S1Point, LimitAngle, SubLimitAngle>().getComplement(new ArcsSet(1.0e-10));
        assertEquals(1.0e-10, empty.getTolerance(), 1.0e-20);
        assertEquals(0.0, empty.getSize(), 1.0e-10);
        assertTrue(empty.asList().isEmpty());
        assertNull(empty.getInteriorPoint());
    }

    @Test
    void testTiny() {
        ArcsSet tiny = new ArcsSet(0.0, Precision.SAFE_MIN / 2, 1.0e-10);
        assertEquals(1.0e-10, tiny.getTolerance(), 1.0e-20);
        assertEquals(Precision.SAFE_MIN / 2, tiny.getSize(), 1.0e-10);
        assertEquals(1, tiny.asList().size());
        assertEquals(0.0, tiny.asList().get(0).getInf(), 1.0e-10);
        assertEquals(Precision.SAFE_MIN / 2, tiny.asList().get(0).getSup(), 1.0e-10);
        assertEquals(Location.BOUNDARY, tiny.checkPoint(tiny.getInteriorPoint()));
    }

    @Test
    void testSpecialConstruction() {
        List<SubLimitAngle> boundary = new ArrayList<>();
        boundary.add(new LimitAngle(new S1Point(0.0), false, 1.0e-10).wholeHyperplane());
        boundary.add(new LimitAngle(new S1Point(MathUtils.TWO_PI - 1.0e-11), true, 1.0e-10).wholeHyperplane());
        ArcsSet set = new ArcsSet(boundary, 1.0e-10);
        assertEquals(MathUtils.TWO_PI, set.getSize(), 1.0e-10);
        assertEquals(1.0e-10, set.getTolerance(), 1.0e-20);
        assertEquals(1, set.asList().size());
        assertEquals(0.0, set.asList().get(0).getInf(), 1.0e-10);
        assertEquals(MathUtils.TWO_PI, set.asList().get(0).getSup(), 1.0e-10);
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));
    }

    @Test
    void testDifference() {

        ArcsSet a   = new ArcsSet(1.0, 6.0, 1.0e-10);
        List<Arc> aList = a.asList();
        assertEquals(1,   aList.size());
        assertEquals(1.0, aList.get(0).getInf(), 1.0e-10);
        assertEquals(6.0, aList.get(0).getSup(), 1.0e-10);

        ArcsSet b   = new ArcsSet(3.0, 5.0, 1.0e-10);
        List<Arc> bList = b.asList();
        assertEquals(1,   bList.size());
        assertEquals(3.0, bList.get(0).getInf(), 1.0e-10);
        assertEquals(5.0, bList.get(0).getSup(), 1.0e-10);

        ArcsSet aMb = (ArcsSet) new RegionFactory<Sphere1D, S1Point, LimitAngle, SubLimitAngle>().difference(a, b);
        for (int k = -2; k < 3; ++k) {
            assertEquals(Location.OUTSIDE,  aMb.checkPoint(new S1Point(0.0 + k * MathUtils.TWO_PI)));
            assertEquals(Location.OUTSIDE,  aMb.checkPoint(new S1Point(0.9 + k * MathUtils.TWO_PI)));
            assertEquals(Location.BOUNDARY, aMb.checkPoint(new S1Point(1.0 + k * MathUtils.TWO_PI)));
            assertEquals(Location.INSIDE,   aMb.checkPoint(new S1Point(1.1 + k * MathUtils.TWO_PI)));
            assertEquals(Location.INSIDE,   aMb.checkPoint(new S1Point(2.9 + k * MathUtils.TWO_PI)));
            assertEquals(Location.BOUNDARY, aMb.checkPoint(new S1Point(3.0 + k * MathUtils.TWO_PI)));
            assertEquals(Location.OUTSIDE,  aMb.checkPoint(new S1Point(3.1 + k * MathUtils.TWO_PI)));
            assertEquals(Location.OUTSIDE,  aMb.checkPoint(new S1Point(4.9 + k * MathUtils.TWO_PI)));
            assertEquals(Location.BOUNDARY, aMb.checkPoint(new S1Point(5.0 + k * MathUtils.TWO_PI)));
            assertEquals(Location.INSIDE,   aMb.checkPoint(new S1Point(5.1 + k * MathUtils.TWO_PI)));
            assertEquals(Location.INSIDE,   aMb.checkPoint(new S1Point(5.9 + k * MathUtils.TWO_PI)));
            assertEquals(Location.BOUNDARY, aMb.checkPoint(new S1Point(6.0 + k * MathUtils.TWO_PI)));
            assertEquals(Location.OUTSIDE,  aMb.checkPoint(new S1Point(6.1 + k * MathUtils.TWO_PI)));
            assertEquals(Location.OUTSIDE,  aMb.checkPoint(new S1Point(6.2 + k * MathUtils.TWO_PI)));
        }

        List<Arc> aMbList = aMb.asList();
        assertEquals(2,   aMbList.size());
        assertEquals(1.0, aMbList.get(0).getInf(), 1.0e-10);
        assertEquals(3.0, aMbList.get(0).getSup(), 1.0e-10);
        assertEquals(5.0, aMbList.get(1).getInf(), 1.0e-10);
        assertEquals(6.0, aMbList.get(1).getSup(), 1.0e-10);
        assertEquals(Location.INSIDE, aMb.checkPoint(aMb.getInteriorPoint()));


    }

    @Test
    void testIntersection() {

        ArcsSet a   = (ArcsSet) new RegionFactory<Sphere1D, S1Point, LimitAngle, SubLimitAngle>().
                       union(new ArcsSet(1.0, 3.0, 1.0e-10), new ArcsSet(5.0, 6.0, 1.0e-10));
        List<Arc> aList = a.asList();
        assertEquals(2,   aList.size());
        assertEquals(1.0, aList.get(0).getInf(), 1.0e-10);
        assertEquals(3.0, aList.get(0).getSup(), 1.0e-10);
        assertEquals(5.0, aList.get(1).getInf(), 1.0e-10);
        assertEquals(6.0, aList.get(1).getSup(), 1.0e-10);

        ArcsSet b   = new ArcsSet(0.0, 5.5, 1.0e-10);
        List<Arc> bList = b.asList();
        assertEquals(1,   bList.size());
        assertEquals(0.0, bList.get(0).getInf(), 1.0e-10);
        assertEquals(5.5, bList.get(0).getSup(), 1.0e-10);

        ArcsSet aMb = (ArcsSet) new RegionFactory<Sphere1D, S1Point, LimitAngle, SubLimitAngle>().
                intersection(a, b);
        for (int k = -2; k < 3; ++k) {
            assertEquals(Location.OUTSIDE,  aMb.checkPoint(new S1Point(0.0 + k * MathUtils.TWO_PI)));
            assertEquals(Location.BOUNDARY, aMb.checkPoint(new S1Point(1.0 + k * MathUtils.TWO_PI)));
            assertEquals(Location.INSIDE,   aMb.checkPoint(new S1Point(1.1 + k * MathUtils.TWO_PI)));
            assertEquals(Location.INSIDE,   aMb.checkPoint(new S1Point(2.9 + k * MathUtils.TWO_PI)));
            assertEquals(Location.BOUNDARY, aMb.checkPoint(new S1Point(3.0 + k * MathUtils.TWO_PI)));
            assertEquals(Location.OUTSIDE,  aMb.checkPoint(new S1Point(3.1 + k * MathUtils.TWO_PI)));
            assertEquals(Location.OUTSIDE,  aMb.checkPoint(new S1Point(4.9 + k * MathUtils.TWO_PI)));
            assertEquals(Location.BOUNDARY, aMb.checkPoint(new S1Point(5.0 + k * MathUtils.TWO_PI)));
            assertEquals(Location.INSIDE,   aMb.checkPoint(new S1Point(5.1 + k * MathUtils.TWO_PI)));
            assertEquals(Location.INSIDE,   aMb.checkPoint(new S1Point(5.4 + k * MathUtils.TWO_PI)));
            assertEquals(Location.BOUNDARY, aMb.checkPoint(new S1Point(5.5 + k * MathUtils.TWO_PI)));
            assertEquals(Location.OUTSIDE,  aMb.checkPoint(new S1Point(5.6 + k * MathUtils.TWO_PI)));
            assertEquals(Location.OUTSIDE,  aMb.checkPoint(new S1Point(6.2 + k * MathUtils.TWO_PI)));
        }

        List<Arc> aMbList = aMb.asList();
        assertEquals(2,   aMbList.size());
        assertEquals(1.0, aMbList.get(0).getInf(), 1.0e-10);
        assertEquals(3.0, aMbList.get(0).getSup(), 1.0e-10);
        assertEquals(5.0, aMbList.get(1).getInf(), 1.0e-10);
        assertEquals(5.5, aMbList.get(1).getSup(), 1.0e-10);
        assertEquals(Location.INSIDE, aMb.checkPoint(aMb.getInteriorPoint()));

    }

    @Test
    void testMultiple() {
        RegionFactory<Sphere1D, S1Point, LimitAngle, SubLimitAngle> factory = new RegionFactory<>();
        ArcsSet set = (ArcsSet)
        factory.intersection(factory.union(factory.difference(new ArcsSet(1.0, 6.0, 1.0e-10),
                                                              new ArcsSet(3.0, 5.0, 1.0e-10)),
                                                              new ArcsSet(0.5, 2.0, 1.0e-10)),
                                                              new ArcsSet(0.0, 5.5, 1.0e-10));
        assertEquals(3.0, set.getSize(), 1.0e-10);
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new S1Point(0.0)));
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new S1Point(4.0)));
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new S1Point(6.0)));
        assertEquals(Region.Location.INSIDE,   set.checkPoint(new S1Point(1.2)));
        assertEquals(Region.Location.INSIDE,   set.checkPoint(new S1Point(5.25)));
        assertEquals(Region.Location.BOUNDARY, set.checkPoint(new S1Point(0.5)));
        assertEquals(Region.Location.BOUNDARY, set.checkPoint(new S1Point(3.0)));
        assertEquals(Region.Location.BOUNDARY, set.checkPoint(new S1Point(5.0)));
        assertEquals(Region.Location.BOUNDARY, set.checkPoint(new S1Point(5.5)));

        List<Arc> list = set.asList();
        assertEquals(2, list.size());
        assertEquals( 0.5, list.get(0).getInf(), 1.0e-10);
        assertEquals( 3.0, list.get(0).getSup(), 1.0e-10);
        assertEquals( 5.0, list.get(1).getInf(), 1.0e-10);
        assertEquals( 5.5, list.get(1).getSup(), 1.0e-10);
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));

    }

    @Test
    void testSinglePoint() {
        ArcsSet set = new ArcsSet(1.0, FastMath.nextAfter(1.0, Double.POSITIVE_INFINITY), 1.0e-10);
        assertEquals(2 * Precision.EPSILON, set.getSize(), Precision.SAFE_MIN);
    }

    @Test
    void testIteration() {
        ArcsSet set = (ArcsSet) new RegionFactory<Sphere1D, S1Point, LimitAngle, SubLimitAngle>().
                difference(new ArcsSet(1.0, 6.0, 1.0e-10), new ArcsSet(3.0, 5.0, 1.0e-10));
        Iterator<double[]> iterator = set.iterator();
        try {
            iterator.remove();
            fail("an exception should have been thrown");
        } catch (UnsupportedOperationException uoe) {
            // expected
        }

        assertTrue(iterator.hasNext());
        double[] a0 = iterator.next();
        assertEquals(2, a0.length);
        assertEquals(1.0, a0[0], 1.0e-10);
        assertEquals(3.0, a0[1], 1.0e-10);

        assertTrue(iterator.hasNext());
        double[] a1 = iterator.next();
        assertEquals(2, a1.length);
        assertEquals(5.0, a1[0], 1.0e-10);
        assertEquals(6.0, a1[1], 1.0e-10);

        assertFalse(iterator.hasNext());
        try {
            iterator.next();
            fail("an exception should have been thrown");
        } catch (NoSuchElementException nsee) {
            // expected
        }

    }

    @Test
    void testEmptyTree() {
        assertEquals(MathUtils.TWO_PI,
                     new ArcsSet(new BSPTree<>(Boolean.TRUE), 1.0e-10).getSize(),
                     1.0e-10);
    }

    @Test
    void testShiftedAngles() {
        for (int k = -2; k < 3; ++k) {
            SubLimitAngle l1  = new LimitAngle(new S1Point(1.0 + k * MathUtils.TWO_PI), false, 1.0e-10).wholeHyperplane();
            SubLimitAngle l2  = new LimitAngle(new S1Point(1.5 + k * MathUtils.TWO_PI), true,  1.0e-10).wholeHyperplane();
            ArcsSet set = new ArcsSet(new BSPTree<>(l1,
                                                    new BSPTree<>(Boolean.FALSE),
                                                    new BSPTree<>(l2,
                                                                  new BSPTree<>(Boolean.FALSE),
                                                                  new BSPTree<>(Boolean.TRUE),
                                                                  null),
                                                    null),
                                      1.0e-10);
            for (double alpha = 1.0e-6; alpha < MathUtils.TWO_PI; alpha += 0.001) {
                if (alpha < 1 || alpha > 1.5) {
                    assertEquals(Location.OUTSIDE, set.checkPoint(new S1Point(alpha)));
                } else {
                    assertEquals(Location.INSIDE,  set.checkPoint(new S1Point(alpha)));
                }
            }
        }

    }

    @Test
    void testInconsistentState() {
        assertThrows(ArcsSet.InconsistentStateAt2PiWrapping.class, () -> {
            SubLimitAngle l1 = new LimitAngle(new S1Point(1.0), false, 1.0e-10).wholeHyperplane();
            SubLimitAngle l2 = new LimitAngle(new S1Point(2.0), true, 1.0e-10).wholeHyperplane();
            SubLimitAngle l3 = new LimitAngle(new S1Point(3.0), false, 1.0e-10).wholeHyperplane();
            new ArcsSet(new BSPTree<>(l1,
                                      new BSPTree<>(Boolean.FALSE),
                                      new BSPTree<>(l2,
                                                    new BSPTree<>(l3, new BSPTree<>(Boolean.FALSE), new BSPTree<>(Boolean.TRUE), null),
                                                    new BSPTree<>(Boolean.TRUE),
                                                    null),
                                      null),
                        1.0e-10);
        });
    }

    @Test
    void testSide() {
        ArcsSet set = (ArcsSet) new RegionFactory<Sphere1D, S1Point, LimitAngle, SubLimitAngle>().
                difference(new ArcsSet(1.0, 6.0, 1.0e-10), new ArcsSet(3.0, 5.0, 1.0e-10));
        for (int k = -2; k < 3; ++k) {
            assertEquals(Side.MINUS, set.split(new Arc(0.5 + k * MathUtils.TWO_PI,
                                                              6.1 + k * MathUtils.TWO_PI,
                                                              set.getTolerance())).getSide());
            assertEquals(Side.PLUS,  set.split(new Arc(0.5 + k * MathUtils.TWO_PI,
                                                              0.8 + k * MathUtils.TWO_PI,
                                                              set.getTolerance())).getSide());
            assertEquals(Side.PLUS,  set.split(new Arc(6.2 + k * MathUtils.TWO_PI,
                                                              6.3 + k * MathUtils.TWO_PI,
                                                              set.getTolerance())).getSide());
            assertEquals(Side.PLUS,  set.split(new Arc(3.5 + k * MathUtils.TWO_PI,
                                                              4.5 + k * MathUtils.TWO_PI,
                                                              set.getTolerance())).getSide());
            assertEquals(Side.BOTH,  set.split(new Arc(2.9 + k * MathUtils.TWO_PI,
                                                              4.5 + k * MathUtils.TWO_PI,
                                                              set.getTolerance())).getSide());
            assertEquals(Side.BOTH,  set.split(new Arc(0.5 + k * MathUtils.TWO_PI,
                                                              1.2 + k * MathUtils.TWO_PI,
                                                              set.getTolerance())).getSide());
            assertEquals(Side.BOTH,  set.split(new Arc(0.5 + k * MathUtils.TWO_PI,
                                                              5.9 + k * MathUtils.TWO_PI,
                                                              set.getTolerance())).getSide());
        }
    }

    @Test
    void testSideEmbedded() {

        ArcsSet s35 = new ArcsSet(3.0, 5.0, 1.0e-10);
        ArcsSet s16 = new ArcsSet(1.0, 6.0, 1.0e-10);

        assertEquals(Side.BOTH,  s16.split(new Arc(3.0, 5.0, 1.0e-10)).getSide());
        assertEquals(Side.BOTH,  s16.split(new Arc(5.0, 3.0 + MathUtils.TWO_PI, 1.0e-10)).getSide());
        assertEquals(Side.MINUS, s35.split(new Arc(1.0, 6.0, 1.0e-10)).getSide());
        assertEquals(Side.PLUS,  s35.split(new Arc(6.0, 1.0 + MathUtils.TWO_PI, 1.0e-10)).getSide());

    }

    @Test
    void testSideOverlapping() {
        ArcsSet s35 = new ArcsSet(3.0, 5.0, 1.0e-10);
        ArcsSet s46 = new ArcsSet(4.0, 6.0, 1.0e-10);

        assertEquals(Side.BOTH,  s46.split(new Arc(3.0, 5.0, 1.0e-10)).getSide());
        assertEquals(Side.BOTH,  s46.split(new Arc(5.0, 3.0 + MathUtils.TWO_PI, 1.0e-10)).getSide());
        assertEquals(Side.BOTH, s35.split(new Arc(4.0, 6.0, 1.0e-10)).getSide());
        assertEquals(Side.BOTH,  s35.split(new Arc(6.0, 4.0 + MathUtils.TWO_PI, 1.0e-10)).getSide());
    }

    @Test
    void testSideHyper() {
        ArcsSet sub = (ArcsSet) new RegionFactory<Sphere1D, S1Point, LimitAngle, SubLimitAngle>().
                getComplement(new ArcsSet(1.0e-10));
        assertTrue(sub.isEmpty());
        assertEquals(Side.HYPER,  sub.split(new Arc(2.0, 3.0, 1.0e-10)).getSide());
    }

    @Test
    void testSplitEmbedded() {

        ArcsSet s35 = new ArcsSet(3.0, 5.0, 1.0e-10);
        ArcsSet s16 = new ArcsSet(1.0, 6.0, 1.0e-10);

        ArcsSet.Split split1 = s16.split(new Arc(3.0, 5.0, 1.0e-10));
        ArcsSet split1Plus  = split1.getPlus();
        ArcsSet split1Minus = split1.getMinus();
        assertEquals(3.0, split1Plus.getSize(), 1.0e-10);
        assertEquals(2,   split1Plus.asList().size());
        assertEquals(1.0, split1Plus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(3.0, split1Plus.asList().get(0).getSup(), 1.0e-10);
        assertEquals(5.0, split1Plus.asList().get(1).getInf(), 1.0e-10);
        assertEquals(6.0, split1Plus.asList().get(1).getSup(), 1.0e-10);
        assertEquals(2.0, split1Minus.getSize(), 1.0e-10);
        assertEquals(1,   split1Minus.asList().size());
        assertEquals(3.0, split1Minus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(5.0, split1Minus.asList().get(0).getSup(), 1.0e-10);

        ArcsSet.Split split2 = s16.split(new Arc(5.0, 3.0 + MathUtils.TWO_PI, 1.0e-10));
        ArcsSet split2Plus  = split2.getPlus();
        ArcsSet split2Minus = split2.getMinus();
        assertEquals(2.0, split2Plus.getSize(), 1.0e-10);
        assertEquals(1,   split2Plus.asList().size());
        assertEquals(3.0, split2Plus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(5.0, split2Plus.asList().get(0).getSup(), 1.0e-10);
        assertEquals(3.0, split2Minus.getSize(), 1.0e-10);
        assertEquals(2,   split2Minus.asList().size());
        assertEquals(1.0, split2Minus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(3.0, split2Minus.asList().get(0).getSup(), 1.0e-10);
        assertEquals(5.0, split2Minus.asList().get(1).getInf(), 1.0e-10);
        assertEquals(6.0, split2Minus.asList().get(1).getSup(), 1.0e-10);

        ArcsSet.Split split3 = s35.split(new Arc(1.0, 6.0, 1.0e-10));
        ArcsSet split3Plus  = split3.getPlus();
        ArcsSet split3Minus = split3.getMinus();
        assertNull(split3Plus);
        assertEquals(2.0, split3Minus.getSize(), 1.0e-10);
        assertEquals(1,   split3Minus.asList().size());
        assertEquals(3.0, split3Minus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(5.0, split3Minus.asList().get(0).getSup(), 1.0e-10);

        ArcsSet.Split split4 = s35.split(new Arc(6.0, 1.0 + MathUtils.TWO_PI, 1.0e-10));
        ArcsSet split4Plus  = split4.getPlus();
        ArcsSet split4Minus = split4.getMinus();
        assertEquals(2.0, split4Plus.getSize(), 1.0e-10);
        assertEquals(1,   split4Plus.asList().size());
        assertEquals(3.0, split4Plus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(5.0, split4Plus.asList().get(0).getSup(), 1.0e-10);
        assertNull(split4Minus);

    }

    @Test
    void testSplitOverlapping() {

        ArcsSet s35 = new ArcsSet(3.0, 5.0, 1.0e-10);
        ArcsSet s46 = new ArcsSet(4.0, 6.0, 1.0e-10);

        ArcsSet.Split split1 = s46.split(new Arc(3.0, 5.0, 1.0e-10));
        ArcsSet split1Plus  = split1.getPlus();
        ArcsSet split1Minus = split1.getMinus();
        assertEquals(1.0, split1Plus.getSize(), 1.0e-10);
        assertEquals(1,   split1Plus.asList().size());
        assertEquals(5.0, split1Plus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(6.0, split1Plus.asList().get(0).getSup(), 1.0e-10);
        assertEquals(1.0, split1Minus.getSize(), 1.0e-10);
        assertEquals(1,   split1Minus.asList().size());
        assertEquals(4.0, split1Minus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(5.0, split1Minus.asList().get(0).getSup(), 1.0e-10);

        ArcsSet.Split split2 = s46.split(new Arc(5.0, 3.0 + MathUtils.TWO_PI, 1.0e-10));
        ArcsSet split2Plus  = split2.getPlus();
        ArcsSet split2Minus = split2.getMinus();
        assertEquals(1.0, split2Plus.getSize(), 1.0e-10);
        assertEquals(1,   split2Plus.asList().size());
        assertEquals(4.0, split2Plus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(5.0, split2Plus.asList().get(0).getSup(), 1.0e-10);
        assertEquals(1.0, split2Minus.getSize(), 1.0e-10);
        assertEquals(1,   split2Minus.asList().size());
        assertEquals(5.0, split2Minus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(6.0, split2Minus.asList().get(0).getSup(), 1.0e-10);

        ArcsSet.Split split3 = s35.split(new Arc(4.0, 6.0, 1.0e-10));
        ArcsSet split3Plus  = split3.getPlus();
        ArcsSet split3Minus = split3.getMinus();
        assertEquals(1.0, split3Plus.getSize(), 1.0e-10);
        assertEquals(1,   split3Plus.asList().size());
        assertEquals(3.0, split3Plus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(4.0, split3Plus.asList().get(0).getSup(), 1.0e-10);
        assertEquals(1.0, split3Minus.getSize(), 1.0e-10);
        assertEquals(1,   split3Minus.asList().size());
        assertEquals(4.0, split3Minus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(5.0, split3Minus.asList().get(0).getSup(), 1.0e-10);

        ArcsSet.Split split4 = s35.split(new Arc(6.0, 4.0 + MathUtils.TWO_PI, 1.0e-10));
        ArcsSet split4Plus  = split4.getPlus();
        ArcsSet split4Minus = split4.getMinus();
        assertEquals(1.0, split4Plus.getSize(), 1.0e-10);
        assertEquals(1,   split4Plus.asList().size());
        assertEquals(4.0, split4Plus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(5.0, split4Plus.asList().get(0).getSup(), 1.0e-10);
        assertEquals(1.0, split4Minus.getSize(), 1.0e-10);
        assertEquals(1,   split4Minus.asList().size());
        assertEquals(3.0, split4Minus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(4.0, split4Minus.asList().get(0).getSup(), 1.0e-10);

    }

    @Test
    void testFarSplit() {
        ArcsSet set = new ArcsSet(FastMath.PI, 2.5 * FastMath.PI, 1.0e-10);
        ArcsSet.Split split = set.split(new Arc(MathUtils.SEMI_PI, 1.5 * FastMath.PI, 1.0e-10));
        ArcsSet splitPlus  = split.getPlus();
        ArcsSet splitMinus = split.getMinus();
        assertEquals(1,   splitMinus.asList().size());
        assertEquals(      FastMath.PI, splitMinus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(1.5 * FastMath.PI, splitMinus.asList().get(0).getSup(), 1.0e-10);
        assertEquals(MathUtils.SEMI_PI, splitMinus.getSize(), 1.0e-10);
        assertEquals(1,   splitPlus.asList().size());
        assertEquals(1.5 * FastMath.PI, splitPlus.asList().get(0).getInf(), 1.0e-10);
        assertEquals(2.5 * FastMath.PI, splitPlus.asList().get(0).getSup(), 1.0e-10);
        assertEquals(      FastMath.PI, splitPlus.getSize(), 1.0e-10);

    }

    @Test
    void testSplitWithinEpsilon() {
        double epsilon = 1.0e-10;
        double a = 6.25;
        double b = a - 0.5 * epsilon;
        ArcsSet set = new ArcsSet(a - 1, a, epsilon);
        Arc arc = new Arc(b, b + FastMath.PI, epsilon);
        ArcsSet.Split split = set.split(arc);
        assertEquals(set.getSize(), split.getPlus().getSize(),  epsilon);
        assertNull(split.getMinus());
    }

    @Test
    void testSideSplitConsistency() {
        double  epsilon = 1.0e-6;
        double  a       = 4.725;
        ArcsSet set     = new ArcsSet(a, a + 0.5, epsilon);
        Arc     arc     = new Arc(a + 0.5 * epsilon, a + 1, epsilon);
        ArcsSet.Split split = set.split(arc);
        assertNotNull(split.getMinus());
        assertNull(split.getPlus());
        assertEquals(Side.MINUS, set.split(arc).getSide());
    }

    @Test
    void testShuffledTreeNonRepresentable() {
        doTestShuffledTree(FastMath.toRadians( 85.0), FastMath.toRadians( 95.0),
                           FastMath.toRadians(265.0), FastMath.toRadians(275.0),
                           1.0e-10);
    }

    @Test
    void testShuffledTreeRepresentable() {
        doTestShuffledTree(1.0, 2.0,
                           4.0, 5.0,
                           1.0e-10);
    }

    private void doTestShuffledTree(final double a0, final double a1,
                                    final double a2, final double a3,
                                    final double tol) {


        // intervals set [ a0 ; a1 ] U [ a2 ; a3 ]
        final ArcsSet setA =
                        new ArcsSet(new BSPTree<>(new LimitAngle(new S1Point(a0), false, tol).wholeHyperplane(),
                                                       new BSPTree<>(Boolean.FALSE),
                                                       new BSPTree<>(new LimitAngle(new S1Point(a1), true, tol).wholeHyperplane(),
                                                                     new BSPTree<>(new LimitAngle(new S1Point(a2), false, tol).wholeHyperplane(),
                                                                                   new BSPTree<>(Boolean.FALSE),
                                                                                   new BSPTree<>(new LimitAngle(new S1Point(a3), true, tol).wholeHyperplane(),
                                                                                                 new BSPTree<>(Boolean.FALSE),
                                                                                                 new BSPTree<>(Boolean.TRUE),
                                                                                                 null),
                                                                                   null),
                                                                     new BSPTree<>(Boolean.TRUE),
                                                                     null),
                                                       null),
                                         1.0e-10);
        assertEquals((a1 - a0) + (a3 - a2), setA.getSize(), 1.0e-10);
        assertEquals(2, setA.asList().size());
        assertEquals(a0, setA.asList().get(0).getInf(), 1.0e-15);
        assertEquals(a1, setA.asList().get(0).getSup(), 1.0e-15);
        assertEquals(a2, setA.asList().get(1).getInf(), 1.0e-15);
        assertEquals(a3, setA.asList().get(1).getSup(), 1.0e-15);

        // same intervals set [ a0 ; a1 ] U [ a2 ; a3 ], but with a different tree organization
        final ArcsSet setB =
                        new ArcsSet(new BSPTree<>(new LimitAngle(new S1Point(a2), false, tol).wholeHyperplane(),
                                                       new BSPTree<>(new LimitAngle(new S1Point(a0), false, tol).wholeHyperplane(),
                                                                     new BSPTree<>(Boolean.FALSE),
                                                                     new BSPTree<>(new LimitAngle(new S1Point(a1), true, tol).wholeHyperplane(),
                                                                                   new BSPTree<>(Boolean.FALSE),
                                                                                   new BSPTree<>(Boolean.TRUE),
                                                                                   null),
                                                                     null),
                                                       new BSPTree<>(new LimitAngle(new S1Point(a3), true, tol).wholeHyperplane(),
                                                                     new BSPTree<>(Boolean.FALSE),
                                                                     new BSPTree<>(Boolean.TRUE),
                                                                     null),
                                                       null),
                                         1.0e-10);
        assertEquals((a1 - a0) + (a3 - a2), setB.getSize(), 1.0e-10);
        assertEquals(2, setB.asList().size());
        assertEquals(a0, setB.asList().get(0).getInf(), 1.0e-15);
        assertEquals(a1, setB.asList().get(0).getSup(), 1.0e-15);
        assertEquals(a2, setB.asList().get(1).getInf(), 1.0e-15);
        assertEquals(a3, setB.asList().get(1).getSup(), 1.0e-15);

        final ArcsSet intersection = (ArcsSet) new RegionFactory<Sphere1D, S1Point, LimitAngle, SubLimitAngle>().
                intersection(setA, setB);
        assertEquals((a1 - a0) + (a3 - a2), intersection.getSize(), 1.0e-10);
        assertEquals(2, intersection.asList().size());
        assertEquals(a0, intersection.asList().get(0).getInf(), 1.0e-15);
        assertEquals(a1, intersection.asList().get(0).getSup(), 1.0e-15);
        assertEquals(a2, intersection.asList().get(1).getInf(), 1.0e-15);
        assertEquals(a3, intersection.asList().get(1).getSup(), 1.0e-15);

    }

}
