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

import org.hipparchus.geometry.partitioning.BSPTree;
import org.hipparchus.geometry.partitioning.Region;
import org.hipparchus.geometry.partitioning.RegionFactory;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class IntervalsSetTest {

    @Test
    void testInterval() {
        IntervalsSet set = new IntervalsSet(Arrays.asList(new OrientedPoint(new Vector1D(2.3), false, 1.0e-10).wholeHyperplane(),
                                                          new OrientedPoint(new Vector1D(5.7), true, 1.0e-10).wholeHyperplane()),
                                            1.0e-10);
        assertEquals(3.4, set.getSize(), 1.0e-10);
        assertEquals(4.0, set.getBarycenter().getX(), 1.0e-10);
        assertEquals(Region.Location.BOUNDARY, set.checkPoint(new Vector1D(2.3)));
        assertEquals(Region.Location.BOUNDARY, set.checkPoint(new Vector1D(5.7)));
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new Vector1D(1.2)));
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new Vector1D(8.7)));
        assertEquals(Region.Location.INSIDE,   set.checkPoint(new Vector1D(3.0)));
        assertEquals(2.3, set.getInf(), 1.0e-10);
        assertEquals(5.7, set.getSup(), 1.0e-10);
        assertEquals(Region.Location.INSIDE, set.checkPoint(set.getInteriorPoint()));

        OrientedPoint op = set.getTree(false).getCut().getHyperplane();
        assertEquals(0.0, op.emptyHyperplane().getSize(), 1.0e-10);
        assertEquals(1.0e-10, op.getTolerance(), 1.0e-20);
        assertTrue(Double.isInfinite(op.wholeSpace().getSize()));
        assertEquals(2.3, op.getLocation().getX(), 1.0e-10);
        assertEquals(-0.7, op.getOffset(new Vector1D(3.0)), 1.0e-10);
        assertSame(op.getLocation(), op.project(new Vector1D(3.0)));

        // beware that after changing point orientation, we completely mess up interval tree
        // so we cannot use it anymore
        op.revertSelf();
        assertEquals(+0.7, op.getOffset(new Vector1D(3.0)), 1.0e-10);

    }

    @Test
    void testNoBoundaries() {
        IntervalsSet set = new IntervalsSet(new ArrayList<>(), 1.0e-10);
        assertEquals(Region.Location.INSIDE, set.checkPoint(new Vector1D(-Double.MAX_VALUE)));
        assertEquals(Region.Location.INSIDE, set.checkPoint(new Vector1D(+Double.MAX_VALUE)));
        assertTrue(Double.isInfinite(set.getSize()));
        assertEquals(Region.Location.INSIDE, set.checkPoint(set.getInteriorPoint()));
    }

    @Test
    void testInfinite() {
        IntervalsSet set = new IntervalsSet(9.0, Double.POSITIVE_INFINITY, 1.0e-10);
        assertEquals(Region.Location.BOUNDARY, set.checkPoint(new Vector1D(9.0)));
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new Vector1D(8.4)));
        for (double e = 1.0; e <= 6.0; e += 1.0) {
            assertEquals(Region.Location.INSIDE,
                                set.checkPoint(new Vector1D(FastMath.pow(10.0, e))));
        }
        assertTrue(Double.isInfinite(set.getSize()));
        assertEquals(9.0, set.getInf(), 1.0e-10);
        assertTrue(Double.isInfinite(set.getSup()));

        set = (IntervalsSet) new RegionFactory<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint>().getComplement(set);
        assertEquals(9.0, set.getSup(), 1.0e-10);
        assertTrue(Double.isInfinite(set.getInf()));
        assertEquals(Region.Location.INSIDE, set.checkPoint(set.getInteriorPoint()));
    }

    @Test
    void testMultiple() {
        RegionFactory<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint> factory = new RegionFactory<>();
        IntervalsSet set = (IntervalsSet) factory.intersection(factory.union(factory.difference(new IntervalsSet(1.0, 6.0, 1.0e-10),
                                                                                                new IntervalsSet(3.0, 5.0, 1.0e-10)),
                                                                             new IntervalsSet(9.0, Double.POSITIVE_INFINITY, 1.0e-10)),
                                                               new IntervalsSet(Double.NEGATIVE_INFINITY, 11.0, 1.0e-10));
        assertEquals(5.0, set.getSize(), 1.0e-10);
        assertEquals(5.9, set.getBarycenter().getX(), 1.0e-10);
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new Vector1D(0.0)));
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new Vector1D(4.0)));
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new Vector1D(8.0)));
        assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new Vector1D(12.0)));
        assertEquals(Region.Location.INSIDE,   set.checkPoint(new Vector1D(1.2)));
        assertEquals(Region.Location.INSIDE,   set.checkPoint(new Vector1D(5.9)));
        assertEquals(Region.Location.INSIDE,   set.checkPoint(new Vector1D(9.01)));
        assertEquals(Region.Location.BOUNDARY, set.checkPoint(new Vector1D(5.0)));
        assertEquals(Region.Location.BOUNDARY, set.checkPoint(new Vector1D(11.0)));
        assertEquals( 1.0, set.getInf(), 1.0e-10);
        assertEquals(11.0, set.getSup(), 1.0e-10);

        List<Interval> list = set.asList();
        assertEquals(3, list.size());
        assertEquals( 1.0, list.get(0).getInf(), 1.0e-10);
        assertEquals( 3.0, list.get(0).getSup(), 1.0e-10);
        assertEquals( 5.0, list.get(1).getInf(), 1.0e-10);
        assertEquals( 6.0, list.get(1).getSup(), 1.0e-10);
        assertEquals( 9.0, list.get(2).getInf(), 1.0e-10);
        assertEquals(11.0, list.get(2).getSup(), 1.0e-10);

        assertEquals(Region.Location.INSIDE, set.checkPoint(set.getInteriorPoint()));

    }

    @Test
    void testInteriorPoint() {

        final IntervalsSet wholeSpace = new IntervalsSet(1.0e-10);
        final RegionFactory<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint> factory = new RegionFactory<>();
        final Region<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint> empty = factory.getComplement(wholeSpace);
        assertNull(empty.getInteriorPoint());

        final IntervalsSet is1 = new IntervalsSet(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0);
        final Vector1D     v1  = is1.getInteriorPoint();
        assertEquals(0.0, v1.getX(), 1.0e-10);
        assertEquals(Region.Location.INSIDE, is1.checkPoint(v1));

        final IntervalsSet is2 = new IntervalsSet(Double.NEGATIVE_INFINITY, 100.0, 1.0);
        final Vector1D     v2  = is2.getInteriorPoint();
        assertEquals(-999999900.0, v2.getX(), 1.0e-10);
        assertEquals(Region.Location.INSIDE, is2.checkPoint(v2));

        final IntervalsSet is3 = new IntervalsSet(-100.0, Double.POSITIVE_INFINITY, 1.0);
        final Vector1D     v3  = is3.getInteriorPoint();
        assertEquals(999999900.0, v3.getX(), 1.0e-10);
        assertEquals(Region.Location.INSIDE, is3.checkPoint(v3));

        final IntervalsSet is4 = new IntervalsSet(11.0, 13.0, 1.0);
        final Vector1D     v4  = is4.getInteriorPoint();
        assertEquals(12.0, v4.getX(), 1.0e-10);
        assertEquals(Region.Location.INSIDE, is4.checkPoint(v4));

        final IntervalsSet is5  = (IntervalsSet) factory.union(factory.union(new IntervalsSet(11.0, 13.0, 1.0),
                                                                             new IntervalsSet( 9.0,  9.5, 1.0)),
                                                               new IntervalsSet(100.0, 120.0, 1.0));
        final Vector1D     v5  = is5.getInteriorPoint();
        assertEquals(110.0, v5.getX(), 1.0e-10);
        assertEquals(Region.Location.INSIDE, is5.checkPoint(v5));

    }

    @Test
    void testSinglePoint() {
        IntervalsSet set = new IntervalsSet(1.0, 1.0, 1.0e-10);
        assertEquals(0.0, set.getSize(), Precision.SAFE_MIN);
        assertEquals(1.0, set.getBarycenter().getX(), Precision.EPSILON);
        try {
            set.iterator().remove();
            fail("an exception should have been thrown");
        } catch (UnsupportedOperationException uoe) {
            // expected
        }
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
        final IntervalsSet setA =
                        new IntervalsSet(new BSPTree<>(new OrientedPoint(new Vector1D(a0), false, tol).wholeHyperplane(),
                                                       new BSPTree<>(Boolean.FALSE),
                                                       new BSPTree<>(new OrientedPoint(new Vector1D(a1), true, tol).wholeHyperplane(),
                                                                     new BSPTree<>(new OrientedPoint(new Vector1D(a2), false, tol).wholeHyperplane(),
                                                                                   new BSPTree<>(Boolean.FALSE),
                                                                                   new BSPTree<>(new OrientedPoint(new Vector1D(a3), true, tol).wholeHyperplane(),
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
        final IntervalsSet setB =
                        new IntervalsSet(new BSPTree<>(new OrientedPoint(new Vector1D(a2), false, tol).wholeHyperplane(),
                                                       new BSPTree<>(new OrientedPoint(new Vector1D(a0), false, tol).wholeHyperplane(),
                                                                     new BSPTree<>(Boolean.FALSE),
                                                                     new BSPTree<>(new OrientedPoint(new Vector1D(a1), true, tol).wholeHyperplane(),
                                                                                   new BSPTree<>(Boolean.FALSE),
                                                                                   new BSPTree<>(Boolean.TRUE),
                                                                                   null),
                                                                     null),
                                                       new BSPTree<>(new OrientedPoint(new Vector1D(a3), true, tol).wholeHyperplane(),
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

        final IntervalsSet intersection = (IntervalsSet) new RegionFactory<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint>().
                                          intersection(setA, setB);
        assertEquals((a1 - a0) + (a3 - a2), intersection.getSize(), 1.0e-10);
        assertEquals(2, intersection.asList().size());
        assertEquals(a0, intersection.asList().get(0).getInf(), 1.0e-15);
        assertEquals(a1, intersection.asList().get(0).getSup(), 1.0e-15);
        assertEquals(a2, intersection.asList().get(1).getInf(), 1.0e-15);
        assertEquals(a3, intersection.asList().get(1).getSup(), 1.0e-15);

    }

}
