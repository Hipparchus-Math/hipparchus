/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.geometry.euclidean.oned;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class OrientedPointTest {

    @Test
    void testOrientedPoint() {
        OrientedPoint pt = new OrientedPoint(new Vector1D(1.0), true, 1.0e-10);
        assertEquals(1.0, pt.getLocation().getX(), 1.0e-10);
        assertTrue(pt.isDirect());
        assertEquals(1.0e-10, pt.getTolerance(), 1.0e-10);
        assertSame(pt, pt.copySelf());
    }

    @Test
    void testGetOffset() {
        OrientedPoint direct = new OrientedPoint(new Vector1D(1.0), true, 1.0e-10);
        assertEquals(0.5, direct.getOffset(new Vector1D(1.5)), 1.0e-10);
        assertEquals(-0.5, direct.getOffset(new Vector1D(0.5)), 1.0e-10);

        OrientedPoint indirect = new OrientedPoint(new Vector1D(1.0), false, 1.0e-10);
        assertEquals(-0.5, indirect.getOffset(new Vector1D(1.5)), 1.0e-10);
        assertEquals(0.5, indirect.getOffset(new Vector1D(0.5)), 1.0e-10);
    }

    @Test
    void testMoveToOffset() {
        OrientedPoint direct = new OrientedPoint(new Vector1D(1.0), true, 1.0e-10);
        assertEquals(1.5, direct.moveToOffset(new Vector1D(10.0), 0.5).getX(), 1.0e-10);
        assertEquals(0.5, direct.moveToOffset(new Vector1D(10.0), -0.5).getX(), 1.0e-10);

        OrientedPoint indirect = new OrientedPoint(new Vector1D(1.0), false, 1.0e-10);
        assertEquals(0.5, indirect.moveToOffset(new Vector1D(10.0), 0.5).getX(), 1.0e-10);
        assertEquals(1.5, indirect.moveToOffset(new Vector1D(10.0), -0.5).getX(), 1.0e-10);
    }

    @Test
    void testArbitraryPoint() {
        Vector1D location = new Vector1D(1.0);
        OrientedPoint pt = new OrientedPoint(location, true, 1.0e-10);
        assertSame(location, pt.arbitraryPoint());
    }

    @Test
    void testWholeHyperplane() {
        OrientedPoint pt = new OrientedPoint(new Vector1D(1.0), true, 1.0e-10);
        SubOrientedPoint sub = pt.wholeHyperplane();
        assertNotNull(sub);
        assertSame(pt, sub.getHyperplane());
    }

    @Test
    void testEmptyHyperplane() {
        OrientedPoint pt = new OrientedPoint(new Vector1D(1.0), true, 1.0e-10);
        SubOrientedPoint sub = pt.emptyHyperplane();
        assertNotNull(sub);
        assertSame(pt, sub.getHyperplane());
    }

    @Test
    void testWholeSpace() {
        OrientedPoint pt = new OrientedPoint(new Vector1D(1.0), true, 1.0e-10);
        IntervalsSet set = pt.wholeSpace();
        assertTrue(Double.isInfinite(set.getSize()));
        assertEquals(1.0e-10, set.getTolerance(), 1.0e-10);
    }

    @Test
    void testSameOrientationAs() {
        OrientedPoint pt1 = new OrientedPoint(new Vector1D(1.0), true, 1.0e-10);
        OrientedPoint pt2 = new OrientedPoint(new Vector1D(2.0), true, 1.0e-10);
        OrientedPoint pt3 = new OrientedPoint(new Vector1D(3.0), false, 1.0e-10);
        assertTrue(pt1.sameOrientationAs(pt2));
        assertFalse(pt1.sameOrientationAs(pt3));
    }

    @Test
    void testProject() {
        Vector1D location = new Vector1D(1.0);
        OrientedPoint pt = new OrientedPoint(location, true, 1.0e-10);
        assertEquals(1.0, pt.project(new Vector1D(5.0)).getX(), 1.0e-10);
    }

    @Test
    void testRevertSelf() {
        OrientedPoint pt = new OrientedPoint(new Vector1D(1.0), true, 1.0e-10);
        assertTrue(pt.isDirect());
        pt.revertSelf();
        assertFalse(pt.isDirect());
        pt.revertSelf();
        assertTrue(pt.isDirect());
    }
}
