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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.geometry.partitioning.Region;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Precision;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class IntervalTest {

    @Test
    void testInterval() {
        Interval interval = new Interval(2.3, 5.7);
        assertEquals(3.4, interval.getSize(), 1.0e-10);
        assertEquals(4.0, interval.getBarycenter(), 1.0e-10);
        assertEquals(Region.Location.BOUNDARY, interval.checkPoint(2.3, 1.0e-10));
        assertEquals(Region.Location.BOUNDARY, interval.checkPoint(5.7, 1.0e-10));
        assertEquals(Region.Location.OUTSIDE,  interval.checkPoint(1.2, 1.0e-10));
        assertEquals(Region.Location.OUTSIDE,  interval.checkPoint(8.7, 1.0e-10));
        assertEquals(Region.Location.INSIDE,   interval.checkPoint(3.0, 1.0e-10));
        assertEquals(2.3, interval.getInf(), 1.0e-10);
        assertEquals(5.7, interval.getSup(), 1.0e-10);
    }

    @Test
    void testTolerance() {
        Interval interval = new Interval(2.3, 5.7);
        assertEquals(Region.Location.OUTSIDE,  interval.checkPoint(1.2, 1.0));
        assertEquals(Region.Location.BOUNDARY, interval.checkPoint(1.2, 1.2));
        assertEquals(Region.Location.OUTSIDE,  interval.checkPoint(8.7, 2.9));
        assertEquals(Region.Location.BOUNDARY, interval.checkPoint(8.7, 3.1));
        assertEquals(Region.Location.INSIDE,   interval.checkPoint(3.0, 0.6));
        assertEquals(Region.Location.BOUNDARY, interval.checkPoint(3.0, 0.8));
    }

    @Test
    void testInfinite() {
        Interval interval = new Interval(9.0, Double.POSITIVE_INFINITY);
        assertEquals(Region.Location.BOUNDARY, interval.checkPoint(9.0, 1.0e-10));
        assertEquals(Region.Location.OUTSIDE,  interval.checkPoint(8.4, 1.0e-10));
        for (double e = 1.0; e <= 6.0; e += 1.0) {
            assertEquals(Region.Location.INSIDE,
                                interval.checkPoint(FastMath.pow(10.0, e), 1.0e-10));
        }
        assertTrue(Double.isInfinite(interval.getSize()));
        assertEquals(9.0, interval.getInf(), 1.0e-10);
        assertTrue(Double.isInfinite(interval.getSup()));

    }

    @Test
    void testWholeLine() {
        Interval interval = new Interval(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        assertTrue(Double.isInfinite(interval.getSize()));
        assertEquals(Region.Location.INSIDE, interval.checkPoint(-Double.MAX_VALUE, 1.0e-10));
        assertEquals(Region.Location.INSIDE, interval.checkPoint(+Double.MAX_VALUE, 1.0e-10));
    }

    @Test
    void testSinglePoint() {
        Interval interval = new Interval(1.0, 1.0);
        assertEquals(0.0, interval.getSize(), Precision.SAFE_MIN);
        assertEquals(1.0, interval.getBarycenter(), Precision.EPSILON);
    }

    // MATH-1256
    @Test
    void testStrictOrdering() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new Interval(0, -1);
        });
    }
}
