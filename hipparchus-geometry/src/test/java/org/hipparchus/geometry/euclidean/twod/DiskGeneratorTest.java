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
package org.hipparchus.geometry.euclidean.twod;

import org.hipparchus.geometry.enclosing.EnclosingBall;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.UnitSphereRandomVectorGenerator;
import org.hipparchus.random.Well1024a;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;


class DiskGeneratorTest {

    @Test
    void testSupport0Point() {
        List<Vector2D> support = Arrays.asList(new Vector2D[0]);
        EnclosingBall<Euclidean2D, Vector2D> disk = new DiskGenerator().ballOnSupport(support);
        assertTrue(disk.getRadius() < 0);
        assertEquals(0, disk.getSupportSize());
        assertEquals(0, disk.getSupport().length);
    }

    @Test
    void testSupport1Point() {
        List<Vector2D> support = Arrays.asList(new Vector2D(1, 2));
        EnclosingBall<Euclidean2D, Vector2D> disk = new DiskGenerator().ballOnSupport(support);
        assertEquals(0.0, disk.getRadius(), 1.0e-10);
        assertTrue(disk.contains(support.get(0)));
        assertTrue(disk.contains(support.get(0), 0.5));
        assertFalse(disk.contains(new Vector2D(support.get(0).getX() + 0.1,
                                                      support.get(0).getY() - 0.1),
                                         0.001));
        assertTrue(disk.contains(new Vector2D(support.get(0).getX() + 0.1,
                                                     support.get(0).getY() - 0.1),
                                        0.5));
        assertEquals(0, support.get(0).distance(disk.getCenter()), 1.0e-10);
        assertEquals(1, disk.getSupportSize());
        assertSame(support.get(0), disk.getSupport()[0]);
    }

    @Test
    void testSupport2Points() {
        List<Vector2D> support = Arrays.asList(new Vector2D(1, 0),
                                               new Vector2D(3, 0));
        EnclosingBall<Euclidean2D, Vector2D> disk = new DiskGenerator().ballOnSupport(support);
        assertEquals(1.0, disk.getRadius(), 1.0e-10);
        int i = 0;
        for (Vector2D v : support) {
            assertTrue(disk.contains(v));
            assertEquals(1.0, v.distance(disk.getCenter()), 1.0e-10);
            assertSame(v, disk.getSupport()[i++]);
        }
        assertTrue(disk.contains(new Vector2D(2, 0.9)));
        assertFalse(disk.contains(Vector2D.ZERO));
        assertEquals(0.0, new Vector2D(2, 0).distance(disk.getCenter()), 1.0e-10);
        assertEquals(2, disk.getSupportSize());
    }

    @Test
    void testSupport3Points() {
        List<Vector2D> support = Arrays.asList(new Vector2D(1, 0),
                                               new Vector2D(3, 0),
                                               new Vector2D(2, 2));
        EnclosingBall<Euclidean2D, Vector2D> disk = new DiskGenerator().ballOnSupport(support);
        assertEquals(5.0 / 4.0, disk.getRadius(), 1.0e-10);
        int i = 0;
        for (Vector2D v : support) {
            assertTrue(disk.contains(v));
            assertEquals(5.0 / 4.0, v.distance(disk.getCenter()), 1.0e-10);
            assertSame(v, disk.getSupport()[i++]);
        }
        assertTrue(disk.contains(new Vector2D(2, 0.9)));
        assertFalse(disk.contains(new Vector2D(0.9,  0)));
        assertFalse(disk.contains(new Vector2D(3.1,  0)));
        assertTrue(disk.contains(new Vector2D(2.0, -0.499)));
        assertFalse(disk.contains(new Vector2D(2.0, -0.501)));
        assertEquals(0.0, new Vector2D(2.0, 3.0 / 4.0).distance(disk.getCenter()), 1.0e-10);
        assertEquals(3, disk.getSupportSize());
    }

    @Test
    void testRandom() {
        final RandomGenerator random = new Well1024a(0x12faa818373ffe90l);
        final UnitSphereRandomVectorGenerator sr = new UnitSphereRandomVectorGenerator(2, random);
        for (int i = 0; i < 500; ++i) {
            double d = 25 * random.nextDouble();
            double refRadius = 10 * random.nextDouble();
            Vector2D refCenter = new Vector2D(d, new Vector2D(sr.nextVector()));
            List<Vector2D> support = new ArrayList<Vector2D>();
            for (int j = 0; j < 3; ++j) {
                support.add(new Vector2D(1.0, refCenter, refRadius, new Vector2D(sr.nextVector())));
            }
            EnclosingBall<Euclidean2D, Vector2D> disk = new DiskGenerator().ballOnSupport(support);
            assertEquals(0.0, refCenter.distance(disk.getCenter()), 3e-9 * refRadius);
            assertEquals(refRadius, disk.getRadius(), 7e-10 * refRadius);
        }

    }
}
