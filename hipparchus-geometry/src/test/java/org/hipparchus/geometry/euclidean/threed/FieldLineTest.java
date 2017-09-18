/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
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
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Decimal64;
import org.hipparchus.util.Decimal64Field;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class FieldLineTest {

    @Test
    public void testContains() throws MathIllegalArgumentException, MathRuntimeException {
        FieldVector3D<Decimal64> p1 = FieldVector3D.getPlusK(Decimal64Field.getInstance());
        FieldLine<Decimal64> l = new FieldLine<>(p1,
                                                 new FieldVector3D<>(new Decimal64(0),
                                                                     new Decimal64(0),
                                                                     new Decimal64(2)),
                                                 1.0e-10);
        Assert.assertTrue(l.contains(p1));
        Assert.assertTrue(l.contains(new FieldVector3D<>(1.0, p1, 0.3, l.getDirection())));
        Assert.assertTrue(l.contains(new Vector3D(1.0, p1.toVector3D(), 0.3, l.getDirection().toVector3D())));
        FieldVector3D<Decimal64> u = l.getDirection().orthogonal();
        FieldVector3D<Decimal64> v = FieldVector3D.crossProduct(l.getDirection(), u);
        for (double a = 0; a < 2 * FastMath.PI; a += 0.3) {
            Decimal64 alpha = new Decimal64(a);
            Assert.assertTrue(! l.contains(p1.add(new FieldVector3D<>(alpha.cos(), u, alpha.sin(), v))));
        }
    }

    @Test
    public void testSimilar() throws MathIllegalArgumentException, MathRuntimeException {
        FieldVector3D<Decimal64> p1  = createVector(1.2, 3.4, -5.8);
        FieldVector3D<Decimal64> p2  = createVector(3.4, -5.8, 1.2);
        FieldLine<Decimal64>     lA  = new FieldLine<>(p1, p2, 1.0e-10);
        FieldLine<Decimal64>     lB  = new FieldLine<>(p2, p1, 1.0e-10);
        Assert.assertTrue(lA.isSimilarTo(lB));
        Assert.assertTrue(! lA.isSimilarTo(new FieldLine<>(p1, p1.add(lA.getDirection().orthogonal()), 1.0e-10)));
    }

    @Test
    public void testPointDistance() throws MathIllegalArgumentException {
        FieldLine<Decimal64> l = new FieldLine<>(createVector(0, 1, 1),
                                                 createVector(0, 2, 2),
                                                 1.0e-10);
        Assert.assertEquals(FastMath.sqrt(3.0 / 2.0),
                            l.distance(createVector(1, 0, 1)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0,
                            l.distance(createVector(0, -4, -4)).getReal(),
                            1.0e-10);
    }

    @Test
    public void testLineDistance() throws MathIllegalArgumentException {
        FieldLine<Decimal64> l = new FieldLine<>(createVector(0, 1, 1),
                                                 createVector(0, 2, 2),
                                                 1.0e-10);
        Assert.assertEquals(1.0,
                            l.distance(new FieldLine<>(createVector(1, 0, 1),
                                                       createVector(1, 0, 2),
                                                       1.0e-10)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.5,
                            l.distance(new FieldLine<>(createVector(-0.5, 0, 0),
                                                       createVector(-0.5, -1, -1),
                                                       1.0e-10)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.0,
                            l.distance(l).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.0,
                            l.distance(new FieldLine<>(createVector(0, -4, -4),
                                                       createVector(0, -5, -5),
                                                       1.0e-10)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.0,
                            l.distance(new FieldLine<>(createVector(0, -4, -4),
                                                       createVector(0, -3, -4),
                                                       1.0e-10)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.0,
                            l.distance(new FieldLine<>(createVector(0, -4, -4),
                                                       createVector(1, -4, -4),
                                                       1.0e-10)).getReal(),
                            1.0e-10);
        Assert.assertEquals(FastMath.sqrt(8),
                            l.distance(new FieldLine<>(createVector(0, -4, 0),
                                                       createVector(1, -4, 0),
                                                       1.0e-10)).getReal(),
                            1.0e-10);
    }

    @Test
    public void testClosest() throws MathIllegalArgumentException {
        FieldLine<Decimal64> l = new FieldLine<>(createVector(0, 1, 1),
                                                 createVector(0, 2, 2),
                                                 1.0e-10);
        Assert.assertEquals(0.0,
                            l.closestPoint(new FieldLine<>(createVector(1, 0, 1),
                                                           createVector(1, 0, 2),
                                                           1.0e-10)).distance(createVector(0, 0, 0)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.5,
                            l.closestPoint(new FieldLine<>(createVector(-0.5, 0, 0),
                                                           createVector(-0.5, -1, -1),
                                                           1.0e-10)).distance(createVector(-0.5, 0, 0)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.0,
                            l.closestPoint(l).distance(createVector(0, 0, 0)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.0,
                            l.closestPoint(new FieldLine<>(createVector(0, -4, -4),
                                                           createVector(0, -5, -5),
                                                           1.0e-10)).distance(createVector(0, 0, 0)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.0,
                            l.closestPoint(new FieldLine<>(createVector(0, -4, -4),
                                                           createVector(0, -3, -4),
                                                           1.0e-10)).distance(createVector(0, -4, -4)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.0,
                            l.closestPoint(new FieldLine<>(createVector(0, -4, -4),
                                                           createVector(1, -4, -4),
                                                           1.0e-10)).distance(createVector(0, -4, -4)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.0,
                            l.closestPoint(new FieldLine<>(createVector(0, -4, 0),
                                                           createVector(1, -4, 0),
                                                           1.0e-10)).distance(createVector(0, -2, -2)).getReal(),
                            1.0e-10);
    }

    @Test
    public void testIntersection() throws MathIllegalArgumentException {
        FieldLine<Decimal64> l = new FieldLine<>(createVector(0, 1, 1), createVector(0, 2, 2), 1.0e-10);
        Assert.assertNull(l.intersection(new FieldLine<>(createVector(1, 0, 1), createVector(1, 0, 2), 1.0e-10)));
        Assert.assertNull(l.intersection(new FieldLine<>(createVector(-0.5, 0, 0), createVector(-0.5, -1, -1), 1.0e-10)));
        Assert.assertEquals(0.0,
                            l.intersection(l).distance(createVector(0, 0, 0)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.0,
                            l.intersection(new FieldLine<>(createVector(0, -4, -4), createVector(0, -5, -5), 1.0e-10)).distance(createVector(0, 0, 0)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.0,
                            l.intersection(new FieldLine<>(createVector(0, -4, -4), createVector(0, -3, -4), 1.0e-10)).distance(createVector(0, -4, -4)).getReal(),
                            1.0e-10);
        Assert.assertEquals(0.0,
                            l.intersection(new FieldLine<>(createVector(0, -4, -4), createVector(1, -4, -4), 1.0e-10)).distance(createVector(0, -4, -4)).getReal(),
                            1.0e-10);
        Assert.assertNull(l.intersection(new FieldLine<>(createVector(0, -4, 0), createVector(1, -4, 0), 1.0e-10)));
    }

    @Test
    public void testRevert() {

        // setup
        FieldLine<Decimal64> line = new FieldLine<>(createVector(1653345.6696423641, 6170370.041579291, 90000),
                                                    createVector(1650757.5050732433, 6160710.879908984, 0.9),
                                                    1.0e-10);
        FieldVector3D<Decimal64> expected = line.getDirection().negate();

        // action
        FieldLine<Decimal64> reverted = line.revert();

        // verify
        Decimal64[] e = expected.toArray();
        Decimal64[] r = reverted.getDirection().toArray();
        Assert.assertEquals(e.length, e.length);
        for (int i = 0; i < e.length; ++i) {
            Assert.assertEquals(e[i].getReal(), r[i].getReal(), 1.0e-10);
        }

    }

    private FieldVector3D<Decimal64> createVector(double x, double y, double z) {
        return new FieldVector3D<>(new Decimal64(x), new Decimal64(y), new Decimal64(z));
    }

}
