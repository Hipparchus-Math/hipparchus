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
import org.hipparchus.geometry.spherical.oned.Sphere1D;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.Assert;
import org.junit.Test;

public class S2PointTest {


    @Test
    public void testS2Point() {
        for (int k = -2; k < 3; ++k) {
            S2Point p = new S2Point(1.0 + k * MathUtils.TWO_PI, 1.4);
            Assert.assertEquals(1.0 + k * MathUtils.TWO_PI, p.getTheta(), 1.0e-10);
            Assert.assertEquals(1.4, p.getPhi(), 1.0e-10);
            Assert.assertEquals(FastMath.cos(1.0) * FastMath.sin(1.4), p.getVector().getX(), 1.0e-10);
            Assert.assertEquals(FastMath.sin(1.0) * FastMath.sin(1.4), p.getVector().getY(), 1.0e-10);
            Assert.assertEquals(FastMath.cos(1.4), p.getVector().getZ(), 1.0e-10);
            Assert.assertFalse(p.isNaN());
        }
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testNegativePolarAngle() {
        new S2Point(1.0, -1.0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testTooLargePolarAngle() {
        new S2Point(1.0, 3.5);
    }

    @Test
    public void testNaN() {
        Assert.assertTrue(S2Point.NaN.isNaN());
        Assert.assertTrue(S2Point.NaN.equals(new S2Point(Double.NaN, 1.0)));
        Assert.assertFalse(new S2Point(1.0, 1.3).equals(S2Point.NaN));
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    public void testEquals() {
        S2Point a = new S2Point(1.0, 1.0);
        S2Point b = new S2Point(1.0, 1.0);
        Assert.assertEquals(a.hashCode(), b.hashCode());
        Assert.assertFalse(a == b);
        Assert.assertTrue(a.equals(b));
        Assert.assertTrue(a.equals(a));
        Assert.assertFalse(a.equals('a'));
        Assert.assertTrue(S2Point.NaN.equals(S2Point.NaN));
        Assert.assertTrue(S2Point.NaN.equals(new S2Point(Double.NaN, 0.0)));
        Assert.assertTrue(S2Point.NaN.equals(new S2Point(0.0, Double.NaN)));
    }

    @Test
    public void testEqualsIeee754() {
        S2Point a = new S2Point(1.0, 1.0);
        S2Point b = new S2Point(1.0, 1.0);
        Assert.assertEquals(a.hashCode(), b.hashCode());
        Assert.assertFalse(a == b);
        Assert.assertTrue(a.equalsIeee754(b));
        Assert.assertTrue(a.equalsIeee754(a));
        Assert.assertFalse(a.equalsIeee754('a'));
        Assert.assertFalse(S2Point.NaN.equalsIeee754(S2Point.NaN));
        Assert.assertFalse(S2Point.NaN.equalsIeee754(new S2Point(Double.NaN, 0.0)));
        Assert.assertFalse(S2Point.NaN.equalsIeee754(new S2Point(0.0, Double.NaN)));
    }

    @Test
    public void testDistance() {
        S2Point a = new S2Point(1.0, 0.5 * FastMath.PI);
        S2Point b = new S2Point(a.getTheta() + 0.5 * FastMath.PI, a.getPhi());
        Assert.assertEquals(0.5 * FastMath.PI, a.distance(b), 1.0e-10);
        Assert.assertEquals(FastMath.PI, a.distance(a.negate()), 1.0e-10);
        Assert.assertEquals(0.5 * FastMath.PI, S2Point.MINUS_I.distance(S2Point.MINUS_K), 1.0e-10);
        Assert.assertEquals(0.0, new S2Point(1.0, 0).distance(new S2Point(2.0, 0)), 1.0e-10);
    }

    @Test
    public void testNegate() {
        RandomGenerator generator = new Well1024a(0x79d1bc2e0999d238l);
        for (int i = 0; i < 100000; ++i) {
            S2Point p = new S2Point(MathUtils.TWO_PI * generator.nextDouble(),
                                    FastMath.PI * generator.nextDouble());
            S2Point np = new S2Point(p.negate().getTheta(), p.negate().getPhi());
            Assert.assertEquals(FastMath.PI, p.distance(np), 1.4e-15);
        }
    }

    @Test
    public void testSpace() {
        S2Point a = new S2Point(1.0, 1.0);
        Assert.assertTrue(a.getSpace() instanceof Sphere2D);
        Assert.assertEquals(2, a.getSpace().getDimension());
        Assert.assertTrue(a.getSpace().getSubSpace() instanceof Sphere1D);
    }

}
