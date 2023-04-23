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

import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.Assert;
import org.junit.Test;

public class S1PointTest {

    @Test
    public void testS1Point() {
        for (int k = -2; k < 3; ++k) {
            S1Point p = new S1Point(1.0 + k * MathUtils.TWO_PI);
            Assert.assertEquals(FastMath.cos(1.0), p.getVector().getX(), 1.0e-10);
            Assert.assertEquals(FastMath.sin(1.0), p.getVector().getY(), 1.0e-10);
            Assert.assertFalse(p.isNaN());
        }
    }

    @Test
    public void testNaN() {
        Assert.assertTrue(S1Point.NaN.isNaN());
        Assert.assertTrue(S1Point.NaN.equals(new S1Point(Double.NaN)));
        Assert.assertFalse(new S1Point(1.0).equals(S1Point.NaN));
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    public void testEquals() {
        S1Point a = new S1Point(1.0);
        S1Point b = new S1Point(1.0);
        Assert.assertEquals(a.hashCode(), b.hashCode());
        Assert.assertFalse(a == b);
        Assert.assertTrue(a.equals(b));
        Assert.assertTrue(a.equals(a));
        Assert.assertFalse(a.equals('a'));
        Assert.assertTrue(S1Point.NaN.equals(S1Point.NaN));
        Assert.assertTrue(S1Point.NaN.equals(new S1Point(Double.NaN)));
    }

    @Test
    public void testEqualsIeee754() {
        S1Point a = new S1Point(1.0);
        S1Point b = new S1Point(1.0);
        Assert.assertEquals(a.hashCode(), b.hashCode());
        Assert.assertFalse(a == b);
        Assert.assertTrue(a.equalsIeee754(b));
        Assert.assertTrue(a.equalsIeee754(a));
        Assert.assertFalse(a.equalsIeee754('a'));
        Assert.assertFalse(S1Point.NaN.equalsIeee754(S1Point.NaN));
        Assert.assertFalse(S1Point.NaN.equalsIeee754(new S1Point(Double.NaN)));
    }

    @Test
    public void testDistance() {
        S1Point a = new S1Point(1.0);
        S1Point b = new S1Point(a.getAlpha() + 0.5 * FastMath.PI);
        Assert.assertEquals(0.5 * FastMath.PI, a.distance(b), 1.0e-10);
    }

    @Test
    public void testSpace() {
        S1Point a = new S1Point(1.0);
        Assert.assertTrue(a.getSpace() instanceof Sphere1D);
        Assert.assertEquals(1, a.getSpace().getDimension());
        try {
            a.getSpace().getSubSpace();
            Assert.fail("an exception should have been thrown");
        } catch (MathRuntimeException muoe) {
            // expected
        }
    }

}
