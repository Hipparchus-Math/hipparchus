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
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class S1PointTest {

    @Test
    public void testS1Point() {
        for (int k = -2; k < 3; ++k) {
            S1Point p = new S1Point(1.0 + k * MathUtils.TWO_PI);
            Assertions.assertEquals(FastMath.cos(1.0), p.getVector().getX(), 1.0e-10);
            Assertions.assertEquals(FastMath.sin(1.0), p.getVector().getY(), 1.0e-10);
            Assertions.assertFalse(p.isNaN());
        }
    }

    @Test
    public void testNaN() {
        Assertions.assertTrue(S1Point.NaN.isNaN());
        Assertions.assertEquals(S1Point.NaN, new S1Point(Double.NaN));
        Assertions.assertNotEquals(S1Point.NaN, new S1Point(1.0));
    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    public void testEquals() {
        S1Point a = new S1Point(1.0);
        S1Point b = new S1Point(1.0);
        Assertions.assertEquals(a.hashCode(), b.hashCode());
        Assertions.assertFalse(a == b);
        Assertions.assertEquals(a, b);
        Assertions.assertEquals(a, a);
        Assertions.assertNotEquals('a', a);
        Assertions.assertEquals(S1Point.NaN, S1Point.NaN);
        Assertions.assertEquals(S1Point.NaN, new S1Point(Double.NaN));
    }

    @Test
    public void testEqualsIeee754() {
        S1Point a = new S1Point(1.0);
        S1Point b = new S1Point(1.0);
        Assertions.assertEquals(a.hashCode(), b.hashCode());
        Assertions.assertFalse(a == b);
        Assertions.assertTrue(a.equalsIeee754(b));
        Assertions.assertTrue(a.equalsIeee754(a));
        Assertions.assertFalse(a.equalsIeee754('a'));
        Assertions.assertFalse(S1Point.NaN.equalsIeee754(S1Point.NaN));
        Assertions.assertFalse(S1Point.NaN.equalsIeee754(new S1Point(Double.NaN)));
    }

    @Test
    public void testDistance() {
        S1Point a = new S1Point(1.0);
        S1Point b = new S1Point(a.getAlpha() + 0.5 * FastMath.PI);
        Assertions.assertEquals(0.5 * FastMath.PI, a.distance(b), 1.0e-10);
    }

    @Test
    public void testSpace() {
        S1Point a = new S1Point(1.0);
        Assertions.assertTrue(a.getSpace() instanceof Sphere1D);
        Assertions.assertEquals(1, a.getSpace().getDimension());
        try {
            a.getSpace().getSubSpace();
            Assertions.fail("an exception should have been thrown");
        } catch (MathRuntimeException muoe) {
            // expected
        }
    }

}
