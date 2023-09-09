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
import org.hipparchus.geometry.partitioning.Region;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;
import org.junit.Assert;
import org.junit.Test;

public class ArcTest {

    @Test
    public void testArc() {
        Arc arc = new Arc(2.3, 5.7, 1.0e-10);
        Assert.assertEquals(3.4, arc.getSize(), 1.0e-10);
        Assert.assertEquals(4.0, arc.getBarycenter(), 1.0e-10);
        Assert.assertEquals(Region.Location.BOUNDARY, arc.checkPoint(2.3));
        Assert.assertEquals(Region.Location.BOUNDARY, arc.checkPoint(5.7));
        Assert.assertEquals(Region.Location.OUTSIDE,  arc.checkPoint(1.2));
        Assert.assertEquals(Region.Location.OUTSIDE,  arc.checkPoint(8.5));
        Assert.assertEquals(Region.Location.INSIDE,   arc.checkPoint(8.7));
        Assert.assertEquals(Region.Location.INSIDE,   arc.checkPoint(3.0));
        Assert.assertEquals(2.3, arc.getInf(), 1.0e-10);
        Assert.assertEquals(5.7, arc.getSup(), 1.0e-10);
        Assert.assertEquals(4.0, arc.getBarycenter(), 1.0e-10);
        Assert.assertEquals(3.4, arc.getSize(), 1.0e-10);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testWrongInterval() {
        new Arc(1.2, 0.0, 1.0e-10);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testTooSmallTolerance() {
        new Arc(0.0, 1.0, 0.9 * Sphere1D.SMALLEST_TOLERANCE);
    }

    @Test
    public void testTolerance() {
        Assert.assertEquals(Region.Location.OUTSIDE,  new Arc(2.3, 5.7, 1.0).checkPoint(1.2));
        Assert.assertEquals(Region.Location.BOUNDARY, new Arc(2.3, 5.7, 1.2).checkPoint(1.2));
        Assert.assertEquals(Region.Location.OUTSIDE,  new Arc(2.3, 5.7, 0.7).checkPoint(6.5));
        Assert.assertEquals(Region.Location.BOUNDARY, new Arc(2.3, 5.7, 0.9).checkPoint(6.5));
        Assert.assertEquals(Region.Location.INSIDE,   new Arc(2.3, 5.7, 0.6).checkPoint(3.0));
        Assert.assertEquals(Region.Location.BOUNDARY, new Arc(2.3, 5.7, 0.8).checkPoint(3.0));
    }

    @Test
    public void testFullCircle() {
        Arc arc = new Arc(9.0, 9.0, 1.0e-10);
        // no boundaries on a full circle
        Assert.assertEquals(Region.Location.INSIDE, arc.checkPoint(9.0));
        Assert.assertEquals(.0, arc.getInf(), 1.0e-10);
        Assert.assertEquals(MathUtils.TWO_PI, arc.getSup(), 1.0e-10);
        Assert.assertEquals(2.0 * FastMath.PI, arc.getSize(), 1.0e-10);
        for (double alpha = -20.0; alpha <= 20.0; alpha += 0.1) {
            Assert.assertEquals(Region.Location.INSIDE, arc.checkPoint(alpha));
        }
    }

    @Test
    public void testSmall() {
        Arc arc = new Arc(1.0, FastMath.nextAfter(1.0, Double.POSITIVE_INFINITY), 1.01 * Sphere1D.SMALLEST_TOLERANCE);
        Assert.assertEquals(2 * Precision.EPSILON, arc.getSize(), Precision.SAFE_MIN);
        Assert.assertEquals(1.0, arc.getBarycenter(), Precision.EPSILON);
    }

    /** Check {@link Arc#getOffset(double)}. */
    @Test
    public void testGetOffset() {
        // setup
        double twopi = 2 * FastMath.PI;
        Arc arc = new Arc(twopi + 1, twopi + 2, 1e-6);

        // action & verify
        double tol = FastMath.ulp(twopi);
        Assert.assertEquals(arc.getOffset(0), 1, tol);
        Assert.assertEquals(arc.getOffset(1), 0, tol);
        Assert.assertEquals(arc.getOffset(1.5), -0.5, tol);
        Assert.assertEquals(arc.getOffset(2), 0, tol);
        Assert.assertEquals(arc.getOffset(3), 1, tol);
    }

}
