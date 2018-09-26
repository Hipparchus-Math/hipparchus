/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.analysis.interpolation;

import java.util.Arrays;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.junit.Assert;
import org.junit.Test;

public class GridAxisTest {

    @Test
    public void testLinearAscending() {
        for (int n = 2; n < 12; ++n) {
            checkAscending(createLinear(25, n));
        }
    }

    @Test
    public void testLinearDescending() {
        for (int n = 2; n < 12; ++n) {
            checkDescending(createLinear(25, n));
        }
    }

    @Test
    public void testLinearRandomAccess() {
        final RandomGenerator random = new Well1024a(0x967ab81207d7b2f4l);
        for (int n = 2; n < 12; ++n) {
            checkRandomAccess(random, createLinear(25, n));
        }
    }

    @Test
    public void testQuadraticAscending() {
        for (int n = 2; n < 12; ++n) {
            checkAscending(createQuadratic(25, n));
        }
    }

    @Test
    public void testQuadraticDescending() {
        for (int n = 2; n < 12; ++n) {
            checkDescending(createQuadratic(25, n));
        }
    }

    @Test
    public void testQuadraticRandomAccess() {
        final RandomGenerator random = new Well1024a(0x80fc3b30a2da4549l);
        for (int n = 2; n < 12; ++n) {
            checkRandomAccess(random, createQuadratic(25, n));
        }
    }

    @Test
    public void testIrregularAscending() {
        final RandomGenerator random = new Well1024a(0x66133fa03616fbc9l);
        for (int n = 2; n < 12; ++n) {
            checkAscending(createIrregular(random, 25, n));
        }
    }

    @Test
    public void testIrregularDescending() {
        final RandomGenerator random = new Well1024a(0x72404bbdc66bc2c7l);
        for (int n = 2; n < 12; ++n) {
            checkDescending(createIrregular(random, 25, n));
        }
    }

    @Test
    public void testIrregularRandomAccess() {
        final RandomGenerator random = new Well1024a(0x254fbf2a8207e0f6l);
        for (int n = 2; n < 12; ++n) {
            checkRandomAccess(random, createIrregular(random, 25, n));
        }
    }

    @Test
    public void testTooSmallGrid() {
        try {
            new GridAxis(new double[3], 4);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.INSUFFICIENT_DIMENSION, miae.getSpecifier());
            Assert.assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            Assert.assertEquals(4, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    public void testDuplicate() {
        try {
            new GridAxis(new double[] { 0.0, 1.0, 2.0, 2.0, 3.0 }, 2);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.NOT_STRICTLY_INCREASING_SEQUENCE, miae.getSpecifier());
            Assert.assertEquals(2.0, ((Double)  miae.getParts()[0]).doubleValue(), 1.0e-15);
            Assert.assertEquals(2.0, ((Double)  miae.getParts()[1]).doubleValue(), 1.0e-15);
            Assert.assertEquals(3,   ((Integer) miae.getParts()[2]).intValue());
            Assert.assertEquals(2,   ((Integer) miae.getParts()[3]).intValue());
        }
    }

    @Test
    public void testUnsorted() {
        try {
            new GridAxis(new double[] { 0.0, 1.0, 0.5, 2.0, 3.0 }, 2);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.NOT_STRICTLY_INCREASING_SEQUENCE, miae.getSpecifier());
            Assert.assertEquals(0.5, ((Double)  miae.getParts()[0]).doubleValue(), 1.0e-15);
            Assert.assertEquals(1.0, ((Double)  miae.getParts()[1]).doubleValue(), 1.0e-15);
            Assert.assertEquals(2,   ((Integer) miae.getParts()[2]).intValue());
            Assert.assertEquals(1,   ((Integer) miae.getParts()[3]).intValue());
        }
    }

    private GridAxis createLinear(final int size, final int n) {
        final double[] gridData = new double[size];
        for (int i = 0; i < size; ++i) {
            gridData[i] = 2 * i + 0.5;
        }
        return create(gridData, n);
    }

    private GridAxis createQuadratic(final int size, final int n) {
        final double[] gridData = new double[size];
        for (int i = 0; i < size; ++i) {
            gridData[i] = (i + 0.5) * (i + 3);
        }
        return create(gridData, n);
    }

    private GridAxis createIrregular(final RandomGenerator random, final int size, final int n) {
        final double[] gridData = new double[size];
        for (int i = 0; i < size; ++i) {
            gridData[i] = 50.0 * random.nextDouble();
        }
        Arrays.sort(gridData);
        return create(gridData, n);
    }

    private GridAxis create(final double[] gridData, final int n) {
        final GridAxis gridAxis = new GridAxis(gridData, n);
        Assert.assertEquals(n, gridAxis.getN());
        for (int i = 0; i < 5; ++i) {
            Assert.assertEquals(gridData[i], gridAxis.node(i), 1.0e-15);
        }
        return gridAxis;
    }

    private void checkAscending(final GridAxis gridAxis) {
        final double inf = gridAxis.node(0) - 2.0;
        final double sup = gridAxis.node(gridAxis.size() - 1) + 2.0;
        for (double t = inf; t < sup; t += 0.125) {
            checkInterpolation(t, gridAxis);
        }
    }

    private void checkDescending(final GridAxis gridAxis) {
        final double inf = gridAxis.node(0) - 2.0;
        final double sup = gridAxis.node(gridAxis.size() - 1) + 2.0;
        for (double t = sup; t > inf; t -= 0.125) {
            checkInterpolation(t, gridAxis);
        }
    }

    private void checkRandomAccess(final RandomGenerator random, final GridAxis gridAxis) {
        final double inf = gridAxis.node(0) - 2.0;
        final double sup = gridAxis.node(gridAxis.size() - 1) + 2.0;
        for (int i = 0; i < 1000; ++i) {
            checkInterpolation(inf + random.nextDouble() * (sup - inf), gridAxis);
        }
    }

    private void checkInterpolation(final double t, final GridAxis gridAxis) {
        final int s = gridAxis.size();
        final int n = gridAxis.getN();
        final int o = (n - 1) / 2;
        final int p = n / 2;
        final int i = gridAxis.interpolationIndex(t);
        Assert.assertTrue(i >= 0);
        Assert.assertTrue(i + n - 1 < s);
        if (t < gridAxis.node(0)) {
            // extrapolating below grid
            Assert.assertEquals(0, i);
        } else if (t < gridAxis.node(s - 1)) {

            // interpolating within the grid
            // the nodes should surround the test value
            Assert.assertTrue(gridAxis.node(i) <= t);
            Assert.assertTrue(gridAxis.node(i + n - 1) > t);

            if (t >= gridAxis.node(o) && t <  gridAxis.node(s - p)) {
                // interpolation in the part of the grid where balancing is possible
                // the central nodes should surround the test value
                Assert.assertTrue(gridAxis.node(i + o)     <= t);
                Assert.assertTrue(gridAxis.node(i + o + 1) >  t);
            }

        } else {
            // extrapolating above grid
            Assert.assertEquals(s - n, i);
        }
    }

}
