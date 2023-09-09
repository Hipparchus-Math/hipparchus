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
package org.hipparchus.util;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.junit.Assert;
import org.junit.Test;

public class MerweUnscentedTransformTest {

    /** test state dimension equal to 0 */
    @Test(expected = MathIllegalArgumentException.class)
    public void testWrongStateDimension() {
        new MerweUnscentedTransform(0);
    }

    /** test weight computation */
    @Test
    public void testWeights() {

        // Initialize
        final int stateDim = 2;
        final MerweUnscentedTransform merwe = new MerweUnscentedTransform(stateDim);
        final RealVector wc = merwe.getWc();
        final RealVector wm = merwe.getWm();

        // Verify
        Assert.assertEquals(5,     wc.getDimension());
        Assert.assertEquals(5,     wm.getDimension());
        Assert.assertEquals(-0.25, wc.getEntry(0), Double.MIN_VALUE);
        Assert.assertEquals(-3.0,  wm.getEntry(0), Double.MIN_VALUE);
        Assert.assertEquals(1.0,   wc.getEntry(1), Double.MIN_VALUE);
        Assert.assertEquals(1.0,   wm.getEntry(1), Double.MIN_VALUE);
        Assert.assertEquals(1.0,   wc.getEntry(2), Double.MIN_VALUE);
        Assert.assertEquals(1.0,   wm.getEntry(2), Double.MIN_VALUE);

    }

    /** test unscented transform */
    @Test
    public void testUnscentedTransform() {

        // Initialize
        final int stateDim = 2;
        final MerweUnscentedTransform julier = new MerweUnscentedTransform(stateDim, 0.5, 2.0, 0.0);
        final RealVector state = MatrixUtils.createRealVector(new double[] {1.0, 1.0});
        final RealMatrix covariance = MatrixUtils.createRealDiagonalMatrix(new double[] {0.5, 0.5});

        // Action
        final RealVector[] sigma = julier.unscentedTransform(state, covariance);

        // Verify
        Assert.assertEquals(5, sigma.length);
        checkSigmaPoint(sigma[0], 1.0, 1.0);
        checkSigmaPoint(sigma[1], 1.5, 1.0);
        checkSigmaPoint(sigma[2], 1.0, 1.5);
        checkSigmaPoint(sigma[3], 0.5, 1.0);
        checkSigmaPoint(sigma[4], 1.0, 0.5);

    }

    private static void checkSigmaPoint(final RealVector sigma, final double ref1, final double ref2) {
        Assert.assertEquals(ref1,  sigma.getEntry(0), Double.MIN_VALUE);
        Assert.assertEquals(ref2,  sigma.getEntry(1), Double.MIN_VALUE);
    }

}
