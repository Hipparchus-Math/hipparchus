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
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

class JulierUnscentedTransformTest {

    /** test state dimension equal to 0 */
    @Test
    void testWrongStateDimension() {
        assertThrows(MathIllegalArgumentException.class, () -> new JulierUnscentedTransform(0));
    }

    /** test weight computation */
    @Test
    void testWeights() {

        // Initialize
        final int stateDim = 2;
        final JulierUnscentedTransform julier = new JulierUnscentedTransform(stateDim, 0.0);
        final RealVector wc = julier.getWc();
        final RealVector wm = julier.getWm();

        // Verify
        assertEquals(5,    wc.getDimension());
        assertEquals(5,    wm.getDimension());
        assertEquals(0.0,  wc.getEntry(0), Double.MIN_VALUE);
        assertEquals(0.0,  wm.getEntry(0), Double.MIN_VALUE);
        assertEquals(0.25, wc.getEntry(1), Double.MIN_VALUE);
        assertEquals(0.25, wm.getEntry(1), Double.MIN_VALUE);
        assertEquals(0.25, wc.getEntry(2), Double.MIN_VALUE);
        assertEquals(0.25, wm.getEntry(2), Double.MIN_VALUE);

    }

    /** Test unscented transform */
    @Test
    void testUnscentedTransform() {

        // Initialize
        final int stateDim = 2;
        final JulierUnscentedTransform julier = new JulierUnscentedTransform(stateDim);
        final RealVector state = MatrixUtils.createRealVector(new double[] {1.0, 1.0});
        final RealMatrix covariance = MatrixUtils.createRealDiagonalMatrix(new double[] {0.5, 0.5});

        // Action
        final RealVector[] sigma = julier.unscentedTransform(state, covariance);

        // Verify
        assertEquals(5, sigma.length);
        checkSigmaPoint(sigma[0], 1.0, 1.0);
        checkSigmaPoint(sigma[1], 2.0, 1.0);
        checkSigmaPoint(sigma[2], 1.0, 2.0);
        checkSigmaPoint(sigma[3], 0.0, 1.0);
        checkSigmaPoint(sigma[4], 1.0, 0.0);

    }

    /** Test inverse unscented transform */
    @Test
    void testInverseUnscentedTransform() {
        
        // Initialize
        final int stateDim = 2;
        final JulierUnscentedTransform julier = new JulierUnscentedTransform(stateDim);
        final RealVector[] sigmaPoints = new RealVector[] {MatrixUtils.createRealVector(new double[] {1.0, 1.0}),
                                                           MatrixUtils.createRealVector(new double[] {2.0, 1.0}),
                                                           MatrixUtils.createRealVector(new double[] {1.0, 2.0}),
                                                           MatrixUtils.createRealVector(new double[] {0.0, 1.0}),
                                                           MatrixUtils.createRealVector(new double[] {1.0, 0.0})};
        // Action
        final Pair<RealVector, RealMatrix> inverse = julier.inverseUnscentedTransform(sigmaPoints);
        final RealVector state = inverse.getFirst();
        final RealMatrix covariance = inverse.getSecond();
        
        // Verify
        assertEquals(2, state.getDimension());
        assertEquals(1.0, state.getEntry(0), 0.);
        assertEquals(1.0, state.getEntry(1), 0.);
        
        assertEquals(2, covariance.getColumnDimension());
        assertEquals(2, covariance.getRowDimension());
        assertEquals(0.5, covariance.getEntry(0, 0), 0.);
        assertEquals(0.0, covariance.getEntry(0, 1), 0.);
        assertEquals(0.0, covariance.getEntry(1, 0), 0.);
        assertEquals(0.5, covariance.getEntry(1, 1), 0.);
    }

    private static void checkSigmaPoint(final RealVector sigma, final double ref1, final double ref2) {
        assertEquals(ref1,  sigma.getEntry(0), Double.MIN_VALUE);
        assertEquals(ref2,  sigma.getEntry(1), Double.MIN_VALUE);
    }

}
