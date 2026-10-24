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
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class BFGSUpdaterTest {

    /** Tolerance for order-1 entries after one rank-two update. */
    private static final double TOL = 1.0e-14;

    @Test
    public void testRegularUpdateSatisfiesSecantCondition() {
        final BFGSUpdater bfgs = createUpdater();
        final RealVector s = new ArrayRealVector(new double[] { 1.0, 0.0 });
        final RealVector y = new ArrayRealVector(new double[] { 2.0, 0.5 });
        Assertions.assertEquals(0, bfgs.update(s, y));

        // H+ = I - s sᵀ / sᵀs + y yᵀ / sᵀy
        final RealMatrix expected = MatrixUtils.createRealMatrix(new double[][] {
            { 2.0, 0.5 },
            { 0.5, 1.125 }
        });
        checkMatrix(expected, bfgs.getHessian(), TOL);
        // secant condition H+ s = y
        Assertions.assertEquals(0.0, bfgs.getHessian().operate(s).subtract(y).getNorm(), TOL);
    }

    @Test
    public void testResetHessianFactor() {
        final BFGSUpdater bfgs = createUpdater();
        Assertions.assertEquals(0,
                                bfgs.update(new ArrayRealVector(new double[] { 1.0, 0.0 }),
                                            new ArrayRealVector(new double[] { 2.0, 0.5 })));
        // the scale factor is yᵀy / sᵀy = 4.25 / 2
        bfgs.resetHessianFactor();
        checkMatrix(MatrixUtils.createRealIdentityMatrix(2).scalarMultiply(2.125), bfgs.getHessian(), TOL);
    }

    @Test
    public void testZeroStepResetsHessian() {
        final BFGSUpdater bfgs = createUpdater();
        bfgs.resetHessian(4.0);
        Assertions.assertEquals(1, bfgs.update(new ArrayRealVector(2), new ArrayRealVector(new double[] { 1.0, 1.0 })));
        checkMatrix(MatrixUtils.createRealIdentityMatrix(2), bfgs.getHessian(), 0.0);
    }

    @Test
    public void testNaNCurvatureSkipsUpdate() {
        final BFGSUpdater bfgs = createUpdater();
        Assertions.assertEquals(2,
                                bfgs.update(new ArrayRealVector(new double[] { 1.0, 0.0 }),
                                            new ArrayRealVector(new double[] { Double.NaN, 0.0 })));
        checkMatrix(MatrixUtils.createRealIdentityMatrix(2), bfgs.getHessian(), 0.0);
    }

    @Test
    public void testTinyStepKeepsUnitScaleFactor() {
        // the curvature y/s = 1e10 is large, but the step is so small that sᵀy = 1e-310
        // is below Precision.SAFE_MIN, so the scale factor yᵀy / sᵀy is not used
        final BFGSUpdater bfgs = createUpdater();
        Assertions.assertEquals(0,
                                bfgs.update(new ArrayRealVector(new double[] { 1.0e-160, 0.0 }),
                                            new ArrayRealVector(new double[] { 1.0e-150, 0.0 })));
        bfgs.resetHessianFactor();
        checkMatrix(MatrixUtils.createRealIdentityMatrix(2), bfgs.getHessian(), 0.0);
    }

    private static BFGSUpdater createUpdater() {
        return new BFGSUpdater(MatrixUtils.createRealIdentityMatrix(2), 1.0e-8, true, 1.0e-12);
    }

    private static void checkMatrix(final RealMatrix expected, final RealMatrix actual, final double tol) {
        Assertions.assertEquals(0.0, actual.subtract(expected).getNorm1(), tol);
    }

}
