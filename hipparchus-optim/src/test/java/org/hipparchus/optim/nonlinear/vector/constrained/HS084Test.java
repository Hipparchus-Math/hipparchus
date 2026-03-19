/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
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

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.jupiter.api.Test;

/**
 * HS84 (TP84): 5 variables with 3 interval constraints represented as 6 inequalities.
 */
public class HS084Test {

    private static final int DIM      = 5;
    private static final int NUM_INEQ = 6;

    // a(1..21), aligned with original TP84 Fortran data (index 0 unused)
    private static final double[] A = {
        0.0,
        -24345.0,
        -8720288.849,
         150512.5253,
        -156.6950325,
         476470.3222,
         729482.8271,
        -145421.402,
         2931.1506,
        -40.427932,
         5106.192,
         15711.36,
        -155011.1084,
         4360.53352,
         12.9492344,
         10236.884,
         13176.786,
        -326669.5104,
         7390.68412,
        -27.8986976,
         16643.076,
         30988.146
    };

    private static final double[] LB = { 0.0, 1.2, 20.0, 9.0, 6.5 };
    private static final double[] UB = { 1000.0, 2.4, 60.0, 9.3, 7.0 };

    // MODE=1 initial guess in Fortran TP84
    private static final double[] X0 = { 2.52, 2.0, 37.5, 9.25, 6.8 };

    /**
     * f(x) = -(a1 + x1*(a2 + a3*x2 + a4*x3 + a5*x4 + a6*x5)).
     */
    private static class TP84Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);

            return -(A[1] + x1 * (A[2] + A[3] * x2 + A[4] * x3 + A[5] * x4 + A[6] * x5));
        }

        @Override
        public RealVector gradient(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);

            final double[] g = new double[DIM];
            g[0] = -(A[2] + A[3] * x2 + A[4] * x3 + A[5] * x4 + A[6] * x5);
            g[1] = -A[3] * x1;
            g[2] = -A[4] * x1;
            g[3] = -A[5] * x1;
            g[4] = -A[6] * x1;
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(final RealVector x) {
            final RealMatrix h = new Array2DRowRealMatrix(DIM, DIM);
            h.setEntry(0, 1, -A[3]);
            h.setEntry(1, 0, -A[3]);
            h.setEntry(0, 2, -A[4]);
            h.setEntry(2, 0, -A[4]);
            h.setEntry(0, 3, -A[5]);
            h.setEntry(3, 0, -A[5]);
            h.setEntry(0, 4, -A[6]);
            h.setEntry(4, 0, -A[6]);
            return h;
        }
    }

    /**
     * Interval constraints from TP84 in inequality form:
     *
     *   0 <= V1 <= 294000
     *   0 <= V2 <= 294000
     *   0 <= V3 <= 277200
     */
    private static class TP84Ineq extends InequalityConstraint {

        TP84Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);

            final double v1 = x1 * (A[7] + A[8] * x2 + A[9] * x3 + A[10] * x4 + A[11] * x5);
            final double v2 = x1 * (A[12] + A[13] * x2 + A[14] * x3 + A[15] * x4 + A[16] * x5);
            final double v3 = x1 * (A[17] + A[18] * x2 + A[19] * x3 + A[20] * x4 + A[21] * x5);

            return new ArrayRealVector(new double[] {
                v1,
                v2,
                v3,
                294000.0 - v1,
                294000.0 - v2,
                277200.0 - v3
            }, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);

            final RealMatrix j = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // V1 derivatives
            j.setEntry(0, 0, A[7] + A[8] * x2 + A[9] * x3 + A[10] * x4 + A[11] * x5);
            j.setEntry(0, 1, A[8] * x1);
            j.setEntry(0, 2, A[9] * x1);
            j.setEntry(0, 3, A[10] * x1);
            j.setEntry(0, 4, A[11] * x1);

            // V2 derivatives
            j.setEntry(1, 0, A[12] + A[13] * x2 + A[14] * x3 + A[15] * x4 + A[16] * x5);
            j.setEntry(1, 1, A[13] * x1);
            j.setEntry(1, 2, A[14] * x1);
            j.setEntry(1, 3, A[15] * x1);
            j.setEntry(1, 4, A[16] * x1);

            // V3 derivatives
            j.setEntry(2, 0, A[17] + A[18] * x2 + A[19] * x3 + A[20] * x4 + A[21] * x5);
            j.setEntry(2, 1, A[18] * x1);
            j.setEntry(2, 2, A[19] * x1);
            j.setEntry(2, 3, A[20] * x1);
            j.setEntry(2, 4, A[21] * x1);

            // Upper interval bounds: -(Vi derivatives)
            for (int c = 0; c < DIM; c++) {
                j.setEntry(3, c, -j.getEntry(0, c));
                j.setEntry(4, c, -j.getEntry(1, c));
                j.setEntry(5, c, -j.getEntry(2, c));
            }

            return j;
        }
    }

    @Test
    public void testHS084Optimization() {
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        final LagrangeSolution sol = optimizer.optimize(
            new InitialGuess(X0),
            new ObjectiveFunction(new TP84Obj()),
            new TP84Ineq(),
            new SimpleBounds(LB, UB)
        );

        // FEX in TP84: -0.528033513306D+07
        HSProblemTestUtils.assertExpectedObjective(-5280335.13306, sol);
    }
}