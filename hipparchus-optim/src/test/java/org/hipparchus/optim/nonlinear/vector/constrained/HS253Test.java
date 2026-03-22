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

/*
 * HS253 (TP253)
 *
 * N     = 3
 * NILI  = 1   (1 linear inequality)
 * NINL  = 0
 * NELI  = 0
 * NENL  = 0
 *
 * Constraint (Fortran G ≥ 0):
 *   G1 = 30 - 3*x1 - 3*x3 ≥ 0
 *
 * Objective:
 *   f(x) = sum_{j=1..8} sqrt[ (A1j - x1)^2 + (A2j - x2)^2 + (A3j - x3)^2 ]
 *
 * Bounds (MODE=1):
 *   x1 ≥ 0
 *   x2 ≥ 0
 *   x3 ≥ 0
 *
 * Initial guess (MODE=1):
 *   x = (0, 2, 0)
 *
 * Solution:
 *   x*  = (5,5,5)
 *   f*  = 69.282032
 */

package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.*;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS253Test {

    private static final int DIM      = 3;
    private static final int NUM_INEQ = 1;

   // Data matrix A(3,8) – CORRETTA
private static final double[][] A = new double[][]{
    {0.0, 10.0, 10.0,  0.0,  0.0, 10.0, 10.0,  0.0},  // row 1 (i=1)
    {0.0,  0.0, 10.0, 10.0,  0.0,  0.0, 10.0, 10.0},  // row 2 (i=2)
    {0.0,  0.0,  0.0,  0.0, 10.0, 10.0, 10.0, 10.0}   // row 3 (i=3)
};

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS253Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() { return DIM; }

        @Override
        public double value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double fx = 0.0;

            for (int j = 0; j < 8; j++) {
                double dx1 = A[0][j] - x1;
                double dx2 = A[1][j] - x2;
                double dx3 = A[2][j] - x3;
                fx += FastMath.sqrt(dx1*dx1 + dx2*dx2 + dx3*dx3);
            }

            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double g1 = 0.0;
            double g2 = 0.0;
            double g3 = 0.0;

            for (int j = 0; j < 8; j++) {

                double dx1 = x1 - A[0][j];
                double dx2 = x2 - A[1][j];
                double dx3 = x3 - A[2][j];

                double r = FastMath.sqrt(dx1*dx1 + dx2*dx2 + dx3*dx3);

                g1 += dx1 / r;
                g2 += dx2 / r;
                g3 += dx3 / r;
            }

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // true Hessian is complicated; return zero (BFGS will build it)
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints (G ≥ 0)
    // -------------------------------------------------------------------------
    private static class HS253Ineq extends InequalityConstraint {

        HS253Ineq() {
            super(new ArrayRealVector(new double[]{0.0}));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x3 = x.getEntry(2);

            // Fortran:
            // G(1) = 30 - 3*x1 - 3*x3  ≥ 0
            double g1 = 30.0 - 3.0 * x1 - 3.0 * x3;

            return new ArrayRealVector(new double[]{g1}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(1, DIM);

            // ∂g1/∂x1 = -3
            // ∂g1/∂x2 = 0
            // ∂g1/∂x3 = -3
            J.setEntry(0, 0, -3.0);
            J.setEntry(0, 1,  0.0);
            J.setEntry(0, 2, -3.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test runner
    // -------------------------------------------------------------------------
    @Test
    public void testHS253() {

        double[] x0 = new double[]{0.0, 2.0, 0.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Bounds:
        // x1 ≥ 0, x2 ≥ 0, x3 ≥ 0
        SimpleBounds bounds = new SimpleBounds(
                new double[]{0.0, 0.0, 0.0},
                new double[]{Double.POSITIVE_INFINITY,
                             Double.POSITIVE_INFINITY,
                             Double.POSITIVE_INFINITY}
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS253Obj()),
                null,               // no equalities
                new HS253Ineq(),    // 1 inequality
                bounds
        );

        double f = sol.getValue();

        double fExpected = 69.282032;
        double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
