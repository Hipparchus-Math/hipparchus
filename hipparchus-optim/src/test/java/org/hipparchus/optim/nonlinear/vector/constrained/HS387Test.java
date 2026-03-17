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
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * HS387 (TP387) – Quadratic inequality constrained problem.
 *
 * N = 15
 * NILI = 0, NINL = 11, NELI = 0, NENL = 0
 *
 * Objective:
 *   f(x) = - Σ_{i=1}^{15} D_i * x_i
 *
 * Inequality constraints (g(x) <= 0):
 *
 *   For i = 1..10:
 *     G_i(x) = B_i - Σ_{j=1}^{15} A_{i,j} * x_j^2 <= 0
 *
 *   For i = 11:
 *     G_11(x) = 0.5 * Σ_{j=1}^{15} j * (x_j - 2)^2 - 61 <= 0
 *
 * No bounds on x.
 *
 * Reference:
 *   LEX = .FALSE.
 *   FEX = -0.82501417D+4
 *   XEX (Fortran reference):
 *     [ 1.0125415D+1, 1.0158505D+1, 1.0309039D+1, 0.99697018D+1,
 *       0.98528372D+1, 1.0368532D+1, 0.99349349D+1, 0.97201160D+1,
 *       0.99994095D+1, 0.99547294D+1, 0.96953850D+1, 1.0080569D+1,
 *       0.98236999D+1, 0.99057993D+1, 0.97760168D+1 ]
 */
public class HS387Test {

    private static final int DIM = 15;
    private static final int NUM_INEQ = 11;

    // B(1..10)
    private static final double[] B = {
        3.85e2, 4.70e2, 5.60e2, 5.65e2, 6.45e2,
        4.30e2, 4.85e2, 4.55e2, 3.90e2, 4.60e2
    };

    // D(1..15)
    private static final double[] D = {
        4.86e2, 6.40e2, 7.58e2, 7.76e2, 4.77e2,
        7.07e2, 1.75e2, 6.19e2, 6.27e2, 6.14e2,
        4.75e2, 3.77e2, 5.24e2, 4.68e2, 5.29e2
    };

    /**
     * A(1..10,1..15) expanded from Fortran DATA (column-major to row-major mapping).
     * Each row here corresponds to fixed i, j=1..15.
     */
    private static final double[][] A = {
        // row 1 = A(1,1..15)
        {100.0, 100.0,  10.0,  5.0, 10.0,  0.0,  0.0, 25.0,   0.0, 10.0,
          55.0,   5.0,  45.0, 20.0,  0.0},
        // row 2 = A(2,1..15)
        { 90.0, 100.0,  10.0, 35.0, 20.0,  5.0,  0.0, 35.0,  55.0, 25.0,
          20.0,   0.0,  40.0, 25.0, 10.0},
        // row 3 = A(3,1..15)
        { 70.0,  50.0,   0.0, 55.0, 25.0,100.0, 40.0, 50.0,   0.0, 30.0,
          60.0,  10.0,  30.0,  0.0, 40.0},
        // row 4 = A(4,1..15)
        { 50.0,   0.0,   0.0, 65.0, 35.0,100.0, 35.0, 60.0,   0.0, 15.0,
           0.0,  75.0,  35.0, 30.0, 65.0},
        // row 5 = A(5,1..15)
        { 50.0,  10.0,  70.0, 60.0, 45.0, 45.0,  0.0, 35.0,  65.0,  5.0,
          75.0, 100.0,  75.0, 10.0,  0.0},
        // row 6 = A(6,1..15)
        { 40.0,   0.0,  50.0, 95.0, 50.0, 35.0, 10.0, 60.0,   0.0, 45.0,
          15.0,  20.0,   0.0,  5.0,  5.0},
        // row 7 = A(7,1..15)
        { 30.0,  60.0,  30.0, 90.0,  0.0, 30.0,  5.0, 25.0,   0.0, 70.0,
          20.0,  25.0,  70.0, 15.0, 15.0},
        // row 8 = A(8,1..15)
        { 20.0,  30.0,  40.0, 25.0, 40.0, 25.0, 15.0, 10.0,  80.0, 20.0,
          30.0,  30.0,   5.0, 65.0, 20.0},
        // row 9 = A(9,1..15)
        { 10.0,  70.0,  10.0, 35.0, 25.0, 65.0,  0.0, 30.0,   0.0,  0.0,
          25.0,   0.0,  15.0, 50.0, 55.0},
        // row 10 = A(10,1..15)
        {  5.0,  10.0, 100.0,  5.0, 20.0,  5.0, 10.0, 35.0,  95.0, 70.0,
          20.0,  10.0,  35.0, 10.0, 30.0}
    };

    // -------------------------------------------------------------------------
    // Objective function
    // -------------------------------------------------------------------------
    private static class HS387Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double fx = 0.0;
            for (int i = 0; i < DIM; i++) {
                fx -= D[i] * x.getEntry(i);
            }
            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[DIM];
            for (int i = 0; i < DIM; i++) {
                g[i] = -D[i];
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Hessian of a linear function is zero
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints (11 constraints)
    // -------------------------------------------------------------------------
    private static class HS387Ineq extends InequalityConstraint {

        HS387Ineq() {
            // All right-hand sides are 0 → g(x) <= 0
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        /**
         * G_i(x) <= 0
         *
         * For i = 0..9 (Fortran 1..10):
         *   G_i(x) = B_i - Σ_j A_{i,j} * x_j^2
         *
         * For i = 10 (Fortran 11):
         *   G_10(x) = 0.5 * Σ_j (j+1) * (x_j - 2)^2 - 61
         */
        @Override
        public RealVector value(RealVector x) {
            double[] g = new double[NUM_INEQ];

            // First 10 quadratic inequalities
            for (int i = 0; i < 10; i++) {
                double c = 0.0;
                for (int j = 0; j < DIM; j++) {
                    double xj = x.getEntry(j);
                    c += A[i][j] * xj * xj;
                }
                g[i] = B[i] - c;
            }

            // 11th inequality
            double cSum = 0.0;
            for (int j = 0; j < DIM; j++) {
                double xj = x.getEntry(j);
                double diff = xj - 2.0;
                cSum += (j + 1) * diff * diff; // DBLE(J)*(X(J)-2)**2
            }
            g[10] = 0.5 * cSum - 61.0;

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // Derivatives for i = 0..9:
            // ∂/∂x_j [B_i - Σ_k A_{i,k} x_k^2] = -2 A_{i,j} x_j
            for (int i = 0; i < 10; i++) {
                for (int j = 0; j < DIM; j++) {
                    double xj = x.getEntry(j);
                    J.setEntry(i, j, -2.0 * A[i][j] * xj);
                }
            }

            // Derivatives for i = 10 (11th constraint):
            // G_11 = 0.5 Σ_j (j+1) (x_j - 2)^2 - 61
            // ∂G_11/∂x_j = (j+1)*(x_j - 2)
            int row = 10;
            for (int j = 0; j < DIM; j++) {
                double xj = x.getEntry(j);
                double dG = (j + 1) * (xj - 2.0);
                J.setEntry(row, j, dG);
            }

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS387_optimization() {

        // Initial guess: X(I) = 0.0
        double[] x0 = new double[DIM];

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new HS387Obj()),
            new HS387Ineq()   // 11 nonlinear inequalities
        );

        double f = sol.getValue();

        // LEX = .FALSE. → FEX is a reference upper bound: FEX >= f (up to tolerance)
        final double fExpected = -0.82501417e4;
        final double tolF = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        assertTrue(fExpected + tolF >= f,
                   "HS387: expected F <= " + (fExpected + tolF) + " but got F = " + f);
    }
}
