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

package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;

import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * HS385 (TP385) – 15-variable problem with 10 nonlinear inequality constraints.
 *
 * From TP385:
 *
 *   N    = 15
 *   NILI = 0
 *   NINL = 10
 *   NELI = 0
 *   NENL = 0
 *
 * Objective (MODE=2):
 *   FX = - sum_{i=1..15} D(i) * X(i)
 *
 * Nonlinear inequalities (MODE=4):
 *   For i = 1..10:
 *     C_i = sum_{j=1..15} A(i,j) * X(j)^2
 *     G(i) = B(i) - C_i   (G <= 0)
 *
 * No bounds:
 *   LXL(i) = FALSE, LXU(i) = FALSE → unbounded variables.
 *
 * Reference:
 *   LEX = .FALSE.
 *   FEX = -0.83152859D+4 = -8315.2859
 */
public class HS385Test {

    private static final int DIM      = 15;
    private static final int NUM_INEQ = 10;

    // B(1..10) from TP385
    private static final double[] B = {
        3.85e2, 4.70e2, 5.60e2, 5.65e2, 6.45e2,
        4.30e2, 4.85e2, 4.55e2, 8.90e2, 4.60e2
    };

    // D(1..15) – same as in TP384
    private static final double[] D = {
        4.86e2, 6.40e2, 7.58e2, 7.76e2, 4.77e2,
        7.07e2, 1.75e2, 6.19e2, 6.27e2, 6.14e2,
        4.75e2, 3.77e2, 5.24e2, 4.68e2, 5.29e2
    };

    /**
     * A(10,15) ricostruita esattamente dal DATA A Fortran (ordine colonna-major).
     * Indici qui: A[row][col] con row = 0..9 (vincolo 1..10), col = 0..14 (x1..x15).
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
        { 10.0,  70.0,  10.0, 35.0, 25.0, 65.0,  0.0, 30.0, 500.0,  0.0,
          25.0,   0.0,  15.0, 50.0, 55.0},
        // row 10 = A(10,1..15)
        {  5.0,  10.0, 100.0,  5.0, 20.0,  5.0, 10.0, 35.0,  95.0, 70.0,
          20.0,  10.0,  35.0, 10.0, 30.0}
    };

    // -------------------------------------------------------------------------
    // Objective: f(x) = - Σ D(i) * x(i)
    // -------------------------------------------------------------------------
    private static class HS385Obj extends TwiceDifferentiableFunction {

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
            // Linear objective → zero Hessian
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints: 10 nonlinear constraints
    // G_i(x) = B(i) - Σ_j A(i,j) * x_j^2  ≤ 0
    // -------------------------------------------------------------------------
    private static class HS385Ineq extends InequalityConstraint {

        HS385Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {
            double[] g = new double[NUM_INEQ];

            for (int i = 0; i < NUM_INEQ; i++) {
                double c = 0.0;
                for (int j = 0; j < DIM; j++) {
                    double xj = x.getEntry(j);
                    c += A[i][j] * xj * xj;
                }
                g[i] = B[i] - c;
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            for (int i = 0; i < NUM_INEQ; i++) {
                for (int j = 0; j < DIM; j++) {
                    double xj = x.getEntry(j);
                    J.setEntry(i, j, -2.0 * A[i][j] * xj);
                }
            }

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS385_optimization() {

        // Initial guess: X(i) = 0.0 (as in MODE=1)
        double[] x0 = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            x0[i] = 0.0;
        }

        

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS385Obj()),
                new HS385Ineq()
               
        );

        double f = sol.getValue();

        final double fExpected = -0.83152859e4;
        final double tolF      = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        // LEX = .FALSE. → FEX is an upper bound: FEX >= f
        assertTrue(fExpected + tolF >= f,
                   "HS385: expected F <= " + (fExpected + tolF) + " but got F = " + f);
    }
}
