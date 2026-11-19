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
 * HS389 (TP389) – Quadratic + linear inequality portfolio-type problem.
 *
 * N    = 15 variables
 * NILI = 4   (4 linear inequality constraints)
 * NINL = 11  (10 quadratic inequality constraints + 1 nonlinear inequality)
 * NELI = 0
 * NENL = 0
 *
 * Objective:
 *   f(x) = - sum_{j=1..15} D_j * x_j
 *
 * Inequality constraints (Fortran G(1..15)):
 *
 * Linear constraints (L = 1..4, Fortran indices I = 12..15):
 *   For L = 1..4:
 *     C_L(x) = sum_j A1(L,j) * x_j
 *     G(L)   = B(11 + L) - C_L(x)
 *
 * Quadratic constraints (i = 1..10):
 *   For i = 1..10:
 *     C_i(x)      = sum_j A(i,j) * x_j^2
 *     G(i+4)(x)   = B(i) - C_i(x)
 *
 * Last nonlinear constraint (index 15):
 *   C_15(x) = sum_j j * (x_j - 2)^2
 *   G(15)   = C_15(x)/2 - 2.D+2 = C_15(x)/2 - 200
 *
 * In the Java wrapper, Gk(x) are returned directly and the base class
 * InequalityConstraint supplies RHS = 0, so the solver enforces Gk(x) <= 0.
 *
 * Reference:
 *   FEX = -0.58097197D+4  (LEX = .FALSE. → only an upper bound, FEX >= f*)
 */
public class HS389Test {

    private static final int DIM      = 15;
    private static final int NUM_INEQ = 15;
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Data (decoded once from Fortran DATA, same as TP388)
    // -------------------------------------------------------------------------

    /** Quadratic coefficients A(i,j) for constraints 5..14 (10×15). */
    private static final double[][] A_MAT = new double[][] {
        {100.0, 100.0, 100.0,  35.0,  10.0,  10.0,   0.0,   0.0,   0.0,  10.0,   0.0,   0.0,   0.0,  10.0,   0.0},
        { 90.0,  50.0,  60.0,  55.0,  20.0,  20.0,  40.0,  40.0,  55.0,  25.0,  10.0,  10.0,  45.0,   0.0,   5.0},
        { 70.0,   0.0,  30.0,  65.0,  30.0,  25.0,  30.0,  60.0,   0.0,  25.0,  10.0,  45.0,  30.0,  10.0,  15.0},
        { 50.0,  10.0,  40.0,  60.0,  40.0,  35.0,   0.0,  60.0,  65.0,  35.0,  70.0,   0.0,  35.0,  30.0,   0.0},
        { 50.0,   0.0,  10.0,  90.0,  10.0,  45.0,  40.0,  35.0,   0.0,  50.0,   0.0,   0.0,  55.0,  35.0,   0.0},
        { 40.0,  60.0, 100.0,  95.0,  50.0,  50.0,  40.0,   0.0,  80.0,   0.0,  50.0,  10.0,   0.0,  20.0,   5.0},
        { 30.0,  30.0,   5.0,  90.0,  20.0,  10.0,  35.0,  75.0,   0.0,  95.0,  10.0,  10.0,  10.0,  25.0,  15.0},
        { 20.0,   0.0,  25.0,  25.0,  25.0,  15.0,   0.0,  15.0,   5.0,  10.0,  40.0,  45.0,  20.0,  20.0,  20.0},
        { 10.0,   0.0,  35.0,  35.0,  10.0,  50.0,  10.0,  20.0,  10.0,  30.0,  65.0,   0.0,  30.0,   0.0,  55.0},
        {  5.0,   0.0,   0.0,   5.0,  35.0,  10.0,  10.0,  30.0,   0.0,  30.0,  50.0,   5.0,  30.0,   5.0,  30.0}
    };

    /** Linear coefficients A1(l,j) for constraints 1..4 (4×15). */
    private static final double[][] A1_MAT = new double[][] {
        // first row: 1, 2, 3, ... 15
        {  1.0,   2.0,   3.0,   4.0,   5.0,   6.0,   7.0,   8.0,   9.0,  10.0, 11.0, 12.0, 13.0, 14.0, 15.0},
        { 45.0,  25.0,  35.0,  40.0,  50.0,  73.0,  17.0,  52.0,  86.0, 14.0, 30.0, 37.0, 17.0, 52.0, 86.0},
        { 53.0,  74.0,  26.0,  58.0,  25.0,  25.0,  26.0,  24.0,  85.0, 35.0, 14.0, 95.0, 26.0, 24.0, 85.0},
        { 12.0,  43.0,  51.0,  58.0,  60.0,  42.0,  60.0,  20.0,  40.0, 80.0, 75.0, 18.0, 60.0, 20.0, 67.0}
    };

    /**
     * B vector as in Fortran DATA B:
     * B(1..10) used in quadratic constraints,
     * B(12..15) used in linear constraints; B(11) = 0.
     */
    private static final double[] B_VEC = new double[] {
        385.0, 470.0, 560.0, 565.0, 645.0,
        430.0, 485.0, 455.0, 390.0, 460.0,
          0.0,  70.0, 361.0, 265.0, 395.0
    };

    /** D vector in f(x) = -sum D_j x_j. */
    private static final double[] D_VEC = new double[] {
        486.0, 640.0, 758.0, 776.0, 477.0,
        707.0, 175.0, 619.0, 627.0, 614.0,
        475.0, 377.0, 524.0, 468.0, 529.0
    };

    // -------------------------------------------------------------------------
    // Objective function
    // -------------------------------------------------------------------------

    private static class HS389Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double fx = 0.0;
            for (int j = 0; j < DIM; j++) {
                fx -= D_VEC[j] * x.getEntry(j);
            }
            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            // df/dx_j = -D_j (independent of x)
            double[] g = new double[DIM];
            for (int j = 0; j < DIM; j++) {
                g[j] = -D_VEC[j];
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Objective is linear → Hessian = 0
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints (15 constraints total)
    // -------------------------------------------------------------------------

    private static class HS389Ineq extends InequalityConstraint {

        HS389Ineq() {
            // RHS = 0 for all inequalities; we return G(x) directly.
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double[] g = new double[NUM_INEQ];

            // 1) Linear constraints #1..4: G(L) = B(11+L) - Σ A1(L,j) x_j
            for (int L = 0; L < 4; L++) {
                double c = 0.0;
                for (int j = 0; j < DIM; j++) {
                    c += A1_MAT[L][j] * x.getEntry(j);
                }
                g[L] = B_VEC[11 + L] - c;
            }

            // 2) Quadratic constraints #5..14: G(4+i) = B(i+1) - Σ A(i,j) x_j²
            for (int i = 0; i < 10; i++) {
                double c = 0.0;
                for (int j = 0; j < DIM; j++) {
                    double xj = x.getEntry(j);
                    c += A_MAT[i][j] * xj * xj;
                }
                g[4 + i] = B_VEC[i] - c;
            }

            // 3) Nonlinear constraint #15:
            //    C(x) = Σ (j+1)*(x_j - 2)²
            //    G(15) = C/2 - 200
            double C = 0.0;
            for (int j = 0; j < DIM; j++) {
                double diff = x.getEntry(j) - 2.0;
                C += (j + 1) * diff * diff;
            }
            g[14] = 0.5 * C - 200.0;

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // 1) Linear constraints 1..4: dG/dx_j = -A1(L,j)
            for (int L = 0; L < 4; L++) {
                for (int j = 0; j < DIM; j++) {
                    J.setEntry(L, j, -A1_MAT[L][j]);
                }
            }

            // 2) Quadratic constraints 5..14:
            //    G = B(i) - Σ A(i,j) x_j² → dG/dx_j = -2*A(i,j)*x_j
            for (int i = 0; i < 10; i++) {
                int row = 4 + i;
                for (int j = 0; j < DIM; j++) {
                    double xj = x.getEntry(j);
                    J.setEntry(row, j, -2.0 * A_MAT[i][j] * xj);
                }
            }

            // 3) Constraint #15:
            //    G(15) = 1/2 Σ (j+1)*(x_j-2)² → dG/dx_j = (j+1)*(x_j-2)
            int row15 = 14;
            for (int j = 0; j < DIM; j++) {
                double val = (j + 1) * (x.getEntry(j) - 2.0);
                J.setEntry(row15, j, val);
            }

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------

    @Test
    public void testHS389_optimization() {

        // Initial guess X(i) = 0.0 (as in Fortran TP389)
        double[] x0 = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            x0[i] = 0.0;
        }

        // No bounds: LXL(i) = LXU(i) = .FALSE.
        double[] lower = new double[DIM];
        double[] upper = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            lower[i] = Double.NEGATIVE_INFINITY;
            upper[i] = Double.POSITIVE_INFINITY;
        }
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new HS389Obj()),
            null,                  // no equality constraints
            new HS389Ineq(),       // 15 inequality constraints
            bounds
        );

        double f = sol.getValue();

        // LEX = .FALSE. → FEX is only an upper bound: FEX >= f*.
        final double val = -0.58097197e4;
       HSProblemTestUtils.assertBetterObjective(val, sol);
    }
}
