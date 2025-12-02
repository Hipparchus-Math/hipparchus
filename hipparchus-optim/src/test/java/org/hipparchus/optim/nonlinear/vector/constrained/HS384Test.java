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

import org.hipparchus.util.FastMath;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * HS384 (TP384) – 15-variable problem with 10 nonlinear inequality constraints.
 *
 * From TP384:
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
 *     G(i) = B(i) - C_i
 *
 * In our wrapper we model them as G(x) <= 0 or >= 0 consistently
 * with the other HS3xx tests by using G(i) exactly as in Fortran.
 *
 * No bounds:
 *   LXL(i) = FALSE, LXU(i) = FALSE → unbounded variables.
 *
 * Reference:
 *   LEX = .FALSE.
 *   FEX = -0.83102590D+4 = -8310.259
 */
public class HS384Test {

    private static final int DIM      = 15;
    private static final int NUM_INEQ = 10;

    // B(1..10)
    private static final double[] B = {
        3.85e2, 4.7e2, 5.6e2, 5.65e2, 6.45e2,
        4.3e2,  4.85e2, 4.55e2, 3.9e2, 8.6e2
    };

    // D(1..15)
    private static final double[] D = {
        4.86e2, 6.4e2, 7.58e2, 7.76e2, 4.77e2,
        7.07e2, 1.75e2, 6.19e2, 6.27e2, 6.14e2,
        4.75e2, 3.77e2, 5.24e2, 4.68e2, 5.29e2
    };

    /**
     * A(10,15), ricostruita dalla DATA Fortran in ordine (i,j).
     *
     * In Fortran:
     *   DATA A/ ... / con A(10,15)
     * riempita in column-major; qui la riportiamo già come A[i][j]
     * con i=0..9 (constraint index 1..10), j=0..14 (var index 1..15).
     */
    private static final double[][] A = {
        // i = 1
        {100.0, 100.0, 10.0,  5.0, 10.0,  0.0,  0.0, 25.0,  0.0, 10.0,
          55.0,  5.0, 45.0, 20.0,  0.0},
        // i = 2
        { 90.0, 100.0, 10.0, 35.0, 20.0,  5.0,  0.0, 35.0, 55.0, 25.0,
          20.0,  0.0, 40.0, 25.0, 10.0},
        // i = 3
        { 70.0,  50.0,  0.0, 55.0, 25.0,100.0, 40.0, 50.0,  0.0, 30.0,
          60.0, 10.0, 30.0,  0.0, 40.0},
        // i = 4
        { 50.0,   0.0,  0.0, 65.0, 35.0,100.0, 35.0, 60.0,  0.0, 15.0,
           0.0, 75.0, 35.0, 30.0, 65.0},
        // i = 5
        { 50.0,  10.0, 70.0, 60.0, 45.0, 45.0,  0.0, 35.0, 65.0,  5.0,
          75.0,100.0, 75.0, 10.0,  0.0},
        // i = 6
        { 40.0,   0.0, 50.0, 95.0, 50.0, 35.0, 10.0, 60.0,  0.0, 45.0,
          15.0, 20.0,  0.0,  5.0,  5.0},
        // i = 7
        { 30.0,  60.0, 30.0, 90.0,  0.0, 30.0,  5.0, 25.0,  0.0, 70.0,
          20.0, 25.0, 70.0, 15.0, 15.0},
        // i = 8
        { 20.0,  30.0, 40.0, 25.0, 40.0, 25.0, 15.0, 10.0, 80.0, 20.0,
          30.0, 30.0,  5.0, 65.0, 20.0},
        // i = 9
        { 10.0,  70.0, 10.0, 35.0, 25.0, 65.0,  0.0, 30.0,  0.0,  0.0,
          25.0,  0.0, 15.0, 50.0, 55.0},
        // i = 10
        {  5.0,  10.0,500.0,  5.0, 20.0,  5.0, 10.0, 35.0, 95.0, 70.0,
          20.0, 10.0, 35.0, 10.0, 30.0}
    };

    // -------------------------------------------------------------------------
    // Objective: f(x) = - Σ D(i) * x(i)
    // -------------------------------------------------------------------------
    private static class HS384Obj extends TwiceDifferentiableFunction {

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
            // GF(i) = -D(i) (costante in MODE=3)
            for (int i = 0; i < DIM; i++) {
                g[i] = -D[i];
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Linear objective → Hessian = 0
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints: 10 nonlinear constraints G_i(x) = B(i) - Σ A(i,j)*x_j^2
    // -------------------------------------------------------------------------
    private static class HS384Ineq extends InequalityConstraint {

        HS384Ineq() {
            // Right-hand sides all 0 (G(x) ≤ 0 or ≥0 handled consistently in solver)
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
                double C = 0.0;
                for (int j = 0; j < DIM; j++) {
                    double xj = x.getEntry(j);
                    C += A[i][j] * xj * xj;
                }
                // Fortran: G(i) = B(i) - C
                g[i] = B[i] - C;
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            for (int i = 0; i < NUM_INEQ; i++) {
                for (int j = 0; j < DIM; j++) {
                    double xj = x.getEntry(j);
                    // GG(i,j) = -2*A(i,j)*X(j)
                    double dG = -2.0 * A[i][j] * xj;
                    J.setEntry(i, j, dG);
                }
            }

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS384_optimization() {

        // Initial guess: X(i) = 0.0 (as in MODE=1)
        double[] x0 = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            x0[i] = 0.0;
        }

        // No bounds (LXL/LXU = .FALSE. in Fortran)
        double[] lower = new double[DIM];
        double[] upper = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            lower[i] = Double.NEGATIVE_INFINITY;
            upper[i] = Double.POSITIVE_INFINITY;
        }
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS384Obj()),
                new HS384Ineq(),   // 10 nonlinear inequalities
                bounds
        );

        double f = sol.getValue();

        // LEX = .FALSE., FEX = -0.83102590D+4 → use as reference upper bound
        final double fExpected = -0.83102590e4;
        final double tolF      = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        assertTrue(fExpected + tolF >= f,
                   "HS384: expected F <= " + (fExpected + tolF) + " but got F = " + f);
    }
}
