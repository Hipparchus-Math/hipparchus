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
 * HS374 (TP374) – 10-variable problem with 35 nonlinear inequality constraints.
 *
 * From TP374:
 *
 *   N     = 10
 *   NILI  = 0
 *   NINL  = 35
 *   NELI  = 0
 *   NENL  = 0
 *
 * Variables: x1..x10
 *
 * Objective:
 *   FX = X(10)
 *
 * Bounds:
 *   All variables free (no bounds) except for the usual test framework bounds;
 *   in TP374 all LXL/LXU are set to .FALSE., so we use no bounds in Java.
 *
 * Auxiliary functions (in Fortran outside TP374):
 *
 *   TP374A(a, x) = sum_{k=1..9} x_k cos(k a)
 *   TP374B(a, x) = sum_{k=1..9} x_k sin(k a)
 *   TP374G(a, x) = TP374A(a,x)^2 + TP374B(a,x)^2
 *
 * Constraints:
 *   Let G(a,x) = TP374G(a,x).
 *
 *   Group 1: i = 1..10
 *     z_i = (π/4) * (0.1 * (i-1))
 *     G_i(x) = G(z_i,x) - (1 - x10)^2 ≤ 0
 *
 *   Group 2: i = 11..20
 *     z_i = (π/4) * (0.1 * (i-11))
 *     G_i(x) = (1 + x10)^2 - G(z_i,x) ≤ 0
 *
 *   Group 3: i = 21..35
 *     z_i = (π/4) * (1.2 + 0.2*(i-21))
 *     G_i(x) = x10^2 - G(z_i,x) ≤ 0
 *
 * Reference:
 *   LEX = .FALSE.
 *   FEX = 0.233264D+0 = 0.233264
 *   So we only require  f(x) <= FEX + tol.
 */
public class HS374Test {

    private static final int DIM      = 10;
    private static final int NUM_INEQ = 35;

    // -------------------------------------------------------------------------
    // Objective function: f(x) = x10
    // -------------------------------------------------------------------------
    private static class HS374Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            // FX = X(10)
            return x.getEntry(9);
        }

        @Override
        public RealVector gradient(RealVector x) {
            // GF(1..9) = 0, GF(10) = 1.0
            double[] g = new double[DIM];
            g[9] = 1.0;
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Objective is linear in x10 → Hessian = 0
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints (35 nonlinear inequalities, G(x) <= 0)
    // -------------------------------------------------------------------------
    private static class HS374Ineq extends InequalityConstraint {

        HS374Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ])); // RHS = 0
        }

        @Override
        public int dim() {
            return DIM;
        }

        // ----------------- Helper functions (TP374A/B/G) ---------------------

        /** TP374A(a,x) = sum_{k=1..9} x_k cos(k a) */
        private double tp374A(double a, RealVector x) {
            double s = 0.0;
            for (int k = 1; k <= 9; k++) {
                s += x.getEntry(k - 1) * FastMath.cos(k * a);
            }
            return s;
        }

        /** TP374B(a,x) = sum_{k=1..9} x_k sin(k a) */
        private double tp374B(double a, RealVector x) {
            double s = 0.0;
            for (int k = 1; k <= 9; k++) {
                s += x.getEntry(k - 1) * FastMath.sin(k * a);
            }
            return s;
        }

        /** TP374G(a,x) = TP374A(a,x)^2 + TP374B(a,x)^2 */
        private double tp374G(double a, RealVector x) {
            double A = tp374A(a, x);
            double B = tp374B(a, x);
            return A * A + B * B;
        }

        // ---------------------------------------------------------------------
        // Constraint values
        // ---------------------------------------------------------------------
        @Override
        public RealVector value(RealVector x) {

            double[] g = new double[NUM_INEQ];

            final double x10 = x.getEntry(9);
            final double pi  = 4.0 * FastMath.atan(1.0);

            // Group 1: i = 1..10 (0-based idx 0..9)
            //   z = (π/4) * (0.1 * (i))
            //   G_i = G(z,x) - (1 - x10)^2
            for (int i = 0; i < 10; i++) {
                double z = (pi / 4.0) * (0.1 * i);
                double val = tp374G(z, x) - FastMath.pow(1.0 - x10, 2.0);
                g[i] = val;
            }

            // Group 2: i = 11..20 (0-based idx 10..19)
            //   z = (π/4) * (0.1 * (i-11))
            //   G_i = (1 + x10)^2 - G(z,x)
            for (int i = 10; i < 20; i++) {
                int idx = i - 10; // 0..9
                double z = (pi / 4.0) * (0.1 * idx);
                double val = FastMath.pow(1.0 + x10, 2.0) - tp374G(z, x);
                g[i] = val;
            }

            // Group 3: i = 21..35 (0-based idx 20..34)
            //   z = (π/4) * (1.2 + 0.2*(i-21))
            //   G_i = x10^2 - G(z,x)
            for (int i = 20; i < 35; i++) {
                int idx = i - 20; // 0..14
                double z = (pi / 4.0) * (1.2 + 0.2 * idx);
                double val = x10 * x10 - tp374G(z, x);
                g[i] = val;
            }

            return new ArrayRealVector(g, false);
        }

        // ---------------------------------------------------------------------
        // Jacobian
        // ---------------------------------------------------------------------
        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);
            final double x10 = x.getEntry(9);
            final double pi  = 4.0 * FastMath.atan(1.0);

            // Re-use A,B in derivatives:
            // TP374G(a,x) = A^2 + B^2
            // ∂G/∂x_k (k=1..9) = 2*A*cos(k a) + 2*B*sin(k a)
            // ∂G/∂x10 = 0

            // ---- Group 1: i = 1..10 (0-based 0..9) ----
            //
            // G_i(x) = G(z,x) - (1 - x10)^2
            //
            // ∂G_i/∂x_k (k=1..9) = ∂G/∂x_k
            // ∂G_i/∂x10         = 2*(1 - x10)
            //
            for (int i = 0; i < 10; i++) {
                double z = (pi / 4.0) * (0.1 * i);

                double A = tp374A(z, x);
                double B = tp374B(z, x);

                // k = 1..9 → x1..x9
                for (int k = 1; k <= 9; k++) {
                    double dGdxk = 2.0 * (A * FastMath.cos(k * z) + B * FastMath.sin(k * z));
                    J.setEntry(i, k - 1, dGdxk);
                }

                // d/dx10 of -(1 - x10)^2 = 2*(1 - x10)
                J.setEntry(i, 9, 2.0 * (1.0 - x10));
            }

            // ---- Group 2: i = 11..20 (0-based 10..19) ----
            //
            // G_i(x) = (1 + x10)^2 - G(z,x)
            //
            // ∂G_i/∂x_k (k=1..9) = - ∂G/∂x_k
            // ∂G_i/∂x10         = 2*(1 + x10)
            //
            for (int i = 10; i < 20; i++) {
                int idx = i - 10; // 0..9
                double z = (pi / 4.0) * (0.1 * idx);

                double A = tp374A(z, x);
                double B = tp374B(z, x);

                for (int k = 1; k <= 9; k++) {
                    double dGdxk = 2.0 * (A * FastMath.cos(k * z) + B * FastMath.sin(k * z));
                    J.setEntry(i, k - 1, -dGdxk);
                }

                // d/dx10 of (1 + x10)^2 = 2*(1 + x10)
                J.setEntry(i, 9, 2.0 * (1.0 + x10));
            }

            // ---- Group 3: i = 21..35 (0-based 20..34) ----
            //
            // G_i(x) = x10^2 - G(z,x)
            //
            // ∂G_i/∂x_k (k=1..9) = - ∂G/∂x_k
            // ∂G_i/∂x10         = 2*x10
            //
            for (int i = 20; i < 35; i++) {
                int idx = i - 20; // 0..14
                double z = (pi / 4.0) * (1.2 + 0.2 * idx);

                double A = tp374A(z, x);
                double B = tp374B(z, x);

                for (int k = 1; k <= 9; k++) {
                    double dGdxk = 2.0 * (A * FastMath.cos(k * z) + B * FastMath.sin(k * z));
                    J.setEntry(i, k - 1, -dGdxk);
                }

                J.setEntry(i, 9, 2.0 * x10);
            }

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS374_optimization() {

        // Initial guess from MODE = 1: X(i) = 0.1 for i=1..10
        double[] x0 = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            x0[i] = 0.1;
        }

        // No bounds (all LXL/LXU = .FALSE. in TP374)
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
                new ObjectiveFunction(new HS374Obj()),
                new HS374Ineq(),   // 35 inequality constraints
                bounds
        );

        double f = sol.getValue();

        // LEX = .FALSE. in TP374 → FEX is an upper bound, not exact:
        final double fExpected = 0.2332640;
        final double tolF      = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        assertTrue(fExpected + tolF >= f,
                   "HS374: expected F <= " + (fExpected + tolF) + " but got F = " + f);
    }
}
