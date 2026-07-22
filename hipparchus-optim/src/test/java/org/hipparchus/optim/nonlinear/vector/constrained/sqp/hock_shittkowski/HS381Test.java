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
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * HS381 (TP381) – Linear objective with linear constraints (13 variables).
 *
 * From TP381:
 *
 *   N    = 13
 *   NILI = 3   (3 linear inequality constraints)
 *   NINL = 0
 *   NELI = 1   (1 linear equality constraint)
 *   NENL = 0
 *
 * Decision variables: x1..x13.
 *
 * Objective:
 *   FX = sum_{i=1..13} R(i) * X(i)
 *
 * with
 *   R = {
 *     0.8, 1.1, 0.85, 3.45, 2.0, 2.1, 3.0,
 *     0.8, 0.45, 0.72, 1.8, 3.0, 0.6
 *   }.
 *
 * Linear constraints (MODE = 4 in TP381):
 *
 *   G1(x) = sum_i S(i) * X(i) - 18.0      (inequality)
 *   G2(x) = sum_i U(i) * X(i) - 10.0      (inequality)
 *   G3(x) = sum_i V(i) * X(i) - 0.9       (inequality)
 *   G4(x) = sum_i X(i) - 10.0             (equality)
 *
 * with coefficient vectors:
 *
 *   S = {
 *     11.6, 13.7,  9.5, 48.5, 31.9, 51.1, 65.5,
 *      0.0,  0.0,  0.0, 21.8, 46.9,  0.0
 *   }
 *
 *   U = {
 *      0.05, 0.07, 0.0, 0.33, 0.0, 1.27, 1.27,
 *     23.35, 35.84, 0.81, 1.79, 7.34, 0.0
 *   }
 *
 *   V = {
 *      0.35, 0.37, 0.10, 0.62, 0.0, 1.03, 1.69,
 *     18.21, 0.01, 0.08, 0.31, 1.59, 22.45
 *   }
 *
 * Bounds:
 *   LXL(i) = TRUE, XL(i) = 0.0, LXU(i) = FALSE
 *   ⇒ 0 ≤ x_i, no upper bound.
 *
 * Reference solution (LEX = .FALSE.):
 *   FEX = 1.0149
 *   XEX ≈ (0.7864989, 0, 0, 0, 0, 0.17105694, 0, 0, 0.020676337,
 *          0, 0, 0, 0.019883712)
 *
 * We use FEX as an upper bound on f: require f ≤ FEX + tol.
 */
public class HS381Test {

    private static final int DIM      = 13;
    private static final int NUM_INEQ = 3;
    private static final int NUM_EQ   = 1;

    // R, S, U, V arrays from DATA statements
    private static final double[] R = {
        0.8, 1.1, 0.85, 3.45, 2.0, 2.1, 3.0,
        0.8, 0.45, 0.72, 1.8, 3.0, 0.6
    };

    private static final double[] S = {
        11.6, 13.7,  9.5, 48.5, 31.9, 51.1, 65.5,
         0.0,  0.0,  0.0, 21.8, 46.9,  0.0
    };

    private static final double[] U = {
        0.05, 0.07, 0.0, 0.33, 0.0, 1.27, 1.27,
        23.35, 35.84, 0.81, 1.79, 7.34, 0.0
    };

    private static final double[] V = {
        0.35, 0.37, 0.10, 0.62, 0.0, 1.03, 1.69,
        18.21, 0.01, 0.08, 0.31, 1.59, 22.45
    };

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS381Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double fx = 0.0;
            for (int i = 0; i < DIM; i++) {
                fx += R[i] * x.getEntry(i);
            }
            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            // Gradient is constant R
            return new ArrayRealVector(R, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Objective is linear → Hessian = 0
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Equality constraint: G4(x) = sum_i x_i - 10 = 0
    // -------------------------------------------------------------------------
    private static class HS381Eq extends EqualityConstraint {

        HS381Eq() {
            // Single equality with RHS = 0
            super(new ArrayRealVector(new double[NUM_EQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {
            double sum = 0.0;
            for (int i = 0; i < DIM; i++) {
                sum += x.getEntry(i);
            }
            double g4 = sum - 1.0;
            return new ArrayRealVector(new double[] { g4 }, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            RealMatrix J = new Array2DRowRealMatrix(NUM_EQ, DIM);
            // ∂G4/∂x_i = 1
            for (int j = 0; j < DIM; j++) {
                J.setEntry(0, j, 1.0);
            }
            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints: G1, G2, G3
    //
    // G1(x) = sum_i S(i)*x_i - 18.0
    // G2(x) = sum_i U(i)*x_i - 10.0
    // G3(x) = sum_i V(i)*x_i -  0.9
    //
    // The original Fortran stores exactly these expressions in G(1..3).
    // -------------------------------------------------------------------------
    private static class HS381Ineq extends InequalityConstraint {

        HS381Ineq() {
            // Right-hand sides are 0 for all inequalities
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {
            double sumS = 0.0;
            double sumU = 0.0;
            double sumV = 0.0;
            for (int i = 0; i < DIM; i++) {
                double xi = x.getEntry(i);
                sumS += S[i] * xi;
                sumU += U[i] * xi;
                sumV += V[i] * xi;
            }

            double g1 = sumS - 18.0;
            double g2 = sumU - 1.0;
            double g3 = sumV - 0.9;

            return new ArrayRealVector(new double[] { g1, g2, g3 }, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // Row 0: ∂G1/∂x_i = S(i)
            for (int j = 0; j < DIM; j++) {
                J.setEntry(0, j, S[j]);
            }

            // Row 1: ∂G2/∂x_i = U(i)
            for (int j = 0; j < DIM; j++) {
                J.setEntry(1, j, U[j]);
            }

            // Row 2: ∂G3/∂x_i = V(i)
            for (int j = 0; j < DIM; j++) {
                J.setEntry(2, j, V[j]);
            }

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS381_optimization() {

        // Initial guess: X(i) = 0.1 (MODE = 1)
        double[] x0 = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            x0[i] = 0.1;
        }

        // Bounds: XL(i) = 0, LXL(i)=TRUE, LXU(i)=FALSE → x_i >= 0
        double[] lower = new double[DIM];
        double[] upper = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            lower[i] = 0.0;
            upper[i] = Double.POSITIVE_INFINITY;
        }
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS381Obj()),
                new HS381Eq(),      // 1 equality
                new HS381Ineq(),    // 3 inequalities
                bounds
        );

        double f = sol.getValue();

        // LEX = .FALSE., FEX = 1.0149 → use as upper bound reference
        final double fExpected = 1.0149;
        final double tolF      = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        assertTrue(fExpected + tolF >= f,
                   "HS381: expected F <= " + (fExpected + tolF) + " but got F = " + f);
    }
}
