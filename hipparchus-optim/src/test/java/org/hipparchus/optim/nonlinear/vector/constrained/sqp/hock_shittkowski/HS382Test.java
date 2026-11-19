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

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.Array2DRowRealMatrix;
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
 * HS382 (TP382) – Linear objective with 3 nonlinear inequality constraints
 * and 1 linear equality constraint (13 variables).
 *
 * From TP382:
 *
 *   N    = 13
 *   NILI = 0
 *   NINL = 3
 *   NELI = 1
 *   NENL = 0
 *
 * Objective:
 *   f(x) = sum_i R(i) * x(i)
 *
 * Inequalities (G1,G2,G3 <= 0):
 *
 *   Let
 *     q1 = sum_i Z1(i) * x(i)^2
 *     q2 = sum_i Z2(i) * x(i)^2
 *     q3 = sum_i Z3(i) * x(i)^2
 *
 *   G1(x) = q1 - 1.645 * sqrt(q1) + sum_i S(i)*x(i) - 18.0
 *   G2(x) = q2 - 1.645 * sqrt(q2) + sum_i U(i)*x(i) -  1.0
 *   G3(x) = q3 - 1.645 * sqrt(q3) + sum_i V(i)*x(i) -  0.9
 *
 * (riproduce MODE=4: attenzione 0.1645D+1 = 1.645 e 0.1D+1 = 1.0).
 *
 * Equality:
 *
 *   G4(x) = sum_i x(i) - 1.0 = 0
 *
 * Bounds:
 *   LXL(i)=TRUE, XL(i)=0; LXU(i)=FALSE → x(i) >= 0, no upper bound.
 *
 * Reference:
 *   LEX = .FALSE., FEX = 1.03831
 *   quindi richiediamo f <= FEX + tol.
 */
public class HS382Test {

    private static final int DIM      = 13;
    private static final int NUM_INEQ = 3;
    private static final int NUM_EQ   = 1;

    // R,S,U,V,Z1,Z2,Z3 from DATA
    private static final double[] R = {
        0.8, 1.1, 0.85, 3.45, 2.0, 2.1, 3.0,
        0.8, 0.45, 0.72, 1.8, 3.0, 0.6
    };

    private static final double[] S = {
        11.6, 13.7, 9.5, 48.5, 31.9, 51.1, 65.5,
        0.0, 0.0, 0.0, 21.8, 46.9, 0.0
    };

    private static final double[] U = {
        0.05, 0.07, 0.0, 0.33, 0.0, 1.27, 1.27,
        23.35, 35.84, 0.81, 1.79, 7.34, 0.0
    };

    private static final double[] V = {
        0.35, 0.37, 0.10, 0.62, 0.0, 1.03, 1.69,
        18.21, 0.01, 0.08, 0.31, 1.59, 22.45
    };

    private static final double[] Z1 = {
        0.4844, 0.3003, 0.1444, 0.0588, 4.9863,
        0.0653, 21.0222, 0.0, 0.0, 0.0, 0.2970, 9.2933, 0.0
    };

    private static final double[] Z2 = {
        0.0001, 0.0, 0.0, 0.0, 0.0,
        0.0040, 0.1404, 1.3631, 0.5138, 0.0289,
        0.0097, 0.3893, 0.0
    };

    private static final double[] Z3 = {
        0.001, 0.0009, 0.0001, 0.0005, 0.0,
        0.0021, 0.0825, 0.2073, 0.0, 0.0004,
        0.0005, 0.0107, 1.0206
    };

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS382Obj extends TwiceDifferentiableFunction {

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
            // In Fortran GF(i) = R(i) in MODE=1, no change in MODE=3 → grad costante
            return new ArrayRealVector(R, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Linear objective → Hessian = 0
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Equality constraint: G4(x) = sum x(i) - 1 = 0
    // -------------------------------------------------------------------------
    private static class HS382Eq extends EqualityConstraint {

        HS382Eq() {
            super(new ArrayRealVector(new double[NUM_EQ])); // RHS = 0
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
            for (int j = 0; j < DIM; j++) {
                J.setEntry(0, j, 1.0);
            }
            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints G1,G2,G3 <= 0, raggruppati in un'unica classe
    // -------------------------------------------------------------------------
    private static class HS382Ineq extends InequalityConstraint {

        HS382Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ])); // RHS = 0
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double q1 = 0.0;
            double q2 = 0.0;
            double q3 = 0.0;
            double sSum = 0.0;
            double uSum = 0.0;
            double vSum = 0.0;

            for (int i = 0; i < DIM; i++) {
                double xi = x.getEntry(i);
                double xi2 = xi * xi;

                q1 += Z1[i] * xi2;
                q2 += Z2[i] * xi2;
                q3 += Z3[i] * xi2;

                sSum += S[i] * xi;
                uSum += U[i] * xi;
                vSum += V[i] * xi;
            }

            double sqrt1 = FastMath.sqrt(q1);
            double sqrt2 = FastMath.sqrt(q2);
            double sqrt3 = FastMath.sqrt(q3);

            // G1,G2,G3 come da MODE=4 (0.1645D+1 = 1.645, 0.1D+1 = 1.0)
            double g1 = q1 - 1.645 * sqrt1 + sSum - 18.0;
            double g2 = q2 - 1.645 * sqrt2 + uSum - 1.0;
            double g3 = q3 - 1.645 * sqrt3 + vSum - 0.9;

            return new ArrayRealVector(new double[] { g1, g2, g3 }, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            double q1 = 0.0;
            double q2 = 0.0;
            double q3 = 0.0;
            for (int i = 0; i < DIM; i++) {
                double xi = x.getEntry(i);
                double xi2 = xi * xi;
                q1 += Z1[i] * xi2;
                q2 += Z2[i] * xi2;
                q3 += Z3[i] * xi2;
            }

            double sqrt1 = FastMath.sqrt(q1);
            double sqrt2 = FastMath.sqrt(q2);
            double sqrt3 = FastMath.sqrt(q3);

            // HELP fattori come in MODE=5:
            // HELP = -1.645D+0 / 2.D+0 / DSQRT(HELP_QUAD)
            double help1 = -1.645 / (2.0 * sqrt1);
            double help2 = -1.645 / (2.0 * sqrt2);
            double help3 = -1.645 / (2.0 * sqrt3);

            for (int i = 0; i < DIM; i++) {
                double xi = x.getEntry(i);

                // Fortran: GG(1,i) = S(i) + HELP * 2 * Z1(i) * X(i)
                double dG1dx = S[i] + help1 * 2.0 * Z1[i] * xi;
                double dG2dx = U[i] + help2 * 2.0 * Z2[i] * xi;
                double dG3dx = V[i] + help3 * 2.0 * Z3[i] * xi;

                J.setEntry(0, i, dG1dx);
                J.setEntry(1, i, dG2dx);
                J.setEntry(2, i, dG3dx);
            }

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // TEST
    // -------------------------------------------------------------------------
    @Test
    public void testHS382_optimization() {

        // Initial guess X(i) = 0.1
        double[] x0 = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            x0[i] = 0.1;
        }

        // Bounds: x >= 0
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
                new ObjectiveFunction(new HS382Obj()),
                new HS382Eq(),     // 1 equality
                new HS382Ineq(),   // 3 inequalities, raggruppate
                bounds
        );

        double f = sol.getValue();

        final double fExpected = 1.03831;
        final double tolF      = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        assertTrue(fExpected + tolF >= f,
                   "HS382: expected F <= " + (fExpected + tolF) + " but got F = " + f);
    }
}
