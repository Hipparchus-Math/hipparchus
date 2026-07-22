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
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * HS377 (TP377) – 10-variable problem with 3 linear equality constraints.
 *
 * From TP377:
 *
 *   N     = 10
 *   NILI  = 0
 *   NINL  = 0
 *   NELI  = 3
 *   NENL  = 0
 *
 * Bounds:
 *   For i = 1..10:
 *     LXL(i) = .TRUE., LXU(i) = .TRUE.
 *     XL(i)  = 0.1D-3 = 1e-4
 *     XU(i)  = 0.1D+2 = 10.0
 *
 * Initial point:
 *   X(i) = 0.1  for all i
 *
 * Objective:
 *
 *   Given A(1..10):
 *     A = (-6.089, -17.164, -34.054, -5.914,
 *          -24.721, -14.986, -24.100, -10.708,
 *          -26.662, -22.179)
 *
 *   SUM = sum_{i=1..10} x_i
 *
 *   FX = sum_{i=1..10} x_i * ( A(i) + log( max( x_i / SUM, 1e-5 ) ) )
 *
 * Gradient (as in TP377 – MODE=3):
 *
 *   SUM = sum x_i
 *   GF(i) = A(i) + log( x_i / SUM )
 *
 * Note: il gradiente Fortran NON include né il clipping max(..,1e-5)
 * né la derivata di SUM: qui lo replichiamo esattamente, senza "correggerlo".
 *
 * Equality constraints (G=0):
 *
 *   G1 = x1 - 2 x2 + 2 x3 + x6 + x10 - 2 = 0
 *   G2 = x4 - 2 x5 + x6 + x7 - 1        = 0
 *   G3 = x3 + x7 + x8 + 2 x9 + x10 - 1  = 0
 *
 * (GG in MODE=1 definisce la Jacobiana costante).
 *
 * Reference:
 *   LEX = .FALSE.
 *   FEX = -795.001
 *   Quindi usiamo FEX come upper bound: f(x) <= FEX + tol.
 */
public class HS377Test {

    private static final int DIM   = 10;
    private static final int NUM_EQ = 3;

    // Coefficienti A(i) dalla DATA
    private static final double[] A = {
        -6.089,   // A(1)
        -17.164,  // A(2)
        -34.054,  // A(3)
        -5.914,   // A(4)
        -24.721,  // A(5)
        -14.986,  // A(6)
        -24.100,  // A(7)
        -10.708,  // A(8)
        -26.662,  // A(9)
        -22.179   // A(10)
    };

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS377Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        /**
         * FX = sum_i x_i * ( A(i) + log( max(x_i / SUM, 1e-5) ) )
         * con SUM = sum_i x_i
         */
        @Override
        public double value(RealVector x) {
            double sum = 0.0;
            for (int i = 0; i < DIM; i++) {
                sum += x.getEntry(i);
            }

            double fx = 0.0;
            for (int i = 0; i < DIM; i++) {
                double xi = x.getEntry(i);
                double ratio = xi / sum;
                double clipped = FastMath.max(ratio, 1.0e-5);
                fx += xi * (A[i] + FastMath.log(clipped));
            }
            return fx;
        }

        /**
         * GF(i) = A(i) + log( x_i / SUM )
         * come nel MODE=3 del Fortran (senza clipping e senza derivata di SUM).
         */
        @Override
        public RealVector gradient(RealVector x) {
            double sum = 0.0;
            for (int i = 0; i < DIM; i++) {
                sum += x.getEntry(i);
            }

            double[] g = new double[DIM];
            for (int i = 0; i < DIM; i++) {
                double xi = x.getEntry(i);
                double ratio = xi / sum;
                g[i] = A[i] + FastMath.log(ratio);
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Non fornita in TP377; usiamo Hessian nullo
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Equality constraints (3 linear constraints)
    // -------------------------------------------------------------------------
    private static class HS377Eq extends EqualityConstraint {

        HS377Eq() {
            // RHS = 0 for all 3 constraints
            super(new ArrayRealVector(new double[NUM_EQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        /**
         * G1 = x1 - 2 x2 + 2 x3 + x6 + x10 - 2
         * G2 = x4 - 2 x5 + x6 + x7 - 1
         * G3 = x3 + x7 + x8 + 2 x9 + x10 - 1
         */
        @Override
        public RealVector value(RealVector x) {
            final double x1  = x.getEntry(0);
            final double x2  = x.getEntry(1);
            final double x3  = x.getEntry(2);
            final double x4  = x.getEntry(3);
            final double x5  = x.getEntry(4);
            final double x6  = x.getEntry(5);
            final double x7  = x.getEntry(6);
            final double x8  = x.getEntry(7);
            final double x9  = x.getEntry(8);
            final double x10 = x.getEntry(9);

            double g1 = x1 - 2.0 * x2 + 2.0 * x3 + x6 + x10 - 2.0;
            double g2 = x4 - 2.0 * x5 + x6 + x7 - 1.0;
            double g3 = x3 + x7 + x8 + 2.0 * x9 + x10 - 1.0;

            return new ArrayRealVector(new double[] { g1, g2, g3 }, false);
        }

        /**
         * Jacobiana costante, come GG in MODE=1:
         *
         * G1: [ 1, -2,  2, 0, 0, 1, 0, 0, 0, 1 ]
         * G2: [ 0,  0,  0, 1,-2, 1, 1, 0, 0, 0 ]
         * G3: [ 0,  1,  1, 0, 0, 0, 1, 1, 2, 1 ]
         */
        @Override
        public RealMatrix jacobian(RealVector x) {
            RealMatrix J = new Array2DRowRealMatrix(NUM_EQ, DIM);

            // Row 0: G1
            J.setEntry(0, 0,  1.0);
            J.setEntry(0, 1, -2.0);
            J.setEntry(0, 2,  2.0);
            J.setEntry(0, 3,  0.0);
            J.setEntry(0, 4,  0.0);
            J.setEntry(0, 5,  1.0);
            J.setEntry(0, 6,  0.0);
            J.setEntry(0, 7,  0.0);
            J.setEntry(0, 8,  0.0);
            J.setEntry(0, 9,  1.0);

            // Row 1: G2
            J.setEntry(1, 0,  0.0);
            J.setEntry(1, 1,  0.0);
            J.setEntry(1, 2,  0.0);
            J.setEntry(1, 3,  1.0);
            J.setEntry(1, 4, -2.0);
            J.setEntry(1, 5,  1.0);
            J.setEntry(1, 6,  1.0);
            J.setEntry(1, 7,  0.0);
            J.setEntry(1, 8,  0.0);
            J.setEntry(1, 9,  0.0);

            // Row 2: G3
            J.setEntry(2, 0,  0.0);
            J.setEntry(2, 1,  1.0);
            J.setEntry(2, 2,  1.0);
            J.setEntry(2, 3,  0.0);
            J.setEntry(2, 4,  0.0);
            J.setEntry(2, 5,  0.0);
            J.setEntry(2, 6,  1.0);
            J.setEntry(2, 7,  1.0);
            J.setEntry(2, 8,  2.0);
            J.setEntry(2, 9,  1.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS377_optimization() {

        // Initial point: X(i) = 0.1
        double[] x0 = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            x0[i] = 0.1;
        }

        // Bounds: 1e-4 <= x_i <= 10.0
        double[] lower = new double[DIM];
        double[] upper = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            lower[i] = 1.0e-4;
            upper[i] = 10.0;
        }
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS377Obj()),
                new HS377Eq(),   // 3 equality constraints
                bounds
        );

        double f = sol.getValue();

        // LEX = .FALSE., FEX = -795.001 → usato come upper bound su f
        final double fExpected = -795.001;
        final double tolF      = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        assertTrue(fExpected + tolF >= f,
                   "HS377: expected F <= " + (fExpected + tolF) + " but got F = " + f);
    }
}
