/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 * https://www.apache.org/licenses/LICENSE-2.0
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

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.hipparchus.util.FastMath;

/**
 * Test case for Hock and Schittkowski problem 366 (HS366).
 * Objective: Quadratic / Non-linear (7 variables, 14 inequality constraints, all bounds finite).
 * Expected minimum value: 704.30560.
 */
public class HS366Test {

    private static final int DIM = 7;
    private static final int NUM_INEQ = 14;

    // Constants array C (38 elements, 0-indexed in Java)
    private static final double[] C = {
        0.59553571e-3, 0.88392857e+0, -0.11756250e+0, 1.1088e+0,
        0.1303533e+0, -0.0066033e+0, 0.66173269e-3, 0.17239878e-1,
        -0.56595559e-2, -0.19120592e-1, 0.5685075e+2, 1.08702e+0,
        0.32175e+0, -0.03762e+0, 0.006198e+0, 0.24623121e+4,
        -0.25125634e+2, 0.16118996e+3, 5.0e+3, -0.48951e+6,
        0.44333333e+2, 0.33e+0, 0.022556e+0, -0.007595e+0,
        0.00061e+0, -0.5e-3, 0.819672e+0, 0.819672e+0,
        0.245e+5, -0.25e+3, 0.10204082e-1, 0.12244898e-4,
        0.625e-4, 0.625e-4, -0.7625e-4, 1.22e+0, 1.0e+0, -1.0e+0
    };

    /**
     * Objective function f(x) for HS366 (MODE=2).
     */
    private static class HS366Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0); // X(1)
            final double x2 = x.getEntry(1); // X(2)
            final double x3 = x.getEntry(2); // X(3)
            final double x5 = x.getEntry(4); // X(5)
            final double x6 = x.getEntry(5); // X(6)

            // FX = 1.715*X(1) + 0.035*X(1)*X(6) + 4.0565*X(3) + 10*X(2) + 3000 - 0.063*X(3)*X(5)
            return 1.715 * x1 + 0.035 * x1 * x6 + 4.0565 * x3 + 10.0 * x2 + 3000.0 - 0.063 * x3 * x5;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[DIM];
            final double x1 = x.getEntry(0);
            final double x3 = x.getEntry(2);
            final double x5 = x.getEntry(4);
            final double x6 = x.getEntry(5);

            // df/dx1: 1.715 + 0.035*X(6)
            g[0] = 1.715 + 0.035 * x6;
            // df/dx2: 10
            g[1] = 10.0;
            // df/dx3: 4.0565 - 0.063*X(5)
            g[2] = 4.0565 - 0.063 * x5;
            // df/dx4: 0
            g[3] = 0.0;
            // df/dx5: -0.063*X(3)
            g[4] = -0.063 * x3;
            // df/dx6: 0.035*X(1)
            g[5] = 0.035 * x1;
            // df/dx7: 0
            g[6] = 0.0;

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            RealMatrix h = new Array2DRowRealMatrix(DIM, DIM);

            // d²f/(dx1 dx6) = 0.035
            h.setEntry(0, 5, 0.035);
            h.setEntry(5, 0, 0.035);

            // d²f/(dx3 dx5) = -0.063
            h.setEntry(2, 4, -0.063);
            h.setEntry(4, 2, -0.063);

            // All other second derivatives are zero
            return h;
        }
    }

    /**
     * Inequality constraints G(i) <= 0 for HS366 (MODE=4).
     * The constraints are defined as 1 - ... <= 0.
     */
    private static class HS366Ineq extends InequalityConstraint {

        HS366Ineq() {
            // All constraints G(i) <= 0, so RHS is a zero vector
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            final double x1 = x.getEntry(0); // X(1)
            final double x2 = x.getEntry(1); // X(2)
            final double x3 = x.getEntry(2); // X(3)
            final double x4 = x.getEntry(3); // X(4)
            final double x5 = x.getEntry(4); // X(5)
            final double x6 = x.getEntry(5); // X(6)
            final double x7 = x.getEntry(6); // X(7)

            double[] g = new double[NUM_INEQ];

            // G(1) = 1 - C(1)*X(6)^2 - C(2)*X(3)/X(1) - C(3)*X(6)
            g[0] = 1.0 - C[0] * FastMath.pow(x6, 2.0) - C[1] * x3 / x1 - C[2] * x6;

            // G(2) = 1 - C(4)*X(1)/X(3) - C(5)*X(1)/X(3)*X(6) - C(6)*X(1)/X(3)*X(6)^2
            g[1] = 1.0 - C[3] * x1 / x3 - C[4] * x1 / x3 * x6 - C[5] * x1 / x3 * FastMath.pow(x6, 2.0);

            // G(3) = 1 - C(7)*X(6)^2 - C(8)*X(5) - C(9)*X(4) - C(10)*X(6)
            g[2] = 1.0 - C[6] * FastMath.pow(x6, 2.0) - C[7] * x5 - C[8] * x4 - C[9] * x6;

            // G(4) = 1 - C(11)/X(5) - C(12)/X(5)*X(6) - C(13)*X(4)/X(5) - C(14)/X(5)*X(6)^2
            g[3] = 1.0 - C[10] / x5 - C[11] / x5 * x6 - C[12] * x4 / x5 - C[13] / x5 * FastMath.pow(x6, 2.0);

            // G(5) = 1 - C(15)*X(7) - C(16)*X(2)/X(3)/X(4) - C(17)*X(2)/X(3)
            g[4] = 1.0 - C[14] * x7 - C[15] * x2 / x3 / x4 - C[16] * x2 / x3;

            // G(6) = 1 - C(18)/X(7) - C(19)*X(2)/X(3)/X(7) - C(20)*X(2)/X(3)/X(4)/X(7)
            g[5] = 1.0 - C[17] / x7 - C[18] * x2 / x3 / x7 - C[19] * x2 / x3 / x4 / x7;

            // G(7) = 1 - C(21)/X(5) - C(22)*X(7)/X(5)
            g[6] = 1.0 - C[20] / x5 - C[21] * x7 / x5;

            // G(8) = 1 - C(23)*X(5) - C(24)*X(7)
            g[7] = 1.0 - C[22] * x5 - C[23] * x7;

            // G(9) = 1 - C(25)*X(3) - C(26)*X(1)
            g[8] = 1.0 - C[24] * x3 - C[25] * x1;

            // G(10) = 1 - C(27)*X(1)/X(3) - C(28)/X(3)
            g[9] = 1.0 - C[26] * x1 / x3 - C[27] / x3;

            // G(11) = 1 - C(29)*X(2)/X(3)/X(4) - C(30)*X(2)/X(3)
            g[10] = 1.0 - C[28] * x2 / x3 / x4 - C[29] * x2 / x3;

            // G(12) = 1 - C(31)*X(4) - C(32)/X(2)*X(3)*X(4)
            g[11] = 1.0 - C[30] * x4 - C[31] / x2 * x3 * x4;

            // G(13) = 1 - C(33)*X(1)*X(6) - C(34)*X(1) - C(35)*X(3)
            g[12] = 1.0 - C[32] * x1 * x6 - C[33] * x1 - C[34] * x3;

            // G(14) = 1 - C(36)/X(1)*X(3) - C(37)/X(1) - C(38)*X(6)
            g[13] = 1.0 - C[35] / x1 * x3 - C[36] / x1 - C[37] * x6;

            return new ArrayRealVector(g, false);
        }

        /**
         * Calculates the Jacobian matrix using central difference (numerical approximation).
         */
        @Override
        public RealMatrix jacobian(RealVector x) {
            final double eps = 1.0e-6;
            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            for (int j = 0; j < DIM; j++) {
                double xj = x.getEntry(j);

                // x + eps e_j
                x.setEntry(j, xj + eps);
                RealVector gp = value(x);

                // x - eps e_j
                x.setEntry(j, xj - eps);
                RealVector gm = value(x);

                // restore x(j)
                x.setEntry(j, xj);

                // central difference derivative
                for (int i = 0; i < NUM_INEQ; i++) {
                    double dij = (gp.getEntry(i) - gm.getEntry(i)) / (2.0 * eps);
                    J.setEntry(i, j, dij);
                }
            }

            return J;
        }
    }

    @Test
    public void testHS366_optimization() {

        // X(1) to X(7) initial guess (1745, 110, 3048, 89, 92.8, 8, 145)
        double[] x0 = {
            1745.0, // X(1)
            110.0,  // X(2)
            3048.0, // X(3)
            89.0,   // X(4)
            92.8,   // X(5)
            8.0,    // X(6)
            145.0   // X(7)
        };

        // Bounds (XL and XU from FORTRAN)
        double[] lower = new double[] {
            1.0,    // X(1) >= 1
            1.0,    // X(2) >= 1
            1.0,    // X(3) >= 1
            85.0,   // X(4) >= 85
            90.0,   // X(5) >= 90
            3.0,    // X(6) >= 3
            145.0   // X(7) >= 145
        };

        double[] upper = new double[] {
            2000.0, // X(1) <= 2000
            120.0,  // X(2) <= 120
            5000.0, // X(3) <= 5000
            93.0,   // X(4) <= 93
            95.0,   // X(5) <= 95
            12.0,   // X(6) <= 12
            162.0   // X(7) <= 162
        };

        SimpleBounds bounds = new SimpleBounds(lower, upper);

        // Optimization setup (using the provided SQPOptimizerS2 for compatibility)
        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new HS366Obj()),
            new HS366Ineq(),
            bounds
        );

        double f = sol.getValue();

        // FEX = 0.70430560D+03 (704.30560)
        final double fExpected = 704.30560;

        // Tolerance for objective function value
        final double tolF = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        // Assert that the found objective value is close to the expected minimum.
        assertEquals(fExpected, f, tolF, "HS366: objective mismatch");



    }
}
