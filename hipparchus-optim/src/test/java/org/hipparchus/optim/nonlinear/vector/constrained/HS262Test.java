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
 * HS262 (TP262)
 *
 * N    = 4
 * NILI = 3  (3 linear inequalities)
 * NINL = 0
 * NELI = 1  (1 linear equality)
 * NENL = 0
 *
 * Objective:
 *   f(x) = -0.5 x1 - x2 - 0.5 x3 - x4
 *
 * Inequalities (Fortran G(i) >= 0):
 *   G1 = 10 - x1 - x2 - x3 - x4                 >= 0
 *   G2 = 10 - 0.2 x1 - 0.5 x2 - x3 - 2 x4      >= 0
 *   G3 = 10 - 2 x1 - x2 - 0.5 x3 - 0.2 x4      >= 0
 *
 * Equality:
 *   G4 = x1 + x2 + x3 - 2 x4 - 6               = 0
 *
 * Bounds:
 *   x1 >= 0, x2 >= 0, x3 >= 0, x4 >= 0   (no upper bounds)
 *
 * Initial guess:
 *   x = (1, 1, 1, 1)
 *
 * Reference solution (Fortran):
 *   x*  ≈ (0, 26/3, 0, 4/3)
 *   f*  = -10
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

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS262Test {

    private static final int DIM       = 4;
    private static final int NUM_INEQ  = 3;
    private static final int NUM_EQ    = 1;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS262Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);
            return -0.5 * x1 - x2 - 0.5 * x3 - x4;
        }

        @Override
        public RealVector gradient(RealVector x) {
            // Constant gradient
            return new ArrayRealVector(new double[]{
                    -0.5,  // df/dx1
                    -1.0,  // df/dx2
                    -0.5,  // df/dx3
                    -1.0   // df/dx4
            }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Linear objective: Hessian is zero
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequalities G >= 0 (3 linear)
    // -------------------------------------------------------------------------
    private static class HS262Ineq extends InequalityConstraint {

        HS262Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ])); // RHS = 0
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            double g1 = 10.0 - x1 - x2 - x3 - x4;
            double g2 = 10.0 - 0.2 * x1 - 0.5 * x2 - x3 - 2.0 * x4;
            double g3 = 10.0 - 2.0 * x1 - x2 - 0.5 * x3 - 0.2 * x4;

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // g1 = 10 - x1 - x2 - x3 - x4
            J.setEntry(0, 0, -1.0);
            J.setEntry(0, 1, -1.0);
            J.setEntry(0, 2, -1.0);
            J.setEntry(0, 3, -1.0);

            // g2 = 10 - 0.2 x1 - 0.5 x2 - x3 - 2 x4
            J.setEntry(1, 0, -0.2);
            J.setEntry(1, 1, -0.5);
            J.setEntry(1, 2, -1.0);
            J.setEntry(1, 3, -2.0);

            // g3 = 10 - 2 x1 - x2 - 0.5 x3 - 0.2 x4
            J.setEntry(2, 0, -2.0);
            J.setEntry(2, 1, -1.0);
            J.setEntry(2, 2, -0.5);
            J.setEntry(2, 3, -0.2);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Equality G = 0 (1 linear equality)
    // -------------------------------------------------------------------------
    private static class HS262Eq extends EqualityConstraint {

        HS262Eq() {
            super(new ArrayRealVector(new double[NUM_EQ])); // RHS = 0
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            // g4 = x1 + x2 + x3 - 2 x4 - 6 = 0
            double g4 = x1 + x2 + x3 - 2.0 * x4 - 6.0;

            return new ArrayRealVector(new double[]{g4}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_EQ, DIM);

            // ∂g4/∂x = [1, 1, 1, -2]
            J.setEntry(0, 0,  1.0);
            J.setEntry(0, 1,  1.0);
            J.setEntry(0, 2,  1.0);
            J.setEntry(0, 3, -2.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS262() {

        // Initial guess: x = (1, 1, 1, 1)
        double[] x0 = new double[]{1.0, 1.0, 1.0, 1.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Bounds: x_i >= 0 (no upper bound)
        SimpleBounds bounds = new SimpleBounds(
                new double[]{0.0, 0.0, 0.0, 0.0},
                new double[]{
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY,
                        Double.POSITIVE_INFINITY
                }
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS262Obj()),
                new HS262Eq(),     // 1 equality
                new HS262Ineq(),   // 3 inequalities
                bounds
        );

        double f = sol.getValue();

        double fExpected = -10.0;
        double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
