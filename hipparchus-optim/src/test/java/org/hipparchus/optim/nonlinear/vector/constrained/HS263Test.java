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
 * HS263 (TP263)
 *
 * N    = 4
 * NILI = 0
 * NINL = 2  (2 nonlinear inequalities)
 * NELI = 0
 * NENL = 2  (2 nonlinear equalities)
 *
 * Objective:
 *   f(x) = -x1
 *
 * Inequalities (treated as G >= 0):
 *   G1 = x2 - x1^3             >= 0
 *   G2 = x1^2 - x2             >= 0
 *
 * Equalities:
 *   G3 = x2 - x1^3 - x3^2      = 0
 *   G4 = x1^2 - x2 - x4^2      = 0
 *
 * Initial guess:
 *   x = (10, 10, 10, 10)
 *
 * Reference (Fortran):
 *   x*  = (10, 10, 0, 0)
 *   f*  = -10
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

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS263Test {

    private static final int DIM       = 4;
    private static final int NUM_INEQ  = 2;
    private static final int NUM_EQ    = 2;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS263Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            return -x.getEntry(0); // -x1
        }

        @Override
        public RealVector gradient(RealVector x) {
            return new ArrayRealVector(new double[]{
                    -1.0,  // df/dx1
                     0.0,
                     0.0,
                     0.0
            }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Linear objective
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Nonlinear inequalities (G >= 0)
    // -------------------------------------------------------------------------
    private static class HS263Ineq extends InequalityConstraint {

        HS263Ineq() {
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

            double g1 = x2 - x1 * x1 * x1; // x2 - x1^3
            double g2 = x1 * x1 - x2;      // x1^2 - x2

            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            double x1 = x.getEntry(0);

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // g1 = x2 - x1^3
            J.setEntry(0, 0, -3.0 * x1 * x1); // ∂g1/∂x1
            J.setEntry(0, 1,  1.0);           // ∂g1/∂x2
            J.setEntry(0, 2,  0.0);           // ∂g1/∂x3
            J.setEntry(0, 3,  0.0);           // ∂g1/∂x4

            // g2 = x1^2 - x2
            J.setEntry(1, 0,  2.0 * x1);      // ∂g2/∂x1
            J.setEntry(1, 1, -1.0);           // ∂g2/∂x2
            J.setEntry(1, 2,  0.0);
            J.setEntry(1, 3,  0.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Nonlinear equalities (G = 0)
    // -------------------------------------------------------------------------
    private static class HS263Eq extends EqualityConstraint {

        HS263Eq() {
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

            double g3 = x2 - x1 * x1 * x1 - x3 * x3; // x2 - x1^3 - x3^2
            double g4 = x1 * x1 - x2 - x4 * x4;      // x1^2 - x2 - x4^2

            return new ArrayRealVector(new double[]{g3, g4}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            RealMatrix J = new Array2DRowRealMatrix(NUM_EQ, DIM);

            // g3 = x2 - x1^3 - x3^2
            J.setEntry(0, 0, -3.0 * x1 * x1);  // ∂g3/∂x1
            J.setEntry(0, 1,  1.0);            // ∂g3/∂x2
            J.setEntry(0, 2, -2.0 * x3);       // ∂g3/∂x3
            J.setEntry(0, 3,  0.0);            // ∂g3/∂x4

            // g4 = x1^2 - x2 - x4^2
            J.setEntry(1, 0,  2.0 * x1);       // ∂g4/∂x1
            J.setEntry(1, 1, -1.0);            // ∂g4/∂x2
            J.setEntry(1, 2,  0.0);            // ∂g4/∂x3
            J.setEntry(1, 3, -2.0 * x4);       // ∂g4/∂x4

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS263() {

        // Initial guess: x = (10, 10, 10, 10)
        double[] x0 = new double[]{10.0, 10.0, 10.0, 10.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS263Obj()),
                new HS263Eq(),     // 2 equalities
                new HS263Ineq(),   // 2 inequalities
                null               // no bounds
        );

        double f = sol.getValue();

        double fExpected = -1.0;
        double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
