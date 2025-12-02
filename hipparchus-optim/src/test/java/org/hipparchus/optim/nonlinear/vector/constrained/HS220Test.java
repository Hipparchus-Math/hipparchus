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

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * HS220 (TP220)
 * --------------
 *
 * N    = 2 variables
 * NILI = 0
 * NINL = 1  (1 nonlinear inequality)
 * NELI = 0
 * NENL = 0
 *
 * Initial guess (mode 1):
 *   X(1) = 0.25D+5 = 25000
 *   X(2) = 0.25D+5 = 25000
 *
 * Bounds (mode 1):
 *   LXL(1) = TRUE, XL(1) = 1.0  → x1 >= 1
 *   LXL(2) = TRUE, XL(2) = 0.0  → x2 >= 0
 *   LXU(1..2) = FALSE           → no upper bounds
 *
 * Objective:
 *   f(x) = x1
 *
 * Gradient:
 *   df/dx1 = 1
 *   df/dx2 = 0
 *
 * Nonlinear inequality:
 *   G1(x) = (x1 - 1)^3 - x2 <= 0
 *
 *   ∂G1/∂x1 = 3 * (x1 - 1)^2
 *   ∂G1/∂x2 = -1
 *
 * Exact solution (LEX = .TRUE.):
 *   XEX(1) = 1.0
 *   XEX(2) = 0.0
 *   FEX    = 1.0
 */
public class HS220Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 1;
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective function
    // -------------------------------------------------------------------------
    private static class HS220Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            return x1;
        }

        @Override
        public RealVector gradient(RealVector x) {
            // df/dx = [1, 0]
            return new ArrayRealVector(new double[] {1.0, 0.0}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Linear objective -> Hessian = 0
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraint: G1(x) = (x1 - 1)^3 - x2 <= 0
    // -------------------------------------------------------------------------
    private static class HS220Ineq extends InequalityConstraint {

        HS220Ineq() {
            // RHS = 0 for the single inequality
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double g1 = FastMath.pow(x1 - 1.0, 3.0) - x2;

            return new ArrayRealVector(new double[] {g1}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            final double x1 = x.getEntry(0);

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // ∂G1/∂x1 = 3 * (x1 - 1)^2
            // ∂G1/∂x2 = -1
            J.setEntry(0, 0, 3.0 * FastMath.pow(x1 - 1.0, 2.0));
            J.setEntry(0, 1, -1.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS220_optimization() {

        // Initial guess (Fortran mode 1): X(i) = 0.25D+5 = 25000
        double[] x0 = new double[] {25000.0, 25000.0};

        // Bounds: x1 >= 1, x2 >= 0
        double[] lower = new double[] {1.0, 0.0};
        double[] upper = new double[] {
                Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY
        };
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
//                new InitialGuess(x0),
                new ObjectiveFunction(new HS220Obj()),
                null,              // no equality constraints
                new HS220Ineq(),   // 1 nonlinear inequality
                bounds
        );

        double f = sol.getValue();

        // Exact optimum (LEX = .TRUE.)
        final double fExpected = 1.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
