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
 * HS219 (TP219)
 * --------------
 *
 * N    = 4 variables
 * NILI = 0
 * NINL = 0
 * NELI = 0
 * NENL = 2  (2 nonlinear inequalities)
 *
 * Initial guess:
 *   X(1..4) = 10
 *
 * Bounds:
 *   All variables free (no lower/upper bounds in the Fortran setup).
 *
 * Objective:
 *   f(x) = -x1
 *
 * Gradient:
 *   df/dx1 = -1
 *   df/dx2 =  0
 *   df/dx3 =  0
 *   df/dx4 =  0
 *
 * Inequality constraints:
 *
 *   G1(x) = x2 - x1^3 - x3^2 <= 0
 *     ∂G1/∂x1 = -3*x1^2
 *     ∂G1/∂x2 =  1
 *     ∂G1/∂x3 = -2*x3
 *     ∂G1/∂x4 =  0
 *
 *   G2(x) = x1^2 - x2 - x4^2 <= 0
 *     ∂G2/∂x1 =  2*x1
 *     ∂G2/∂x2 = -1
 *     ∂G2/∂x3 =  0
 *     ∂G2/∂x4 = -2*x4
 *
 * Exact solution (LEX = .TRUE.):
 *   XEX(1) = 1
 *   XEX(2) = 1
 *   XEX(3) = 0
 *   XEX(4) = 0
 *   FEX    = -1
 */
public class HS219Test {

    private static final int DIM      = 4;
    private static final int NUM_INEQ = 2;
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective function
    // -------------------------------------------------------------------------
    private static class HS219Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            return -x1;
        }

        @Override
        public RealVector gradient(RealVector x) {
            // df/dx = [-1, 0, 0, 0]
            return new ArrayRealVector(new double[] {-1.0, 0.0, 0.0, 0.0}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Linear objective -> Hessian = 0
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints: G1, G2
    // -------------------------------------------------------------------------
    private static class HS219Ineq extends InequalityConstraint {

        HS219Ineq() {
            // RHS = 0 for both inequalities
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
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);

            // G1(x) = x2 - x1^3 - x3^2
            final double g1 = x2 - FastMath.pow(x1, 3.0) - x3 * x3;

            // G2(x) = x1^2 - x2 - x4^2
            final double g2 = x1 * x1 - x2 - x4 * x4;

            return new ArrayRealVector(new double[] {g1, g2}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // Row 0: grad G1
            // ∂G1/∂x1 = -3*x1^2
            // ∂G1/∂x2 =  1
            // ∂G1/∂x3 = -2*x3
            // ∂G1/∂x4 =  0
            J.setEntry(0, 0, -3.0 * x1 * x1);
            J.setEntry(0, 1, 1.0);
            J.setEntry(0, 2, -2.0 * x3);
            J.setEntry(0, 3, 0.0);

            // Row 1: grad G2
            // ∂G2/∂x1 =  2*x1
            // ∂G2/∂x2 = -1
            // ∂G2/∂x3 =  0
            // ∂G2/∂x4 = -2*x4
            J.setEntry(1, 0, 2.0 * x1);
            J.setEntry(1, 1, -1.0);
            J.setEntry(1, 2, 0.0);
            J.setEntry(1, 3, -2.0 * x4);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS219_optimization() {

        // Initial guess (Fortran mode 1): X(i) = 10
        double[] x0 = new double[] {10.0, 10.0, 10.0, 10.0};

        // No bounds (all variables free)
        double[] lower = new double[] {
                Double.NEGATIVE_INFINITY,
                Double.NEGATIVE_INFINITY,
                Double.NEGATIVE_INFINITY,
                Double.NEGATIVE_INFINITY
        };
        double[] upper = new double[] {
                Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY
        };
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS219Obj()),
                null,              // no equality constraints
                new HS219Ineq(),   // 2 nonlinear inequalities
                bounds
        );

        double f = sol.getValue();

        // Exact optimum (LEX = .TRUE.)
        final double fExpected = -1.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
