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
 * HS218 (TP218)
 * --------------
 *
 * N    = 2 variables
 * NILI = 0
 * NINL = 1   (1 nonlinear inequality)
 * NELI = 0
 * NENL = 0
 *
 * Initial guess:
 *   X(1) = 9
 *   X(2) = 100
 *
 * Bounds:
 *   x1 free, x2 >= 0
 *
 * Objective:
 *   f(x) = x2
 *
 * Gradient:
 *   df/dx1 = 0
 *   df/dx2 = 1
 *
 * Inequality constraint:
 *   G1(x) = x2 - x1^2 <= 0
 *   ∂G1/∂x1 = -2*x1
 *   ∂G1/∂x2 =  1
 *
 * Exact solution (LEX = .TRUE.):
 *   XEX(1) = 0
 *   XEX(2) = 0
 *   FEX    = 0
 */
public class HS218Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 1;
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective function
    // -------------------------------------------------------------------------
    private static class HS218Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            final double x2 = x.getEntry(1);
            return x2;
        }

        @Override
        public RealVector gradient(RealVector x) {
            // df/dx1 = 0, df/dx2 = 1
            return new ArrayRealVector(new double[] {0.0, 1.0}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Hessian = 0 (linear objective)
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraint: G1(x) <= 0
    // -------------------------------------------------------------------------
    private static class HS218Ineq extends InequalityConstraint {

        HS218Ineq() {
            // Single inequality, RHS = 0
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

            // G1(x) = x2 - x1^2
            final double g1 = x2 - x1 * x1;

            return new ArrayRealVector(new double[] { g1 }, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            final double x1 = x.getEntry(0);
            // ∂G1/∂x1 = -2*x1
            // ∂G1/∂x2 =  1
            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);
            J.setEntry(0, 0, -2.0 * x1);
            J.setEntry(0, 1, 1.0);
            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS218_optimization() {

        // Initial guess from Fortran
        double[] x0 = new double[] {9.0, 100.0};

        // Bounds: x1 free, x2 >= 0
        double[] lower = new double[] {Double.NEGATIVE_INFINITY, 0.0};
        double[] upper = new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS218Obj()),
                null,              // no equality constraints
                new HS218Ineq(),   // 1 nonlinear inequality
                bounds
        );

        double f = sol.getValue();

        // Exact optimum (LEX = .TRUE.)
        final double fExpected = 0.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
