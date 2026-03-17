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
 * HS217 (TP217)
 * --------------
 *
 * N    = 2 variables
 * NILI = 1   (1 linear inequality)
 * NINL = 0
 * NELI = 0
 * NENL = 1   (1 nonlinear equality)
 *
 * Initial guess:
 *   X(1) = 10
 *   X(2) = 10
 *
 * Bounds:
 *   x1 >= 0, x2 free
 *
 * Objective:
 *   f(x) = -x2
 *
 * Gradient:
 *   df/dx1 = 0
 *   df/dx2 = -1
 *
 * Inequality constraint (linear):
 *   G1(x) = 1 + x1 - 2*x2 <= 0
 *   ∂G1/∂x1 =  1
 *   ∂G1/∂x2 = -2
 *
 * Nonlinear equality constraint:
 *   G2(x) = x1^2 + x2^2 - 1 = 0
 *   ∂G2/∂x1 = 2*x1
 *   ∂G2/∂x2 = 2*x2
 *
 * Exact solution (LEX = .TRUE.):
 *   XEX(1) = 0.6
 *   XEX(2) = 0.8
 *   FEX    = -0.8
 */
public class HS217Test {

    private static final int DIM     = 2;
    private static final int NUM_EQ  = 1;
    private static final int NUM_INEQ = 1;

    // -------------------------------------------------------------------------
    // Objective function
    // -------------------------------------------------------------------------
    private static class HS217Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            final double x2 = x.getEntry(1);
            return -x2;
        }

        @Override
        public RealVector gradient(RealVector x) {
            // df/dx1 = 0, df/dx2 = -1
            return new ArrayRealVector(new double[] {0.0, -1.0}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Hessian = 0 (linear objective)
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Equality constraint: G2(x) = 0
    // -------------------------------------------------------------------------
    private static class HS217Eq extends EqualityConstraint {

        HS217Eq() {
            // Single equality, RHS = 0
            super(new ArrayRealVector(new double[NUM_EQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            // G2(x) = x1^2 + x2^2 - 1
            final double g2 = x1 * x1 + x2 * x2 - 1.0;

            return new ArrayRealVector(new double[] { g2 }, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            // ∂G2/∂x1 = 2*x1
            // ∂G2/∂x2 = 2*x2
            RealMatrix J = new Array2DRowRealMatrix(NUM_EQ, DIM);
            J.setEntry(0, 0, 2.0 * x1);
            J.setEntry(0, 1, 2.0 * x2);
            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraint: G1(x) <= 0
    // -------------------------------------------------------------------------
    private static class HS217Ineq extends InequalityConstraint {

        HS217Ineq() {
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

            // G1(x) = 1 + x1 - 2*x2
            final double g1 = 1.0 + x1 - 2.0 * x2;

            return new ArrayRealVector(new double[] { g1 }, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            // ∂G1/∂x1 =  1
            // ∂G1/∂x2 = -2
            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);
            J.setEntry(0, 0, 1.0);
            J.setEntry(0, 1, -2.0);
            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS217_optimization() {

        // Initial guess from Fortran
        double[] x0 = new double[] {10.0, 10.0};

        // Bounds: x1 >= 0, x2 free
        double[] lower = new double[] {0.0, Double.NEGATIVE_INFINITY};
        double[] upper = new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS217Obj()),
                new HS217Eq(),     // 1 nonlinear equality
                new HS217Ineq(),   // 1 linear inequality
                bounds
        );

        double f = sol.getValue();

        // Exact optimum (LEX = .TRUE.)
        final double fExpected = -0.8;
        final double tol = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
