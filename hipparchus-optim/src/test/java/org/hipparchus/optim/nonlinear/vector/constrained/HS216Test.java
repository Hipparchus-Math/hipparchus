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


/**
 * HS216 (TP216)
 * -------------------------
 *
 * N    = 2 variables
 * NILI = 0
 * NINL = 0
 * NELI = 0
 * NENL = 1  (1 nonlinear equality constraint)
 *
 * Initial guess:
 *   X(1) = -1.2
 *   X(2) =  1.0
 *
 * Bounds (both variables):
 *   XL(i) = -3.0
 *   XU(i) = 10.0   for i = 1,2
 *
 * Objective:
 *   f(x) = 100 * (x1^2 - x2)^2 + (x1 - 1)^2
 *
 * Gradient:
 *   df/dx1 = 400 * (x1^2 - x2) * x1 + 2 * (x1 - 1)
 *   df/dx2 = -200 * (x1^2 - x2)
 *
 * Nonlinear equality constraint:
 *   G1(x) = x1*(x1 - 4) - 2*x2 + 12 = 0
 *
 * Jacobian:
 *   dG1/dx1 = 2*x1 - 4
 *   dG1/dx2 = -2
 *
 * Exact solution (correcting a typo in the Fortran listing):
 *   Fortran has XEX(1) = 2.0D+01, but to obtain FEX = 1.0 it must be 2.0:
 *
 *   x* = (2, 4)
 *   f* = 1
 *
 * Check:
 *   f(2,4) = 100 * (4 - 4)^2 + (2 - 1)^2 = 1
 *   G1(2,4) = 2*(2 - 4) - 2*4 + 12 = 0
 */
public class HS216Test {

    private static final int DIM = 2;
    private static final int NUM_EQ = 1;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS216Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            return 100.0 * FastMath.pow(x1 * x1 - x2, 2) +
                   FastMath.pow(x1 - 1.0, 2);
        }

        @Override
        public RealVector gradient(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double t = x1 * x1 - x2;
            final double df1 = 400.0 * t * x1 + 2.0 * (x1 - 1.0);
            final double df2 = -200.0 * t;
            return new ArrayRealVector(new double[] { df1, df2 }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Nonlinear equality constraint: G1(x) = 0
    // -------------------------------------------------------------------------
    private static class HS216Eq extends EqualityConstraint {

        HS216Eq() {
            // RHS = 0 for the single equality
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

            // G1(x) = x1*(x1 - 4) - 2*x2 + 12
            double g1 = x1 * (x1 - 4.0) - 2.0 * x2 + 12.0;

            return new ArrayRealVector(new double[] { g1 }, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            final double x1 = x.getEntry(0);
            final RealMatrix J = new Array2DRowRealMatrix(NUM_EQ, DIM);
            J.setEntry(0, 0, 2.0 * x1 - 4.0);
            J.setEntry(0, 1, -2.0);
            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS216_optimization() {

        // Initial guess from Fortran
        double[] x0 = new double[] { -1.2, 1.0 };

        // Bounds: -3 <= x_i <= 10, i = 1,2
        double[] lower = new double[] { -3.0, -3.0 };
        double[] upper = new double[] { 10.0, 10.0 };
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                HSProblemTestUtils.newExternalOption(),
                new InitialGuess(x0),
                new ObjectiveFunction(new HS216Obj()),
                new HS216Eq(),
                bounds
        );

        HSProblemTestUtils.assertExpectedObjective(1.0, sol);
    }
}