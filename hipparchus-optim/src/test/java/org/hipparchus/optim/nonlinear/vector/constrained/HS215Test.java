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
 * HS215 (TP215)
 * -------------------------
 *
 * N     = 2 variables
 * NILI  = 0
 * NINL  = 1   (1 nonlinear inequality constraint)
 * NELI  = 0
 * NENL  = 0
 *
 * Initial guess:
 *   X(1) = 1
 *   X(2) = 1
 *
 * Bounds:
 *   x1 >= 0   (LXL(1)=TRUE, XL(1)=0)
 *   x2 free   (LXL(2)=FALSE)
 *
 * Inequality constraint (G1(x) <= 0):
 *   G1(x) = x2 - x1^2
 *
 * Jacobian:
 *   dG1/dx1 = -2*x1
 *   dG1/dx2 = 1
 *
 * Objective:
 *   f(x) = x2
 *
 * Exact solution (LEX = TRUE):
 *   x* = (0, 0)
 *   f* = 0
 */
public class HS215Test {

    private static final int DIM = 2;
    private static final int NUM_INEQ = 1;

    // -------------------------------------------------------------------------
    // Objective function
    // -------------------------------------------------------------------------
    private static class HS215Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() { return DIM; }

        @Override
        public double value(RealVector x) {
            // FX = X(2)
            return x.getEntry(1);
        }

        @Override
        public RealVector gradient(RealVector x) {
            // GF(1) = 0
            // GF(2) = 1
            return new ArrayRealVector(new double[]{0.0, 1.0}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            return new Array2DRowRealMatrix(DIM, DIM); // Start with zero Hessian
        }
    }

    // -------------------------------------------------------------------------
    // Single nonlinear inequality G(x) <= 0
    // -------------------------------------------------------------------------
    private static class HS215Ineq extends InequalityConstraint {

        HS215Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ])); // RHS = 0
        }

        @Override
        public int dim() { return DIM; }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // G1 = x2 - x1^2
            return new ArrayRealVector(new double[]{x2 - x1 * x1}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            double x1 = x.getEntry(0);

            // dG1/dx1 = -2*x1
            // dG1/dx2 =  1
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
    public void testHS215_optimization() {

        // Initial guess
        double[] x0 = new double[]{1.0, 1.0};

        // Bounds: x1 >= 0, x2 free
        double[] lower = new double[]{0.0, Double.NEGATIVE_INFINITY};
        double[] upper = new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS215Obj()),
                null,              // no equality constraints
                new HS215Ineq(),   // 1 nonlinear inequality
                bounds
        );

        double f = sol.getValue();

        // Exact optimum (LEX = .TRUE.)
        double val = 0.0;
        

        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
