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
 * HS226 (TP226)
 *
 * N    = 2
 * NILI = 0
 * NINL = 2 (nonlinear inequalities)
 * NELI = 0
 * NENL = 0
 *
 * Bounds (from MODE=1):
 *   x1 >= 0, x2 >= 0, no upper bounds.
 *
 * Objective:
 *   f(x) = -x1 * x2
 *
 * Inequality constraints G(1..2):
 *   G1(x) = x1^2 + x2^2
 *   G2(x) = 1 - x1^2 - x2^2
 *
 * Reference (LEX = TRUE):
 *   x*  = (1/sqrt(2), 1/sqrt(2))
 *   f*  = -0.5
 */
public class HS226Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 2;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS226Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            return -x1 * x2;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            return new ArrayRealVector(new double[]{
                -x2,
                -x1
            }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);
            // f = -x1*x2 ⇒ ∂²f/∂x1∂x2 = ∂²f/∂x2∂x1 = -1
            H.setEntry(0, 1, -1.0);
            H.setEntry(1, 0, -1.0);
            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints G1..G2
    // -------------------------------------------------------------------------
    private static class HS226Ineq extends InequalityConstraint {

        HS226Ineq() {
            // RHS = 0 for both inequalities; we return G(x) as in the Fortran code.
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            double g1 = x1 * x1 + x2 * x2;
            double g2 = 1.0 - x1 * x1 - x2 * x2;

            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // From MODE=5:
            // GG(1,1) =  2*x1
            // GG(1,2) =  2*x2
            // GG(2,1) = -2*x1
            // GG(2,2) = -2*x2

            // Row 0: grad G1
            J.setEntry(0, 0,  2.0 * x1);
            J.setEntry(0, 1,  2.0 * x2);

            // Row 1: grad G2
            J.setEntry(1, 0, -2.0 * x1);
            J.setEntry(1, 1, -2.0 * x2);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS226_optimization() {

        // Initial guess from MODE=1
        double[] x0 = new double[]{0.8, 0.05};

        // Bounds: x1 >= 0, x2 >= 0, no upper bounds
        double[] lower = new double[]{0.0, 0.0};
        double[] upper = new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS226Obj()),
                null,
                new HS226Ineq(),
                bounds
        );

        double f = sol.getValue();

        // Reference optimum
        final double val = -0.5;
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
