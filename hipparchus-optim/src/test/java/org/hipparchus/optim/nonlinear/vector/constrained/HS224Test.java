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

/*
 * HS224 (TP224)
 *
 * N    = 2
 * NILI = 4  (4 linear inequality constraints)
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * Bounds:
 *   0 ≤ x1 ≤ 6
 *   0 ≤ x2 ≤ 6
 *
 * Objective:
 *   f(x) = 2*x1^2 + x2^2 - 48*x1 - 40*x2
 *
 * Constraints:
 *   G1: x1 + 3*x2
 *   G2: 18 - x1 - 3*x2
 *   G3: x1 + x2
 *   G4: 8 - x1 - x2
 *
 * Optimum (LEX = TRUE):
 *   x* = (4,4)
 *   f* = -304
 */
public class HS224Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 4;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS224Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            return 2.0*x1*x1 + x2*x2 - 48.0*x1 - 40.0*x2;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double g1 = 4.0*x1 - 48.0;
            double g2 = 2.0*x2 - 40.0;
            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Hessian = [[4,0], [0,2]]
            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);
            H.setEntry(0,0,4.0);
            H.setEntry(1,1,2.0);
            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Inequalities G1..G4
    // -------------------------------------------------------------------------
    private static class HS224Ineq extends InequalityConstraint {

        HS224Ineq() {
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

            double g1 = x1 + 3.0*x2;
            double g2 = 18.0 - x1 - 3.0*x2;
            double g3 = x1 + x2;
            double g4 = 8.0 - x1 - x2;

            return new ArrayRealVector(new double[]{g1, g2, g3, g4}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // Row 1: d/d[x1,x2] of (x1 + 3x2)
            J.setEntry(0,0,1.0);
            J.setEntry(0,1,3.0);

            // Row 2: d/d[x1,x2] of (18 - x1 - 3x2)
            J.setEntry(1,0,-1.0);
            J.setEntry(1,1,-3.0);

            // Row 3: d/d[x1,x2] of (x1 + x2)
            J.setEntry(2,0,1.0);
            J.setEntry(2,1,1.0);

            // Row 4: d/d[x1,x2] of (8 - x1 - x2)
            J.setEntry(3,0,-1.0);
            J.setEntry(3,1,-1.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS224_optimization() {

        double[] x0 = new double[]{0.1, 0.1};

        double[] lower = new double[]{0.0, 0.0};
        double[] upper = new double[]{6.0, 6.0};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS224Obj()),
                null,
                new HS224Ineq(),
                bounds
        );

        double f = sol.getValue();

        final double fExpected = -304.0;
        final double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
