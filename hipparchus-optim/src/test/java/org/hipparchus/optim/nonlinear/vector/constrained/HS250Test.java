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

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/*
 * HS250 (TP250)
 *
 * NILI = 2  (linear inequalities)
 * Objective: f(x) = -x1*x2*x3
 *
 * Fortran constraints G(i) >= 0. Wrapper uses g <= 0, so we negate them.
 *
 * G1 = x1 + 2 x2 + 2 x3
 *     wrapper g1 = -(x1 + 2 x2 + 2 x3)
 *
 * G2 = 72 - x1 - 2 x2 - 2 x3
 *     wrapper g2 = x1 + 2 x2 + 2 x3 - 72
 *
 * Bounds:
 *    0 ≤ x1 ≤ 20
 *    0 ≤ x2 ≤ 11
 *    0 ≤ x3 ≤ 42
 *
 * Reference optimum:
 *    x* = (20, 11, 15)
 *    f* = -3300
 */
public class HS250Test {

    private static final int DIM      = 3;
    private static final int NUM_INEQ = 2;

    // ---------------------------------------------------------
    // Objective
    // ---------------------------------------------------------
    private static class HS250Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            return -x.getEntry(0) * x.getEntry(1) * x.getEntry(2);
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            return new ArrayRealVector(new double[]{
                    -x2 * x3,   // df/dx1
                    -x1 * x3,   // df/dx2
                    -x1 * x2    // df/dx3
            }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);

            // Second derivatives of -x1*x2*x3
            H.setEntry(0, 1, -x3);
            H.setEntry(1, 0, -x3);

            H.setEntry(0, 2, -x2);
            H.setEntry(2, 0, -x2);

            H.setEntry(1, 2, -x1);
            H.setEntry(2, 1, -x1);

            return H;
        }
    }

    // ---------------------------------------------------------
    // Inequalities g(x) <= 0
    // ---------------------------------------------------------
    private static class HS250Ineq extends InequalityConstraint {

        HS250Ineq() {
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
            double x3 = x.getEntry(2);

            double g1 = -(x1 + 2 * x2 + 2 * x3);

            double g2 = x1 + 2 * x2 + 2 * x3 - 72;

            return new ArrayRealVector(new double[]{-g1, -g2}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // g1 = -(x1 + 2x2 + 2x3)
            J.setEntry(0, 0, -1);
            J.setEntry(0, 1, -2);
            J.setEntry(0, 2, -2);

            // g2 = x1 + 2x2 + 2x3 - 72
            J.setEntry(1, 0, 1);
            J.setEntry(1, 1, 2);
            J.setEntry(1, 2, 2);

            return J;
        }
    }

    // ---------------------------------------------------------
    // Bounds

// Box Constraints (Variable Bounds)
SimpleBounds bounds = new SimpleBounds(
        new double[] {
                0.0,  // XL(1)
                0.0,  // XL(2)
                0.0   // XL(3)
        },
        new double[] {
                20.0, // XU(1)
                11.0, // XU(2)
                42.0  // XU(3)
        }
);

    // ---------------------------------------------------------
    // JUnit Test
    // ---------------------------------------------------------
    @Test
    public void testHS250_optimization() {

        double[] x0 = new double[]{10.0, 10.0, 10.0};  // Fortran initial guess = 10,10,10

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS250Obj()),
                null,
                new HS250Ineq(),
               bounds
        );

        double f = sol.getValue();

        final double fExpected = -3300.0;
        final double tol = 1e-6 * (FastMath.abs(fExpected) + 1);

        assertEquals(fExpected, f, tol);
    }
}
