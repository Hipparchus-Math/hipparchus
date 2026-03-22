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
 * HS251 (TP251)
 *
 * N    = 3
 * NILI = 1  (one linear inequality)
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * Objective (MODE=2):
 *   f(x) = -x1 * x2 * x3
 *
 * Fortran inequality (G(i) ≥ 0):
 *   G1(x) = 72 - x1 - 2 x2 - 2 x3 ≥ 0
 *
 * We keep the same ≥ 0 convention used by the original Fortran
 * and by the SQP optimizer: the constraint function returned by
 * {@link InequalityConstraint#value(RealVector)} is G1 itself.
 *
 * Reference solution (MODE=1):
 *   x*  = (20, 11, 15)
 *   f*  = -3456
 * with bounds:
 *   0 ≤ x1 ≤ 42,  0 ≤ x2 ≤ 42,  0 ≤ x3 ≤ 42
 * but with tighter upper bounds from TP251:
 *   0 ≤ x1 ≤ 20
 *   0 ≤ x2 ≤ 11
 *   0 ≤ x3 ≤ 42
 */
public class HS251Test {

    private static final int DIM      = 3;
    private static final int NUM_INEQ = 1;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS251Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            return -x1 * x2 * x3;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            return new ArrayRealVector(new double[] {
                    -x2 * x3, // df/dx1
                    -x1 * x3, // df/dx2
                    -x1 * x2  // df/dx3
            }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            // f(x) = -x1 x2 x3
            // ∂²f/∂x1²  = 0
            // ∂²f/∂x2²  = 0
            // ∂²f/∂x3²  = 0
            // ∂²f/∂x1∂x2 = -x3
            // ∂²f/∂x1∂x3 = -x2
            // ∂²f/∂x2∂x3 = -x1
            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);

            H.setEntry(0, 0, 0.0);
            H.setEntry(1, 1, 0.0);
            H.setEntry(2, 2, 0.0);

            H.setEntry(0, 1, -x3);
            H.setEntry(1, 0, -x3);

            H.setEntry(0, 2, -x2);
            H.setEntry(2, 0, -x2);

            H.setEntry(1, 2, -x1);
            H.setEntry(2, 1, -x1);

            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints (G >= 0, Fortran convention)
    // -------------------------------------------------------------------------
    private static class HS251Ineq extends InequalityConstraint {

        HS251Ineq() {
            // RHS = 0 for all inequalities (G(x) >= 0)
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

            // Fortran:
            // G(1) = 72 - x1 - 2*x2 - 2*x3  (>= 0)
            double g1 = 72.0 - x1 - 2.0 * x2 - 2.0 * x3;

            // SQP optimizer convention: constraints are returned as G(x) >= 0
            return new ArrayRealVector(new double[]{g1}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // g1 = 72 - x1 - 2*x2 - 2*x3 → grad = [-1, -2, -2]
            J.setEntry(0, 0, -1.0);
            J.setEntry(0, 1, -2.0);
            J.setEntry(0, 2, -2.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS251_optimization() {

        // Initial guess (MODE=1): X(i) = 10
        double[] x0 = new double[]{10.0, 10.0, 10.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Box Constraints (Variable Bounds) from TP251:
        // 0 <= x1 <= 42
        // 0 <= x2 <= 42
        // 0 <= x3 <= 42
        SimpleBounds bounds = new SimpleBounds(
                new double[]{
                        0.0, // XL(1)
                        0.0, // XL(2)
                        0.0  // XL(3)
                },
                new double[]{
                        42.0, // XU(1)
                        42.0, // XU(2)
                        42.0  // XU(3)
                }
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS251Obj()),
                null,            // no equalities
                new HS251Ineq(), // 1 inequality
                bounds
        );

        double f = sol.getValue();

        // Reference optimum from Fortran: FEX = -3456
        final double fExpected = -3456.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
