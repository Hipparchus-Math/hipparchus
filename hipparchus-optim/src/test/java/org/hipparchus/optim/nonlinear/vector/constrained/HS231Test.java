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
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/*
 * HS231 (TP231)
 *
 * N    = 2
 * NILI = 2 (two linear inequalities)
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * Objective:
 *   f(x) = 100 (x2 - x1^2)^2 + (1 - x1)^2
 *
 * Fortran constraints (G(x) ≥ 0):
 *   G1(x) =  x1/3 + x2 + 0.1
 *   G2(x) = -x1/3 + x2 + 0.1
 *
 * Reference solution (MODE=1):
 *   x*  = (1, 1)
 *   f*  = 0
 */
public class HS231Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 2;
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS231Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double t  = x2 - x1 * x1;
            return 100.0 * t * t + (1.0 - x1) * (1.0 - x1);
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double t  = x2 - x1 * x1;

            double df1 = -400.0 * x1 * t - 2.0 * (1.0 - x1);
            double df2 =  200.0 * t;

            return new ArrayRealVector(new double[]{df1, df2}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints G(x) ≥ 0
    // -------------------------------------------------------------------------
    private static class HS231Ineq extends InequalityConstraint {

        HS231Ineq() {
            // RHS = 0 for both inequalities
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

            double G1 =  x1 / 3.0 + x2 + 0.1;
            double G2 = -x1 / 3.0 + x2 + 0.1;

            return new ArrayRealVector(new double[]{G1, G2}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // G1 = x1/3 + x2 + 0.1 → grad = [1/3, 1]
            J.setEntry(0, 0,  1.0 / 3.0);
            J.setEntry(0, 1,  1.0);

            // G2 = -x1/3 + x2 + 0.1 → grad = [-1/3, 1]
            J.setEntry(1, 0, -1.0 / 3.0);
            J.setEntry(1, 1,  1.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS231_optimization() {

        // Initial guess (MODE=1): X(1) = -1.2, X(2) = 1
        double[] x0 = new double[]{-1.2, 1.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS231Obj()),
                null,                // no equalities
                new HS231Ineq(),     // 2 linear inequalities
                null                 // no bounds
        );

        double f = sol.getValue();

        final double fExpected = 0.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
