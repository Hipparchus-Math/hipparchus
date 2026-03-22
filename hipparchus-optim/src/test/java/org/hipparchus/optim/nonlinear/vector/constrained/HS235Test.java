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
 * HS235 (TP235)
 *
 * N    = 3
 * NILI = 0
 * NINL = 0
 * NELI = 0
 * NENL = 1 (one nonlinear constraint, originally equality in Fortran)
 *
 * Objective (MODE=2):
 *   f(x) = (x2 - x1^2)^2 + 0.01 (x1 - 1)^2
 *
 * Fortran constraint (G(x) >= 0; originally equality):
 *   G1(x) = x1 + x3^2 + 1
 *
 * Reference solution (MODE=1):
 *   x*  = (-1, 1, 0)
 *   f*  = 0.04
 */
public class HS235Test {

    private static final int DIM      = 3;
    private static final int NUM_INEQ = 1; // model the Fortran constraint as G >= 0
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS235Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double t  = x2 - x1 * x1;
            return t * t + 0.01 * (x1 - 1.0) * (x1 - 1.0);
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            double t = x2 - x1 * x1;

            // From Fortran:
            // GF(1) = -4*x1*(x2 - x1^2) + 0.02*(x1 - 1)
            // GF(2) =  2*(x2 - x1^2)
            // GF(3) =  0
            double df1 = -4.0 * x1 * t + 0.02 * (x1 - 1.0);
            double df2 =  2.0 * t;
            double df3 =  0.0;

            return new ArrayRealVector(new double[]{df1, df2, df3}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Zero Hessian; let BFGS/SQP handle curvature
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Nonlinear constraint G(x) = 0 (originally equality in Fortran)
    // -------------------------------------------------------------------------
    private static class HS235eq extends EqualityConstraint {

        HS235eq() {
            // One constraint, RHS = 0
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x3 = x.getEntry(2);

            // Fortran: G(1) = x1 + x3^2 + 1  (used as equality G = 0)
            // Here we keep G(x) as-is with convention G >= 0.
            double G1 = x1 + x3 * x3 + 1.0;

            return new ArrayRealVector(new double[]{G1}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            double x1 = x.getEntry(0); // unused, but kept for clarity
            double x3 = x.getEntry(2);

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // G1 = x1 + x3^2 + 1
            // dG1/dx1 = 1
            // dG1/dx2 = 0
            // dG1/dx3 = 2*x3
            J.setEntry(0, 0, 1.0);
            J.setEntry(0, 1, 0.0);
            J.setEntry(0, 2, 2.0 * x3);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS235_optimization() {

        // Initial guess (MODE=1): X(1)=-2, X(2)=3, X(3)=1
        double[] x0 = new double[]{-2.0, 3.0, 1.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS235Obj()),
                null,               // no explicit equalities (we model Fortran equality as inequality)
                new HS235eq(),    // 1 nonlinear constraint
                null                // no bounds
        );

        double f = sol.getValue();

        final double fExpected = 0.04;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
