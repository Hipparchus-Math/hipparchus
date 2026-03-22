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
 * HS241 (TP241)
 *
 * N    = 3
 * NILI = 0  (no linear inequalities)
 * NINL = 0  (no nonlinear inequalities)
 * NELI = 0  (no linear equalities)
 * NENL = 0  (no nonlinear equalities)
 *
 * Fortran defines residuals F(1..5):
 *
 *   F1 = x1^2 + x2^2 + x3^2 - 1
 *   F2 = x1^2 + x2^2 + (x3 - 2)^2 - 1
 *   F3 = x1 + x2 + x3 - 1
 *   F4 = x1 + x2 - x3 + 1
 *   F5 = x1^3 + 3*x2^2 + (5*x3 - x1 + 1)^2 - 36
 *
 * Objective (MODE = 2):
 *   f(x) = sum_{i=1}^5 F_i(x)^2
 *
 * Gradient (MODE = 3):
 *   g = 2 * J^T * F, where J(i, :) = ∂F_i/∂x
 *
 * Reference solution (MODE = 1):
 *   x*  = (0, 0, 1)
 *   f*  = 0
 */
public class HS241Test {

    private static final int DIM      = 3;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS241Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            // Residuals F(1..5)
            double f1 = x1 * x1 + x2 * x2 + x3 * x3 - 1.0;
            double f2 = x1 * x1 + x2 * x2 + (x3 - 2.0) * (x3 - 2.0) - 1.0;
            double f3 = x1 + x2 + x3 - 1.0;
            double f4 = x1 + x2 - x3 + 1.0;
            double tmp = 5.0 * x3 - x1 + 1.0;
            double f5 = x1 * x1 * x1 + 3.0 * x2 * x2 + tmp * tmp - 36.0;

            return f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4 + f5 * f5;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            // Residuals F(1..5) (same as in value)
            double f1 = x1 * x1 + x2 * x2 + x3 * x3 - 1.0;
            double f2 = x1 * x1 + x2 * x2 + (x3 - 2.0) * (x3 - 2.0) - 1.0;
            double f3 = x1 + x2 + x3 - 1.0;
            double f4 = x1 + x2 - x3 + 1.0;
            double tmp = 5.0 * x3 - x1 + 1.0;
            double f5 = x1 * x1 * x1 + 3.0 * x2 * x2 + tmp * tmp - 36.0;

            // Jacobian rows DF(i, j) = dF_i / dx_j, following Fortran
            // F1: [2*x1, 2*x2, 2*x3]
            double df1x1 = 2.0 * x1;
            double df1x2 = 2.0 * x2;
            double df1x3 = 2.0 * x3;

            // F2: [2*x1, 2*x2, 2*(x3-2)]
            double df2x1 = df1x1;
            double df2x2 = df1x2;
            double df2x3 = 2.0 * (x3 - 2.0);

            // F3: [1, 1, 1]
            double df3x1 = 1.0;
            double df3x2 = 1.0;
            double df3x3 = 1.0;

            // F4: [1, 1, -1]
            double df4x1 = 1.0;
            double df4x2 = 1.0;
            double df4x3 = -1.0;

            // F5:
            // DF(5,1) = 3*x1^2 - 2*(5*x3 - x1 + 1)
            // DF(5,2) = 6*x2
            // DF(5,3) = 10*(5*x3 - x1 + 1)
            double df5x1 = 3.0 * x1 * x1 - 2.0 * (5.0 * x3 - x1 + 1.0);
            double df5x2 = 6.0 * x2;
            double df5x3 = 10.0 * (5.0 * x3 - x1 + 1.0);

            // Gradient g = 2 * Σ_i F_i * ∇F_i
            double g1 = 0.0;
            double g2 = 0.0;
            double g3 = 0.0;

            g1 += 2.0 * f1 * df1x1;
            g2 += 2.0 * f1 * df1x2;
            g3 += 2.0 * f1 * df1x3;

            g1 += 2.0 * f2 * df2x1;
            g2 += 2.0 * f2 * df2x2;
            g3 += 2.0 * f2 * df2x3;

            g1 += 2.0 * f3 * df3x1;
            g2 += 2.0 * f3 * df3x2;
            g3 += 2.0 * f3 * df3x3;

            g1 += 2.0 * f4 * df4x1;
            g2 += 2.0 * f4 * df4x2;
            g3 += 2.0 * f4 * df4x3;

            g1 += 2.0 * f5 * df5x1;
            g2 += 2.0 * f5 * df5x2;
            g3 += 2.0 * f5 * df5x3;

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Let SQP/BFGS build an approximation; start from zero matrix.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test (unconstrained)
    // -------------------------------------------------------------------------
    @Test
    public void testHS241_optimization() {

        // Initial guess (MODE=1): X = (1, 2, 0)
        double[] x0 = new double[]{1.0, 2.0, 0.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS241Obj()),
                null,  // no equalities
                null,  // no inequalities
                null   // no bounds
        );

        double f = sol.getValue();

        // Reference optimum from Fortran: FEX = 0.0 at x* = (0,0,1)
        final double fExpected = 0.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
