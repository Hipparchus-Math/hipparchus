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

package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * HS211 (TP211) – Modified Rosenbrock-type problem.
 *
 * Fortran TP211:
 *
 *   N    = 2
 *   NILI = 0
 *   NINL = 0
 *   NELI = 0
 *   NENL = 0
 *
 *   Initial guess:
 *     x1 = -1.2
 *     x2 =  1.0
 *
 *   Objective:
 *
 *     f(x) = 100 * (x2 - x1^3)^2 + (1 - x1)^2
 *
 *   Gradient (from Fortran):
 *
 *     df/dx1 = -200*(x2 - x1^3)*3*x1^2 - 2*(1 - x1)
 *     df/dx2 =  200*(x2 - x1^3)
 *
 *   No constraints, no bounds.
 *
 *   LEX = .TRUE., FEX = 0, XEX = (1, 1)
 */
public class HS211Test {

    /** Problem dimension. */
    private static final int DIM = 2;

    // -------------------------------------------------------------------------
    // Objective: f(x) = 100*(x2 - x1^3)^2 + (1 - x1)^2
    // -------------------------------------------------------------------------

    private static class HS211Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double t = x2 - x1 * x1 * x1; // x2 - x1^3

            return 100.0 * t * t + (1.0 - x1) * (1.0 - x1);
        }

        @Override
        public RealVector gradient(RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double t = x2 - x1 * x1 * x1; // x2 - x1^3

            // From Fortran:
            // GF(1) = -200*(x2 - x1^3)*3*x1^2 - 2*(1 - x1)
            // GF(2) =  200*(x2 - x1^3)

            final double df1 = -200.0 * t * 3.0 * x1 * x1 - 2.0 * (1.0 - x1);
            final double df2 =  200.0 * t;

            return new ArrayRealVector(new double[] { df1, df2 }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Start SQP with zero Hessian; BFGS (or other updater) will build it.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test driver using SQPOptimizerS2
    // -------------------------------------------------------------------------

    @Test
    public void testHS211_optimization() {

        // Initial guess from Fortran (mode 1):
        final double[] x0 = new double[DIM];
        x0[0] = -1.2;  // X(1)
        x0[1] =  1.0;  // X(2)

        // No bounds, no constraints for this problem.
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new HS211Obj()),
            null,   // no equality constraints
            null,   // no inequality constraints
            null    // no bounds
        );

        final double f = sol.getValue();

        // Fortran: LEX = .TRUE., FEX = 0.0
        final double fExpected = 0.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
