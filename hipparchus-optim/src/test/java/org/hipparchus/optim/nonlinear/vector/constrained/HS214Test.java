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
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * HS214 (TP214) – 2D non-linear scalar test problem.
 *
 * Fortran TP214:
 *
 *   N    = 2
 *   NILI = 0
 *   NINL = 0
 *   NELI = 0
 *   NENL = 0
 *
 *   Initial guess:
 *     X(1) = -1.2
 *     X(2) =  1.0
 *
 *   Objective:
 *
 *     A(x) = 10 * (x1 - x2)^2 + (x1 - 1)^2
 *     f(x) = A(x)^(1/4)
 *
 *   In Fortran:
 *
 *     FX = (10.0D0*(X(1)-X(2))**2 + (X(1)-1.0D0)**2)**0.25D0
 *
 *   Gradient (from Fortran):
 *
 *     GF(1) = (0.25 / A^0.75) * (22 * x1 - 20 * x2 - 2)
 *     GF(2) = (0.25 / A^0.75) * (20 * (x2 - x1))
 *
 *   No constraints, no bounds.
 *
 *   LEX = .TRUE., FEX = 0, XEX = (1, 1)
 */
public class HS214Test {

    /** Dimension. */
    private static final int DIM = 2;

    // -------------------------------------------------------------------------
    // Objective f(x)
    // -------------------------------------------------------------------------

    private static class HS214Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double dx = x1 - x2;
            final double t  = x1 - 1.0;

            // A = 10*(x1 - x2)^2 + (x1 - 1)^2
            final double A = 10.0 * dx * dx + t * t;

            // f = A^0.25
            return FastMath.pow(A, 0.25);
        }

        @Override
        public RealVector gradient(RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double dx = x1 - x2;
            final double t  = x1 - 1.0;

            // A = 10*(x1 - x2)^2 + (x1 - 1)^2
            final double A = 10.0 * dx * dx + t * t;

            // Direct transcription of Fortran GF(1), GF(2):
            //
            // GF(1)=((0.25D+0/(10.0D0*(X(1)-X(2))**2+(X(1)-1.0D0)**2)**0.75D0)
            //         *(22.0D0*X(1)-20.D0*X(2)-2.0D0))
            //
            // GF(2)=((0.25D+0/(10.0D0*(X(1)-X(2))**2+(X(1)-1.0D0)**2)**0.75D0)
            //         *20.0D0*(X(2)-X(1)))

            final double denom = 0.25 / FastMath.pow(A, 0.75);

            final double df1 = denom * (22.0 * x1 - 20.0 * x2 - 2.0);
            final double df2 = denom * (20.0 * (x2 - x1));

            return new ArrayRealVector(new double[] { df1, df2 }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Start with zero Hessian; BFGS will build curvature information.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test driver using SQPOptimizerS2
    // -------------------------------------------------------------------------

    @Test
    public void testHS214_optimization() {

        // Initial guess from Fortran (mode 1):
        final double[] x0 = new double[DIM];
        x0[0] = -1.2;   // X(1)
        x0[1] =  1.0;   // X(2)

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new HS214Obj()),
            null,   // no equality constraints
            null,   // no inequality constraints
            null    // no bounds
        );

        final double f = sol.getValue();

        // Fortran: LEX = .TRUE., FEX = 0.0 at x* = (1,1)
        final double fExpected = 0.0;
        final double tol = 1.0e-3 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
