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
 * HS212 (TP212) – 2D nonconvex problem.
 *
 * Fortran TP212:
 *
 *   N    = 2
 *   NILI = 0
 *   NINL = 0
 *   NELI = 0
 *   NENL = 0
 *
 *   Initial guess:
 *     x1 =  2.0
 *     x2 =  0.0
 *
 *   Objective:
 *
 *     Let
 *       S  = x1 + x2
 *       D  = x1 - x2
 *       R  = (x1 - 2)^2 + x2^2 - 1
 *       W  = 4*S + D*R
 *
 *     f(x) = (4*S)^2 + W^2
 *
 *   Expanded exactly as in Fortran:
 *
 *     FX = (4*(x1 + x2))^2
 *        + (4*(x1 + x2) + (x1 - x2)*((x1 - 2)^2 + x2^2 - 1))^2
 *
 *   Gradient (from Fortran):
 *
 *     GF(1) = 32*(x1 + x2)
 *             + 2*W*(4 + ((x1 - 2)^2 + x2^2 - 1)
 *                    + (x1 - x2)*2*(x1 - 2))
 *
 *     GF(2) = 32*(x1 + x2)
 *             + 2*W*(4 - (x1 - 2)^2 + x2^2 - 1
 *                    + (x1 - x2)*2*x2)
 *
 *   No constraints, no bounds.
 *
 *   LEX = .TRUE., FEX = 0, XEX = (0, 0)
 */
public class HS212Test {

    /** Dimension. */
    private static final int DIM = 2;

    // -------------------------------------------------------------------------
    // Objective f(x)
    // -------------------------------------------------------------------------

    private static class HS212Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double S  = x1 + x2;
            final double D  = x1 - x2;
            final double R  = (x1 - 2.0) * (x1 - 2.0) + x2 * x2 - 1.0;
            final double W  = 4.0 * S + D * R;

            final double term1 = 4.0 * S;
            final double fx = term1 * term1 + W * W;

            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double S  = x1 + x2;
            final double D  = x1 - x2;
            final double R  = (x1 - 2.0) * (x1 - 2.0) + x2 * x2 - 1.0;
            final double W  = 4.0 * S + D * R;

            // Direct transcription from Fortran:

            // GF(1)=32*(x1+x2)
            //        +2*(4*(x1+x2)+(x1-x2)*((x1-2)^2+x2^2-1))
            //           * (4 + ((x1-2)^2 + x2^2 -1)
            //                + (x1-x2)*2*(x1-2))

            final double df1 = 32.0 * (x1 + x2)
                    + 2.0 * W
                      * (4.0
                         + ((x1 - 2.0) * (x1 - 2.0) + x2 * x2 - 1.0)
                         + (x1 - x2) * 2.0 * (x1 - 2.0));

            // GF(2)=32*(x1+x2)
            //        +2*(4*(x1+x2)+(x1-x2)*((x1-2)^2+x2^2-1))
            //           * (4 - (x1-2)^2 + x2^2 -1
            //                + (x1-x2)*2*x2)

            final double df2 = 32.0 * (x1 + x2)
                    + 2.0 * W
                      * (4.0
                         - (x1 - 2.0) * (x1 - 2.0)
                         + x2 * x2
                         - 1.0
                         + (x1 - x2) * 2.0 * x2);

            return new ArrayRealVector(new double[] { df1, df2 }, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Start with zero Hessian; let the SQP BFGS updater build curvature.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test driver using SQPOptimizerS2
    // -------------------------------------------------------------------------

    @Test
    public void testHS212_optimization() {

        // Initial guess from Fortran (mode 1):
        final double[] x0 = new double[DIM];
        x0[0] = 2.0;   // X(1)
        x0[1] = 0.0;   // X(2)

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new HS212Obj()),
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
