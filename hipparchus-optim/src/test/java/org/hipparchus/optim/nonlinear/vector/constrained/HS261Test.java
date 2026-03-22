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

/*
 * HS261 (TP261)
 *
 * N     = 4
 * NILI  = 0
 * NINL  = 0
 * NELI  = 0
 * NENL  = 0
 *
 * Fortran TP261:
 *
 * F(1) = (exp(x1) - x2)^2
 * F(2) = 10 (x2 - x3)^3
 * F(3) = tan(x3 - x4)^2
 * F(4) = x1^4
 * F(5) = x4 - 1
 *
 * FX   = sum_{i=1..5} F(i)^2
 *
 * A = exp(x1) - x2
 * B = tan(x3 - x4)
 * C = B / cos(x3 - x4)^2
 *
 * Gradient (MODE=3 in Fortran):
 *   GF(1) = 4 * exp(x1) * A^3 + 8 * x1^7
 *   GF(2) = -4 * A^3 + 600 (x2 - x3)^5
 *   GF(3) = 4 * B^2 * C - 600 (x2 - x3)^5
 *   GF(4) = -4 * B^2 * C + 2 (x4 - 1)
 *
 * Bounds (MODE=1):
 *   0 ≤ xi ≤ 10, i = 1..4
 *
 * Initial guess:
 *   x0 = (0, 0, 0, 0)
 *
 * Reference solution:
 *   x*  = (0, 1, 1, 1)
 *   f*  = 0
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

public class HS261Test {

    private static final int DIM = 4;

    // -------------------------------------------------------------------------
    // Objective (sum of squares)
    // -------------------------------------------------------------------------
    private static class HS261Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            double A = FastMath.exp(x1) - x2;
            double B = FastMath.tan(x3 - x4);

            double f1 = A * A;
            double f2 = 10.0 * FastMath.pow(x2 - x3, 3);
            double f3 = B * B;
            double f4 = FastMath.pow(x1, 4);
            double f5 = x4 - 1.0;

            return f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4 + f5 * f5;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            double A = FastMath.exp(x1) - x2;
            double B = FastMath.tan(x3 - x4);
            double cosArg = FastMath.cos(x3 - x4);
            double C = B / (cosArg * cosArg);

            double g1 = 4.0 * FastMath.exp(x1) * FastMath.pow(A, 3)
                      + 8.0 * FastMath.pow(x1, 7);

            double g2 = -4.0 * FastMath.pow(A, 3)
                      + 600.0 * FastMath.pow(x2 - x3, 5);

            double g3 = 4.0 * B * B * C
                      - 600.0 * FastMath.pow(x2 - x3, 5);

            double g4 = -4.0 * B * B * C
                      + 2.0 * (x4 - 1.0);

            return new ArrayRealVector(new double[]{g1, g2, g3, g4}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Let BFGS approximate the Hessian
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS261() {

        double[] x0 = new double[]{0.0, 0.0, 0.0, 0.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Box bounds: 0 ≤ xi ≤ 10
        SimpleBounds bounds = new SimpleBounds(
                new double[]{0.0, 0.0, 0.0, 0.0},
                new double[]{10.0, 10.0, 10.0, 10.0}
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS261Obj()),
                null,   // no equalities
                null,   // no inequalities
                bounds
        );

        double f = sol.getValue();

        double fExpected = 0.0;
        double tol = 1e-5 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
