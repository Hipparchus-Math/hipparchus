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
 * HS260 (TP260)
 *
 * N     = 4
 * NILI  = 0
 * NINL  = 0
 * NELI  = 0
 * NENL  = 0
 *
 * Fortran TP260:
 *
 * F(1) = 10 (x2 - x1^2)
 * F(2) = 1 - x1
 * F(3) = sqrt(90) (x4 - x3^2)
 * F(4) = 1 - x3
 * F(5) = sqrt(9.9) ((x2 - 1) + (x4 - 1))
 * F(6) = sqrt(0.2) (x2 - 1)
 * F(7) = sqrt(0.2) (x4 - 1)
 *
 * FX   = sum_{i=1..7} F(i)^2
 *
 * Bounds: none (unconstrained).
 *
 * Initial guess:
 *   x0 = (-3, -1, -3, -1)
 *
 * Reference solution:
 *   x*  = (1, 1, 1, 1)
 *   f*  = 0
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

public class HS260Test {

    private static final int DIM = 4;

    // -------------------------------------------------------------------------
    // Objective (sum of squares)
    // -------------------------------------------------------------------------
    private static class HS260Obj extends TwiceDifferentiableFunction {

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

            double f1 = 10.0 * (x2 - x1 * x1);
            double f2 = 1.0 - x1;
            double f3 = FastMath.sqrt(90.0) * (x4 - x3 * x3);
            double f4 = 1.0 - x3;
            double f5 = FastMath.sqrt(9.9) * ((x2 - 1.0) + (x4 - 1.0));
            double f6 = FastMath.sqrt(0.2) * (x2 - 1.0);
            double f7 = FastMath.sqrt(0.2) * (x4 - 1.0);

            return f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4
                 + f5 * f5 + f6 * f6 + f7 * f7;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            double sqrt90  = FastMath.sqrt(90.0);
            double sqrt9_9 = FastMath.sqrt(9.9);
            double sqrt0_2 = FastMath.sqrt(0.2);

            double f1 = 10.0 * (x2 - x1 * x1);
            double f2 = 1.0 - x1;
            double f3 = sqrt90 * (x4 - x3 * x3);
            double f4 = 1.0 - x3;
            double f5 = sqrt9_9 * ((x2 - 1.0) + (x4 - 1.0));
            double f6 = sqrt0_2 * (x2 - 1.0);
            double f7 = sqrt0_2 * (x4 - 1.0);

            // Rebuild DF exactly as in Fortran after MODE=3 modifications:
            // Row j, column i: DF(j,i)
            double[][] DF = new double[7][4];

            // Initialize to 0, then set entries:

            // From MODE=1 constants:
            DF[0][1] = 10.0;              // DF(1,2)
            DF[1][0] = -1.0;              // DF(2,1)
            DF[2][3] = sqrt90;            // DF(3,4)
            DF[3][2] = -1.0;              // DF(4,3)
            DF[4][1] = sqrt9_9;           // DF(5,2)
            DF[4][3] = sqrt9_9;           // DF(5,4)
            DF[5][1] = sqrt0_2;           // DF(6,2)
            DF[6][3] = sqrt0_2;           // DF(7,4)

            // MODE=3 modifications:
            DF[0][0] = -20.0 * x1;        // DF(1,1)
            DF[2][2] = -sqrt90 * 2.0 * x3;// DF(3,3)

            double[] F = new double[]{f1, f2, f3, f4, f5, f6, f7};
            double[] g = new double[4];

            // GF(i) = sum_j 2 * F(j) * DF(j,i)
            for (int i = 0; i < 4; i++) {
                double gi = 0.0;
                for (int j = 0; j < 7; j++) {
                    gi += 2.0 * F[j] * DF[j][i];
                }
                g[i] = gi;
            }

            return new ArrayRealVector(g, false);
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
    public void testHS260() {

        double[] x0 = new double[]{-3.0, -1.0, -3.0, -1.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS260Obj()),
                null,   // no equalities
                null,   // no inequalities
                null    // no bounds
        );

        double f = sol.getValue();

        double fExpected = 0.0;
        double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
