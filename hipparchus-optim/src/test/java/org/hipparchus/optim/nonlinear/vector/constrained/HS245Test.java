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
 * HS245 (TP245)
 *
 * N    = 3
 * NILI = 0
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * Least-squares problem with 10 residuals F(i):
 *
 *   D_i = i
 *   Z_i = D_i / 10
 *
 *   F_i(x) = exp(-Z_i * x1)
 *            - exp(-Z_i * x2)
 *            - x3 * (exp(-Z_i) - exp(-D_i))
 *
 * Objective (MODE = 2):
 *   FX = sum_{i=1..10} F_i(x)^2
 *
 * Bounds (MODE = 1):
 *   0   <= x1 <= 12
 *   0   <= x2 <= 12
 *   0   <= x3 <= 20
 *
 * Reference solution:
 *   x* = (1, 10, 1)
 *   f* = 0
 */
public class HS245Test {

    private static final int DIM = 3;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS245Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double fx = 0.0;

            for (int i = 1; i <= 10; i++) {
                double Di = (double) i;
                double Zi = Di / 10.0;

                double term1 = FastMath.exp(-Di * x1 / 10.0);
                double term2 = FastMath.exp(-Di * x2 / 10.0);
                double term3 = x3 * (FastMath.exp(-Di / 10.0) - FastMath.exp(-Di));

                double Fi = term1 - term2 - term3;
                fx += Fi * Fi;
            }

            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double g1 = 0.0;
            double g2 = 0.0;
            double g3 = 0.0;

            for (int i = 1; i <= 10; i++) {
                double Di = (double) i;

                double exp1 = FastMath.exp(-Di * x1 / 10.0);
                double exp2 = FastMath.exp(-Di * x2 / 10.0);
                double expZ = FastMath.exp(-Di / 10.0);
                double expD = FastMath.exp(-Di);

                double Fi = exp1 - exp2 - x3 * (expZ - expD);

                double dFi_dx1 = -Di / 10.0 * exp1;
                double dFi_dx2 =  Di / 10.0 * exp2;
                double dFi_dx3 = -(expZ - expD);

                g1 += 2.0 * Fi * dFi_dx1;
                g2 += 2.0 * Fi * dFi_dx2;
                g3 += 2.0 * Fi * dFi_dx3;
            }

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Let the SQP/BFGS machinery approximate the Hessian.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS245_optimization() {

        // Initial guess (MODE=1):
        //   X(1) = 0, X(2) = 10, X(3) = 20
        double[] x0 = new double[]{0.0, 10.0, 20.0};

        // Bounds:
        //   0 <= x1 <= 12
        //   0 <= x2 <= 12
        //   0 <= x3 <= 20
        double[] lower = new double[]{0.0, 0.0, 0.0};
        double[] upper = new double[]{12.0, 12.0, 20.0};
         SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS245Obj()),
                null,   // no equalities
                null,   // no inequalities
                bounds  // bounds from Fortran XL, XU
        );

        double f = sol.getValue();

        final double fExpected = 0.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
