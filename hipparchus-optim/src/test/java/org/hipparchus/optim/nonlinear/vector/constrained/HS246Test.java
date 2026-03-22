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
 * HS246 (TP246)
 *
 * N    = 3
 * NILI = 0
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * Least-squares problem with 3 residuals F(j):
 *
 *   F1(x) = 10 * ( x3 - ((x1 + x2)/2)^2 )
 *   F2(x) = 1 - x1
 *   F3(x) = 1 - x2
 *
 * Objective (MODE = 2):
 *   FX = F1^2 + F2^2 + F3^2
 *
 * Fortran derivatives (MODE = 3):
 *   DF(1,1) = -10 * (x1 + x2)
 *   DF(1,2) = -10 * (x1 + x2)
 *   DF(1,3) =  10
 *   DF(2,1) = -1
 *   DF(3,2) = -1
 *   others  =  0
 *
 * Gradient in Fortran:
 *   GF(i) = sum_j 2 * F(j) * DF(j,i)
 *
 * Reference solution (MODE = 1):
 *   x*  = (1, 1, 1)
 *   f*  = 0
 */
public class HS246Test {

    private static final int DIM = 3;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS246Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double avg = 0.5 * (x1 + x2);

            double F1 = 10.0 * (x3 - avg * avg);
            double F2 = 1.0 - x1;
            double F3 = 1.0 - x2;

            return F1 * F1 + F2 * F2 + F3 * F3;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double avg = 0.5 * (x1 + x2);

            // Residuals (as in MODE=2)
            double F1 = 10.0 * (x3 - avg * avg);
            double F2 = 1.0 - x1;
            double F3 = 1.0 - x2;

            // Derivatives DF(j,i) exactly as in Fortran MODE=3
            double DF11 = -10.0 * (x1 + x2); // dF1/dx1
            double DF12 = -10.0 * (x1 + x2); // dF1/dx2
            double DF13 =  10.0;             // dF1/dx3

            double DF21 = -1.0;              // dF2/dx1
            double DF22 =  0.0;              // dF2/dx2
            double DF23 =  0.0;              // dF2/dx3

            double DF31 =  0.0;              // dF3/dx1
            double DF32 = -1.0;              // dF3/dx2
            double DF33 =  0.0;              // dF3/dx3

            // GF(i) = sum_j 2 * F(j) * DF(j,i)
            double g1 = 2.0 * (F1 * DF11 + F2 * DF21 + F3 * DF31);
            double g2 = 2.0 * (F1 * DF12 + F2 * DF22 + F3 * DF32);
            double g3 = 2.0 * (F1 * DF13 + F2 * DF23 + F3 * DF33);

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
    public void testHS246_optimization() {

        // Initial guess (MODE=1):
        //   X(1) = -1.2, X(2) = 2, X(3) = 0
        double[] x0 = new double[]{-1.2, 2.0, 0.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // No constraints, no bounds
        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS246Obj()),
                null,   // no equalities
                null,   // no inequalities
                null    // no bounds
        );

        double f = sol.getValue();

        final double fExpected = 0.0; // FEX from Fortran
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
