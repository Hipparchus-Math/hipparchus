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
 * HS256 (TP256)
 *
 * N     = 4
 * NILI  = 0
 * NINL  = 0
 * NELI  = 0
 * NENL  = 0
 *
 * Objective (Fortran TP256):
 *
 *   f(x) =
 *       (x1 + 10*x2)^2
 *     + 5 * (x3 - x4)^2
 *     + (x2 - 2*x3)^4
 *     + 10 * (x1 - x4)^4
 *
 * Gradient (MODE=3 in Fortran):
 *
 *   g1 =  2*(x1 + 10*x2) + 40*(x1 - x4)^3
 *   g2 = 20*(x1 + 10*x2) +  4*(x2 - 2*x3)^3
 *   g3 = 10*(x3 - x4)    -  8*(x2 - 2*x3)^3
 *   g4 = -10*(x3 - x4)   - 40*(x1 - x4)^3
 *
 * Initial guess (MODE=1):
 *   x = (3, -1, 0, 1)
 *
 * Reference solution (MODE=1):
 *   x* = (0, 0, 0, 0)
 *   f* = 0
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

public class HS256Test {

    private static final int DIM = 4;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS256Obj extends TwiceDifferentiableFunction {

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

            double term1 = (x1 + 10.0 * x2);
            term1 = term1 * term1;

            double term2 = x3 - x4;
            term2 = 5.0 * term2 * term2;

            double term3 = x2 - 2.0 * x3;
            term3 = term3 * term3 * term3 * term3; // (x2 - 2*x3)^4

            double term4 = x1 - x4;
            term4 = 10.0 * term4 * term4 * term4 * term4; // 10*(x1 - x4)^4

            return term1 + term2 + term3 + term4;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            double t12 = x1 + 10.0 * x2;       // (x1 + 10*x2)
            double t34 = x3 - x4;              // (x3 - x4)
            double t23 = x2 - 2.0 * x3;        // (x2 - 2*x3)
            double t14 = x1 - x4;              // (x1 - x4)

            // From Fortran:
            // GF(1)=(2*(x1+10*x2)+40*(x1-x4)^3)
            // GF(2)=(20*(x1+10*x2)+4*(x2-2*x3)^3)
            // GF(3)=(10*(x3-x4)-8*(x2-2*x3)^3)
            // GF(4)=(-10*(x3-x4)-40*(x1-x4)^3)
            double g1 = 2.0 * t12 + 40.0 * t14 * t14 * t14;
            double g2 = 20.0 * t12 + 4.0 * t23 * t23 * t23;
            double g3 = 10.0 * t34 - 8.0 * t23 * t23 * t23;
            double g4 = -10.0 * t34 - 40.0 * t14 * t14 * t14;

            return new ArrayRealVector(new double[]{g1, g2, g3, g4}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // True Hessian is complicated; let BFGS approximate it.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS256_optimization() {

        // Initial guess (MODE=1): X = (3, -1, 0, 1)
        double[] x0 = new double[]{3.0, -1.0, 0.0, 1.0};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Unconstrained problem: no equalities, no inequalities, no bounds
        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS256Obj()),
                null,   // no equalities
                null,   // no inequalities
                null    // no bounds
        );

        double f = sol.getValue();

        // Reference optimum: x* = (0,0,0,0), f* = 0
        double fExpected = 0.0;
        double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
