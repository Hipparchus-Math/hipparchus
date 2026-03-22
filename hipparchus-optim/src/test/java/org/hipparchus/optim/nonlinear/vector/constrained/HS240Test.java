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
 * HS240 (TP240)
 *
 * N    = 3
 * NILI = 0  (no linear inequalities)
 * NINL = 0  (no nonlinear inequalities)
 * NELI = 0  (no linear equalities)
 * NENL = 0  (no nonlinear equalities)
 *
 * Objective (MODE=2):
 *   f(x) = (x1 - x2 + x3)^2
 *        + (-x1 + x2 + x3)^2
 *        + (x1 + x2 - x3)^2
 *
 * Gradient (MODE=3 in Fortran, simplified):
 *   df/dx1 =  6*x1 - 2*x2 - 2*x3
 *   df/dx2 = -2*x1 + 6*x2 - 2*x3
 *   df/dx3 = -2*x1 - 2*x2 + 6*x3
 *
 * Hessian (constant):
 *   H = [  6  -2  -2 ]
 *       [ -2   6  -2 ]
 *       [ -2  -2   6 ]
 *
 * Reference solution (MODE=1):
 *   x*  = (0, 0, 0)
 *   f*  = 0
 */
public class HS240Test {

    private static final int DIM      = 3;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS240Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double t1 =  x1 - x2 + x3;
            double t2 = -x1 + x2 + x3;
            double t3 =  x1 + x2 - x3;

            return t1 * t1 + t2 * t2 + t3 * t3;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            double g1 =  6.0 * x1 - 2.0 * x2 - 2.0 * x3;
            double g2 = -2.0 * x1 + 6.0 * x2 - 2.0 * x3;
            double g3 = -2.0 * x1 - 2.0 * x2 + 6.0 * x3;

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Constant Hessian
            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);
            H.setEntry(0, 0,  6.0);
            H.setEntry(0, 1, -2.0);
            H.setEntry(0, 2, -2.0);

            H.setEntry(1, 0, -2.0);
            H.setEntry(1, 1,  6.0);
            H.setEntry(1, 2, -2.0);

            H.setEntry(2, 0, -2.0);
            H.setEntry(2, 1, -2.0);
            H.setEntry(2, 2,  6.0);

            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Test (unconstrained)
    // -------------------------------------------------------------------------
    @Test
    public void testHS240_optimization() {

        // Initial guess from Fortran MODE=1:
        //   X(1)=1.0E+2, X(2)=-1.0, X(3)=2.5
        double[] x0 = new double[]{100.0, -1.0, 2.5};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS240Obj()),
                null,  // no equalities
                null,  // no inequalities
                null   // no bounds
        );

        double f = sol.getValue();

        // Reference optimum from Fortran: FEX = 0.0
        final double fExpected = 0.0;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
