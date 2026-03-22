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
 * HS255 (TP255)
 *
 * N     = 4
 * NILI  = 0
 * NINL  = 0
 * NELI  = 0
 * NENL  = 0
 *
 * Fortran TP255:
 *
 *   MODE=2:
 *     FX = 100*(x2 - x1^2) + (1 - x1)^2
 *        + 90*(x4 - x3^2) + (1 - x3)^2
 *        + 10.1*((x2 - 1)^2 + (x4 - 1)^2)
 *        + 19.8*(x2 - 1)*(x4 - 1)
 *     FX = 0.5 * FX**2
 *
 *   MODE=3:
 *     FX = same inner expression F(x) (without 0.5 * ...^2)
 *     GF(1) = FX * (-198*x1 - 2)
 *     GF(2) = FX * (20.2*x2 + 19.8*x4 + 60)
 *     GF(3) = FX * (-178*x3 - 2)
 *     GF(4) = FX * (19.8*x2 + 20.2*x4 + 50)
 *
 * Hence in our notation:
 *
 *   f(x) = 0.5 * F(x)^2
 *   ∇f(x) = F(x) * ∇F(x)
 *
 * Bounds:
 *   -10 ≤ xi ≤ 10, i = 1..4
 *
 * Initial guess:
 *   x0 = (-3,  1, -3, 1)
 *
 * Reference solution:
 *   x* = (1, 1, 1, 1), f* = 0
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

public class HS255Test {

    private static final int DIM = 4;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS255Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        /** Inner Fortran expression F(x) (before squaring). */
        private double computeInnerF(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            double term1 = 100.0 * (x2 - x1 * x1);          // 100*(x2 - x1^2)
            double term2 = (1.0 - x1) * (1.0 - x1);         // (1 - x1)^2
            double term3 = 90.0 * (x4 - x3 * x3);           // 90*(x4 - x3^2)
            double term4 = (1.0 - x3) * (1.0 - x3);         // (1 - x3)^2

            double t2m1 = x2 - 1.0;
            double t4m1 = x4 - 1.0;

            double term5 = 10.1 * (t2m1 * t2m1 + t4m1 * t4m1); // 10.1*((x2-1)^2 + (x4-1)^2)
            double term6 = 19.8 * t2m1 * t4m1;                 // 19.8*(x2-1)*(x4-1)

            return term1 + term2 + term3 + term4 + term5 + term6;
        }

        @Override
        public double value(RealVector x) {
            double F = computeInnerF(x);
            return 0.5 * F * F;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);
            double x4 = x.getEntry(3);

            // F(x) as in Fortran MODE=3 (FX before the final square)
            double F = computeInnerF(x);

            // From Fortran:
            // GF(1)=FX*(-198*X1 - 2)
            // GF(2)=FX*(20.2*X2 + 19.8*X4 + 60)
            // GF(3)=FX*(-178*X3 - 2)
            // GF(4)=FX*(19.8*X2 + 20.2*X4 + 50)

            double g1 = F * (-198.0 * x1 - 2.0);
            double g2 = F * (20.2 * x2 + 19.8 * x4 + 60.0);
            double g3 = F * (-178.0 * x3 - 2.0);
            double g4 = F * (19.8 * x2 + 20.2 * x4 + 50.0);

            return new ArrayRealVector(new double[]{g1, g2, g3, g4}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Let BFGS build an approximation
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS255() {

        // Initial guess: X = (-3, 1, -3, 1)
        double[] x0 = new double[]{-3.0, 1.0, -3.0, 1.0};

        // Bounds: -10 ≤ xi ≤ 10
        SimpleBounds bounds = new SimpleBounds(
                new double[]{-10.0, -10.0, -10.0, -10.0},
                new double[]{ 10.0,  10.0,  10.0,  10.0}
        );

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
               new InitialGuess(x0),
                new ObjectiveFunction(new HS255Obj()),
                null,   // no equalities
                null   // no inequalities
               // bounds
        );

        double f = sol.getValue();

        double fExpected = 0.0;
        double tol = 1e-4 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
