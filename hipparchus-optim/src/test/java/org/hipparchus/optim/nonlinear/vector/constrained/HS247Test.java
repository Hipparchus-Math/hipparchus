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
 * HS247 (TP247)
 *
 * N    = 3
 * NILI = 0
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * Variables:
 *   x = (x1, x2, x3)
 *
 * Bounds from Fortran (MODE=1):
 *   LXL(1)=TRUE, XL(1)=0.1      → x1 ≥ 0.1
 *   LXL(3)=TRUE, XL(3)=-2.5     → x3 ≥ -2.5
 *   LXU(3)=TRUE, XU(3)= 7.5     → x3 ≤  7.5
 *   x2 is free
 *
 * Objective (MODE=2):
 *
 *   XPI   = 2 * asin(1) = 2π
 *   THETA = (1 / (2*XPI)) * atan(x2 / x1)
 *   if (x1 < 0) then THETA = THETA + 0.5
 *
 *   r = sqrt(x1^2 + x2^2)
 *
 *   f(x) = 100 * ( (x3 - 10*THETA)^2 + (r - 1)^2 ) + x3^2
 *
 * Gradient (MODE=3 in Fortran):
 *
 *   THETA  = (1 / (2*XPI)) * atan(x2 / x1)
 *   DTHETA(1) = -x2 / ( (1 + (x2/x1)^2) * x1^2 )
 *   DTHETA(2) =  1  / ( (1 + (x2/x1)^2) * x1   )
 *   DTHETA(3) =  0
 *   if (x1 < 0) then THETA = THETA + 0.5    (DTHETA unchanged)
 *
 *   Let r = sqrt(x1^2 + x2^2).
 *
 *   GF(1) = 100 * ( 20 * (x3 - 10*THETA) * DTHETA(1)
 *                 +  2 * (r - 1) / r * x1 )
 *
 *   GF(2) = 100 * ( 20 * (x3 - 10*THETA) * DTHETA(2)
 *                 +  2 * (r - 1) / r * x2 )
 *
 *   GF(3) = 100 * ( 2 * (x3 - 10*THETA) ) + 2 * x3
 *
 * Reference solution (MODE=1):
 *   x*  = (1, 0, 0)
 *   f*  = 0
 */
public class HS247Test {

    private static final int DIM = 3;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS247Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            // XPI = 2 * asin(1) = 2π
            double xpi = 2.0 * FastMath.asin(1.0);

            // THETA = 1/(2*XPI) * atan(x2/x1)
            double theta = (1.0 / (2.0 * xpi)) * FastMath.atan(x2 / x1);
            if (x1 < 0.0) {
                theta += 0.5;
            }

            double r = FastMath.sqrt(x1 * x1 + x2 * x2);

            double term1 = x3 - 10.0 * theta;
            double term2 = r - 1.0;

            return 100.0 * (term1 * term1 + term2 * term2) + x3 * x3;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            double x3 = x.getEntry(2);

            // XPI = 2 * asin(1) = 2π
            double xpi = 2.0 * FastMath.asin(1.0);

            // Base THETA and its partial derivatives (before the x1<0 branch)
            double ratio   = x2 / x1;
            double denom   = (1.0 + ratio * ratio) * x1 * x1; // for DTHETA(1)
            double theta   = (1.0 / (2.0 * xpi)) * FastMath.atan(ratio);

            double dtheta1 = -x2 / denom;
            double dtheta2 = 1.0 / ((1.0 + ratio * ratio) * x1);
            double dtheta3 = 0.0;

            // Apply the same shift as Fortran (no change in derivative)
            if (x1 < 0.0) {
                theta += 0.5;
            }

            double r = FastMath.sqrt(x1 * x1 + x2 * x2);

            double term1 = x3 - 10.0 * theta;

            // GF(1), GF(2), GF(3) as in Fortran
            double g1 = 100.0 * (20.0 * term1 * dtheta1
                                +  2.0 * (r - 1.0) / r * x1);

            double g2 = 100.0 * (20.0 * term1 * dtheta2
                                +  2.0 * (r - 1.0) / r * x2);

            double g3 = 100.0 * (2.0 * term1) + 2.0 * x3;

            return new ArrayRealVector(new double[]{g1, g2, g3}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Let BFGS / SQP machinery approximate the Hessian.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS247_optimization() {

        // Initial guess from Fortran (MODE=1):
        //   X(1) = 0.1, X(2) = 0, X(3) = 0
        double[] x0 = new double[]{0.1, 0.0, 0.0};

        // Bounds:
        //   x1 ≥ 0.1
        //   x2 free
        //   -2.5 ≤ x3 ≤ 7.5
        double[] lower = new double[]{
                0.1,                     // x1 lower (XL(1))
                Double.NEGATIVE_INFINITY, // x2 free
                -2.5                     // x3 lower
        };
        double[] upper = new double[]{
                Double.POSITIVE_INFINITY, // x1 no upper
                Double.POSITIVE_INFINITY, // x2 no upper
                7.5                       // x3 upper
        };
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // No equalities, no inequalities; only bounds
        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS247Obj()),
                null,      // no equalities
                null,      // no inequalities
                bounds     // bounds from Fortran XL/XU
        );

        double f = sol.getValue();

        final double fExpected = 0.0; // FEX from Fortran
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
