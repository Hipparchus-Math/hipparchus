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
 * HS239 (TP239)
 *
 * N    = 2
 * NILI = 0
 * NINL = 1  (one nonlinear inequality)
 * NELI = 0
 * NENL = 0
 *
 * Objective and gradient come from TP236239.
 *
 * Fortran constraint G(x) ≥ 0 (MODE=4):
 *   G1(x) = x1 * x2 - 7.D+2  = x1 * x2 - 700
 *
 * Bounds (MODE=1):
 *   LXL(i) = .TRUE., XL(i) = 0  →  x_i >= 0
 *   LXU(i) = .TRUE., XU(1)=75, XU(2)=65  → x1 <= 75, x2 <= 65
 *
 * Reference solution (MODE=1):
 *   x*  = (75, 65)
 *   f*  = -58.903436
 */
public class HS239Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 1;
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective = TP236239 (stessa implementazione usata in HS238)
    // -------------------------------------------------------------------------
    private static class HS239Obj extends TwiceDifferentiableFunction {

        // Coefficienti B(1..20) di TP236239 (Fortran → Java)
        private static final double[] B = new double[] {
                7.5196366677e+01,   // B1
               -3.8112755343e+00,   // B2
                1.2693663450e-01,   // B3
               -2.0567665000e-03,   // B4
                1.0345000000e-05,   // B5
               -6.8306567613e+00,   // B6
                3.0234479300e-02,   // B7
               -1.2813448000e-03,   // B8
                3.5255900000e-05,   // B9
               -2.2660000000e-07,   // B10
                2.5645812530e-01,   // B11
               -3.4604030000e-03,   // B12
                1.3513900000e-05,   // B13
               -2.81064434908e+01,  // B14
               -5.2375000000e-06,   // B15
               -6.3000000000e-09,   // B16
                7.0000000000e-10,   // B17
                3.4054620000e-04,   // B18
               -1.6638000000e-06,   // B19
               -2.8673112392e+00   // B20
        };

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double x1_2 = x1 * x1;
            final double x1_3 = x1_2 * x1;
            final double x1_4 = x1_3 * x1;

            final double x2_2 = x2 * x2;
            final double x2_3 = x2_2 * x2;
            final double x2_4 = x2_3 * x2;

            final double x1x2     = x1 * x2;
            final double x1_2x2   = x1_2 * x2;
            final double x1_3x2   = x1_3 * x2;
            final double x1_4x2   = x1_4 * x2;
            final double x1_2x2_2 = x1_2 * x2_2;
            final double x1_3x2_2 = x1_3 * x2_2;
            final double x1_3x2_3 = x1_3 * x2_3;
            final double x1x2_2   = x1 * x2_2;
            final double x1x2_3   = x1 * x2_3;

            final double expTerm = FastMath.exp(5.0e-4 * x1x2);

            double fx = 0.0;
            fx += B[0];
            fx += B[1] * x1;
            fx += B[2] * x1_2;
            fx += B[3] * x1_3;
            fx += B[4] * x1_4;
            fx += B[5] * x2;
            fx += B[6] * x1x2;
            fx += B[7] * x1_2x2;
            fx += B[8] * x1_3x2;
            fx += B[9] * x1_4x2;
            fx += B[10] * x2_2;
            fx += B[11] * x2_3;
            fx += B[12] * x2_4;
            fx += B[13] * (1.0 / (x2 + 1.0));
            fx += B[14] * x1_2x2_2;
            fx += B[15] * x1_3x2_2;
            fx += B[16] * x1_3x2_3;
            fx += B[17] * x1x2_2;
            fx += B[18] * x1x2_3;
            fx += B[19] * expTerm;

            // Fortran: FX = -FX
            return -fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double x1_2 = x1 * x1;
            final double x1_3 = x1_2 * x1;

            final double x2_2 = x2 * x2;
            final double x2_3 = x2_2 * x2;

            final double expTerm = FastMath.exp(5.0e-4 * x1 * x2);

            // df/dx1 (prima del segno -)
            double dfxdx1 = 0.0;
            dfxdx1 += B[1];
            dfxdx1 += B[2] * 2.0 * x1;
            dfxdx1 += B[3] * 3.0 * x1_2;
            dfxdx1 += B[4] * 4.0 * x1_3;
            dfxdx1 += B[6] * x2;
            dfxdx1 += B[7] * 2.0 * x1 * x2;
            dfxdx1 += B[8] * 3.0 * x1_2 * x2;
            dfxdx1 += B[9] * 4.0 * x1_3 * x2;
            dfxdx1 += B[14] * 2.0 * x1 * x2_2;
            dfxdx1 += B[15] * 3.0 * x1_2 * x2_2;
            dfxdx1 += B[16] * 3.0 * x1_2 * x2_3;
            dfxdx1 += B[17] * x2_2;
            dfxdx1 += B[18] * x2_3;
            dfxdx1 += B[19] * expTerm * (5.0e-4 * x2);

            double g1 = -dfxdx1;

            // df/dx2 (prima del segno -)
            double dfxdx2 = 0.0;
            dfxdx2 += B[5];
            dfxdx2 += B[6] * x1;
            dfxdx2 += B[7] * x1_2;
            dfxdx2 += B[8] * x1_3;
            dfxdx2 += B[9] * x1_3 * x1;     // B10 * x1^4
            dfxdx2 += B[10] * 2.0 * x2;
            dfxdx2 += B[11] * 3.0 * x2_2;
            dfxdx2 += B[12] * 4.0 * x2_3;
            dfxdx2 += B[13] * (-1.0 / ((x2 + 1.0) * (x2 + 1.0)));
            dfxdx2 += B[14] * x1_2 * 2.0 * x2;
            dfxdx2 += B[15] * x1_3 * 2.0 * x2;
            dfxdx2 += B[16] * x1_3 * 3.0 * x2_2;
            dfxdx2 += B[17] * x1 * 2.0 * x2;
            dfxdx2 += B[18] * x1 * 3.0 * x2_2;
            dfxdx2 += B[19] * expTerm * (5.0e-4 * x1);

            double g2 = -dfxdx2;

            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Nessuna Hessiana esplicita in Fortran; lasciamo zero.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraint G(x) >= 0
    // -------------------------------------------------------------------------
    private static class HS239Ineq extends InequalityConstraint {

        HS239Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // Fortran G1 = x1 * x2 - 700  (>= 0)
            double g1 = x1 * x2 - 700.0;

            return new ArrayRealVector(new double[]{g1}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // dG1/dx1 = x2, dG1/dx2 = x1
            J.setEntry(0, 0, x2);
            J.setEntry(0, 1, x1);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS239_optimization() {

        // Initial guess (MODE=1): X(1)=10, X(2)=10
        double[] x0 = new double[]{10.0, 10.0};

        // Bounds: 0 <= x1 <= 75, 0 <= x2 <= 65
        double[] lower = new double[]{0.0, 0.0};
        double[] upper = new double[]{75.0, 65.0};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS239Obj()),
                null,               // no equalities
                new HS239Ineq(),    // 1 inequality
                bounds              // box bounds
        );

        double f = sol.getValue();

        final double fExpected = -58.9034360;
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
