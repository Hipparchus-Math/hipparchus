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
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * HS105 (TP105) – Mixture of three normal distributions, maximum likelihood.
 *
 * N    = 8 variables
 * NILI = 1 (1 linear inequality constraint)
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * Fortran objective:
 *
 *   if all mixture terms are "safe":
 *       V = 1 / sqrt(8 * atan(1))  = 1 / sqrt(2π)
 *       V1 = x1 / x6
 *       V2 = x2 / x7
 *       V3 = (1 - x1 - x2) / x8
 *       V4 = 1 / (2 * x6^2)
 *       V5 = 1 / (2 * x7^2)
 *       V6 = 1 / (2 * x8^2)
 *
 *       For i = 1..235
 *         A_i = V1 * exp(max(-(Y_i - x3)^2 * V4, -10))
 *         B_i = V2 * exp(max(-(Y_i - x4)^2 * V5, -10))
 *         C_i = V3 * exp(max(-(Y_i - x5)^2 * V6, -10))
 *         V11 = (A_i + B_i + C_i) * V
 *         if V11 <= 1e-5 → fallback branch
 *         S += log(V11)
 *
 *       FX = -S
 *
 *   else (fallback branch – label 70):
 *       SUM = Σ_{j=1..8} (x_j - 5)^2
 *       FX  = SUM + 2.09e3
 *
 * Constraint:
 *
 *   G1(x) = 1 - x1 - x2  <= 0 (in Fortran: G(1) = 1 - X(1) - X(2))
 *
 * Bounds:
 *   X1, X2:  1e-3 <= X <= 0.499
 *   X3:      100  <= X <= 180
 *   X4:      130  <= X <= 210
 *   X5:      170  <= X <= 240
 *   X6..X8:  5    <= X <= 25
 *
 * Initial guess:
 *   X = (0.1, 0.2, 100, 125, 175, 11.2, 13.2, 15.8)
 *
 * Reference:
 *   FEX = 0.113841623960D+04
 */
public class HS105Test {

    private static final int DIM      = 8;
    private static final int NUM_INEQ = 1;
    private static final int NUM_EQ   = 0;

    /** Y(1..235) from COMMON /D105/. */
    private static final double[] Y = new double[235];

    static {
        // Direct translation of the Fortran Y initialization:

        // Y(1)=95
        Y[0] = 95.0;
        // Y(2)=105
        Y[1] = 105.0;

        // DO 30 I=3,6       Y(I)=110
        for (int i = 3; i <= 6; i++) {
            Y[i - 1] = 110.0;
        }
        // DO 31 I=7,10      Y(I)=115
        for (int i = 7; i <= 10; i++) {
            Y[i - 1] = 115.0;
        }
        // DO 32 I=11,25     Y(I)=120
        for (int i = 11; i <= 25; i++) {
            Y[i - 1] = 120.0;
        }
        // DO 33 I=26,40     Y(I)=125
        for (int i = 26; i <= 40; i++) {
            Y[i - 1] = 125.0;
        }
        // DO 34 I=41,55     Y(I)=130
        for (int i = 41; i <= 55; i++) {
            Y[i - 1] = 130.0;
        }
        // DO 35 I=56,68     Y(I)=135
        for (int i = 56; i <= 68; i++) {
            Y[i - 1] = 135.0;
        }
        // DO 36 I=69,89     Y(I)=140
        for (int i = 69; i <= 89; i++) {
            Y[i - 1] = 140.0;
        }
        // DO 37 I=90,101    Y(I)=145
        for (int i = 90; i <= 101; i++) {
            Y[i - 1] = 145.0;
        }
        // DO 38 I=102,118   Y(I)=150
        for (int i = 102; i <= 118; i++) {
            Y[i - 1] = 150.0;
        }
        // DO 39 I=119,122   Y(I)=155
        for (int i = 119; i <= 122; i++) {
            Y[i - 1] = 155.0;
        }
        // DO 40 I=123,142   Y(I)=160
        for (int i = 123; i <= 142; i++) {
            Y[i - 1] = 160.0;
        }
        // DO 41 I=143,150   Y(I)=165
        for (int i = 143; i <= 150; i++) {
            Y[i - 1] = 165.0;
        }
        // DO 42 I=151,167   Y(I)=170
        for (int i = 151; i <= 167; i++) {
            Y[i - 1] = 170.0;
        }
        // DO 43 I=168,175   Y(I)=175
        for (int i = 168; i <= 175; i++) {
            Y[i - 1] = 175.0;
        }
        // DO 44 I=176,181   Y(I)=180
        for (int i = 176; i <= 181; i++) {
            Y[i - 1] = 180.0;
        }
        // DO 45 I=182,187   Y(I)=185
        for (int i = 182; i <= 187; i++) {
            Y[i - 1] = 185.0;
        }
        // DO 46 I=188,194   Y(I)=190
        for (int i = 188; i <= 194; i++) {
            Y[i - 1] = 190.0;
        }
        // DO 47 I=195,198   Y(I)=195
        for (int i = 195; i <= 198; i++) {
            Y[i - 1] = 195.0;
        }
        // DO 48 I=199,201   Y(I)=200
        for (int i = 199; i <= 201; i++) {
            Y[i - 1] = 200.0;
        }
        // DO 49 I=202,204   Y(I)=205
        for (int i = 202; i <= 204; i++) {
            Y[i - 1] = 205.0;
        }
        // DO 50 I=205,212   Y(I)=210
        for (int i = 205; i <= 212; i++) {
            Y[i - 1] = 210.0;
        }
        // Y(213)=215
        Y[212] = 215.0;
        // DO 51 I=214,219   Y(I)=220
        for (int i = 214; i <= 219; i++) {
            Y[i - 1] = 220.0;
        }
        // DO 52 I=220,224   Y(I)=230
        for (int i = 220; i <= 224; i++) {
            Y[i - 1] = 230.0;
        }
        // Y(225)=235
        Y[224] = 235.0;
        // DO 53 I=226,232   Y(I)=240
        for (int i = 226; i <= 232; i++) {
            Y[i - 1] = 240.0;
        }
        // Y(233)=245, Y(234)=260, Y(235)=260
        Y[232] = 245.0;
        Y[233] = 260.0;
        Y[234] = 260.0;
    }

    // -------------------------------------------------------------------------
    // Objective function (analytic gradient, zero Hessian)
    // -------------------------------------------------------------------------

    private static class HS105Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);
            final double x6 = x.getEntry(5);
            final double x7 = x.getEntry(6);
            final double x8 = x.getEntry(7);

            // Fallback branch if any mixture term is too small
            double sumPenalty = 0.0;
            boolean underflow = false;

            // Constants
            final double V = 1.0 / FastMath.sqrt(8.0 * FastMath.atan(1.0)); // 1/sqrt(2π)

            final double V1 = x1 / x6;
            final double V2 = x2 / x7;
            final double V3 = (1.0 - x1 - x2) / x8;
            final double V4 = 1.0 / (2.0 * x6 * x6);
            final double V5 = 1.0 / (2.0 * x7 * x7);
            final double V6 = 1.0 / (2.0 * x8 * x8);

            double S = 0.0;

            for (int i = 0; i < 235; i++) {
                final double yi = Y[i];

                final double e1 = FastMath.max(-(yi - x3) * (yi - x3) * V4, -10.0);
                final double e2 = FastMath.max(-(yi - x4) * (yi - x4) * V5, -10.0);
                final double e3 = FastMath.max(-(yi - x5) * (yi - x5) * V6, -10.0);

                final double Ai = V1 * FastMath.exp(e1);
                final double Bi = V2 * FastMath.exp(e2);
                final double Ci = V3 * FastMath.exp(e3);

                final double V11 = (Ai + Bi + Ci) * V;

                if (V11 <= 1.0e-5) {
                    underflow = true;
                    break;
                }

                S += FastMath.log(V11);
            }

            if (underflow) {
                // label 70 in Fortran
                for (int i = 0; i < DIM; i++) {
                    final double diff = x.getEntry(i) - 5.0;
                    sumPenalty += diff * diff;
                }
                return sumPenalty + 2.09e3;
            }

            return -S;
        }

        @Override
        public RealVector gradient(RealVector x) {
            // Gradient is only meaningful in the main branch (no underflow).
            // We implement the analytic gradient exactly as in the Fortran code.

            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);
            final double x6 = x.getEntry(5);
            final double x7 = x.getEntry(6);
            final double x8 = x.getEntry(7);

            double[] g = new double[DIM];

            // First, compute A,B,C for all i as in MODE=2 (without clipping in grad),
            // and DA, DB, DC as in the Fortran MODE=3 block.
            final double[][] DA = new double[235][DIM];
            final double[][] DB = new double[235][DIM];
            final double[][] DC = new double[235][DIM];

            // Precompute A,B,C with *unclipped* exponentials for gradient (as in Fortran)
            final double V0x6 = x6 * x6;
            final double V0x7 = x7 * x7;
            final double V0x8 = x8 * x8;

            final double[] A = new double[235];
            final double[] B = new double[235];
            final double[] C = new double[235];

            for (int i = 0; i < 235; i++) {

                // --- A part ---
                double v2 = Y[i] - x3;
                double v1 = FastMath.exp(-v2 * v2 / (2.0 * V0x6));
                A[i] = x1 * v1 / x6;

                DA[i][0] = v1 / x6;                          // dA/dx1
                DA[i][2] = x1 * v2 / FastMath.pow(x6, 3) * v1;   // dA/dx3
                DA[i][5] = x1 / V0x6 * (v2 * v2 / V0x6 - 1.0) * v1; // dA/dx6

                // --- B part ---
                double v3 = Y[i] - x4;
                double v5 = FastMath.exp(-v3 * v3 / (2.0 * V0x7));
                B[i] = x2 * v5 / x7;

                DB[i][1] = v5 / x7;                          // dB/dx2
                DB[i][3] = x2 * v3 / FastMath.pow(x7, 3) * v5;   // dB/dx4
                DB[i][6] = x2 / V0x7 * (v3 * v3 / V0x7 - 1.0) * v5; // dB/dx7

                // --- C part ---
                double v9 = Y[i] - x5;
                double v8 = FastMath.exp(-v9 * v9 / (2.0 * V0x8));
                double v10 = 1.0 - x1 - x2;
                C[i] = v10 * v8 / x8;

                DC[i][0] = -v8 / x8;                         // dC/dx1
                DC[i][1] = -v8 / x8;                         // dC/dx2
                DC[i][4] = v10 * v9 / FastMath.pow(x8, 3) * v8;   // dC/dx5
                DC[i][7] = v10 / V0x8 * (v9 * v9 / V0x8 - 1.0) * v8; // dC/dx8
            }

            // Now accumulate gradient: GF(J) = - Σ (DA+DB+DC)/(A+B+C)
            for (int j = 0; j < DIM; j++) {
                double t1 = 0.0;
                for (int i = 0; i < 235; i++) {
                    final double denom = A[i] + B[i] + C[i];
                    // If denom is extremely small, we are in the underflow regime;
                    // in that case the Fortran code jumps out before reaching here.
                    if (denom > 0.0) {
                        t1 += (DA[i][j] + DB[i][j] + DC[i][j]) / denom;
                    }
                }
                g[j] = -t1;
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Start with zero Hessian; SQP will update via BFGS.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraint: G1(x) = 1 - x1 - x2 <= 0  (Fortran: NILI=1)
    // -------------------------------------------------------------------------

    private static class HS105Ineq extends InequalityConstraint {

        HS105Ineq() {
            // RHS = 0, we provide g(x) directly.
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {
            double[] g = new double[NUM_INEQ];
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            g[0] = 1.0 - x1 - x2;
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);
            // dG1/dx1 = -1, dG1/dx2 = -1, others zero
            J.setEntry(0, 0, -1.0);
            J.setEntry(0, 1, -1.0);
            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------

    @Test
    public void testHS105_optimization() {

        // Initial guess (Fortran mode 1)
        double[] x0 = new double[DIM];
        x0[0] = 0.1;    // X(1)
        x0[1] = 0.2;    // X(2)
        x0[2] = 100.0;  // X(3)
        x0[3] = 131.0;  // X(4)
        x0[4] = 175.0;  // X(5)
        x0[5] = 11.2;   // X(6)
        x0[6] = 13.2;   // X(7)
        x0[7] = 15.8;   // X(8)

        // Bounds from Fortran
        double[] lower = new double[DIM];
        double[] upper = new double[DIM];

        // X1, X2
        lower[0] = 1.0e-3;
        lower[1] = 1.0e-3;
        upper[0] = 0.499;
        upper[1] = 0.499;

        // X3..X5
        lower[2] = 100.0;
        upper[2] = 180.0;
        lower[3] = 130.0;
        upper[3] = 210.0;
        lower[4] = 170.0;
        upper[4] = 240.0;

        // X6..X8
        for (int i = 5; i < DIM; i++) {
            lower[i] = 5.0;
            upper[i] = 25.0;
        }

        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                //new InitialGuess(x0),
                new ObjectiveFunction(new HS105Obj()),
                null,               // no equality constraints
                new HS105Ineq(),    // 1 inequality constraint
                bounds
        );

        double f = sol.getValue();
        double fExpected = 0.113841623960e4;
        double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        
        assertTrue(f < fExpected+ tol,
                   "HS390: expected F ≈ " + fExpected + " but got F = " + f);
    }
}
