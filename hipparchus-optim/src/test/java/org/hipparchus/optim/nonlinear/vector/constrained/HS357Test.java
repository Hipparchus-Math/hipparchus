/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.
 *
 * https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*
 * HS357 (TP357)
 *
 * N    = 4
 * NILI = 0
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * Bounds:
 *   0 <= x1 <= 150
 *   0 <= x2 <=  50
 *   0 <= x3 <= 100
 *   0 <= x4 <= 100
 *
 * Initial guess (MODE=1):
 *   x0 = (136, 0, 74.8, 75.7)
 *
 * Reference solution:
 *   x*  = (136.00762, 0.031371415, 73.59439, 72.187426)
 *   f*  = 0.35845660
 *
 * Objective (MODE=2):
 *   Complex geometric functional based on 36 measurement points (XPT, YPT)
 *   and parameters P0,Q0,R0,S0. The code below reproduces Fortran exactly,
 *   including the "failure" branch FX = 1.0e20.
 */

/*
 * HS357 (TP357)
 *
 * N    = 4
 * NILI = 0
 * NINL = 0
 * NELI = 0
 * NENL = 0
 *
 * Bounds:
 *   0 <= x1 <= 150
 *   0 <= x2 <=  50
 *   0 <= x3 <= 100
 *   0 <= x4 <= 100
 *
 * Initial guess (MODE=1):
 *   x0 = (136, 0, 74.8, 75.7)
 *
 * Reference solution:
 *   x*  = (136.00762, 0.031371415, 73.59439, 72.187426)
 *   f*  = 0.35845660
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

public class HS357Test {

    private static final int DIM = 4;

    private static class HS357Obj extends TwiceDifferentiableFunction {

        // XPT(1..36)
        private static final double[] XPT = {
                0.113e3, 0.1101e3, 0.1062e3, 0.1013e3,
                0.954e2, 0.888e2,  0.816e2,  0.74e2,  0.661e2,
                0.584e2, 0.51e2,   0.443e2,  0.387e2, 0.345e2,
                0.324e2, 0.329e2,  0.364e2,  0.428e2, 0.509e2,
                0.59e2,  0.658e2,  0.715e2,  0.765e2, 0.811e2,
                0.856e2, 0.902e2,  0.946e2,  0.989e2, 0.103e3,
                0.1067e3,0.1099e3, 0.1125e3, 0.1144e3,
                0.1155e3,0.1157e3, 0.1149e3
        };

        // YPT(1..36)
        private static final double[] YPT = {
                0.402e2, 0.468e2, 0.533e2, 0.594e2,
                0.65e2,  0.699e2, 0.739e2, 0.769e2, 0.789e2,
                0.798e2, 0.797e2, 0.785e2, 0.765e2, 0.736e2,
                0.702e2, 0.66e2,  0.609e2, 0.543e2, 0.458e2,
                0.361e2, 0.265e2, 0.181e2, 0.114e2, 0.62e1,
                0.26e1,  0.3e0,  -0.7e0,  -0.6e0,   0.7e0,
                0.31e1,  0.64e1,  0.105e2, 0.155e2, 0.21e2,
                0.271e2, 0.336e2
        };

        // P0,Q0,R0,S0 from DATA
        private static final double P0 = 0.9e2;   // 90.0
        private static final double Q0 = 0.0;
        private static final double R0 = 0.0;
        private static final double S0 = 0.0;

        // DALPHA = 0.3141527D+1 / 0.18D+2
        private static final double DALPHA = 0.3141527e1 / 0.18e2;

        @Override
        public int dim() {
            return DIM;
        }

        /** Fortran TP357(MODE=2) translated 1:1 on a raw array X. */
        private double f(final double[] X) {

            final double P1 = X[0];
            final double Q1 = X[1];
            final double R1 = X[2];
            final double S1 = X[3];

            double sum = 0.0;

            // DO 54 I=2,36
            for (int i = 2; i <= 36; i++) {

                double alpha = DALPHA * (i - 1);
                double ca = FastMath.cos(alpha);
                double sa = FastMath.sin(alpha);

                double PI = P1 * ca - Q1 * sa + P0 * (1.0 - ca) + Q0 * sa;
                double QI = P1 * sa + Q1 * ca + Q0 * (1.0 - ca) - P0 * sa;

                double A = R0 * S1 - S0 * R1 - Q1 * R0 + P1 * S0
                         + PI * Q1 - P1 * QI + QI * R1 - PI * S1;

                double B = -R0 * R1 - S0 * S1 + P1 * R0 + Q1 * S0
                         - P1 * PI - Q1 * QI + PI * R1 + QI * S1;

                double C = -R1 * R0 - S1 * S0 + PI * R0 + QI * S0
                         + P1 * R1 + Q1 * S1
                         - (P1 * P1 + Q1 * Q1 + PI * PI + QI * QI) / 2.0;

                double AABB = A * A + B * B;
                double AJ = 1.0;
                double TEST;
                double PH;

                if (AABB < 1.0e-30) {
                    // IF (AABB.LT.0.1D-29)
                    if (FastMath.abs(A) < 1.0e-29) {
                        A = 1.0e-29;
                    }
                    // label 50: PH=-DATAN(B/A)
                    TEST = 0.0;
                    PH = -FastMath.atan(B / A);
                } else {
                    // TEST=C/DSQRT(AABB)
                    TEST = C / FastMath.sqrt(AABB);
                    // IF(DABS(TEST).GT.0.1D+1) GOTO 51
                    if (FastMath.abs(TEST) > 10.0) {
                        // label 51: FX=0.1D+21
                        return 1.0e20;
                    }
                    // 52 PH=DASIN(TEST)-DATAN(B/A)
                    PH = FastMath.asin(TEST) - FastMath.atan(B / A);
                }

                // 55–53–52 loop (at most due tentativi)
                while (true) {

                    double SP = FastMath.sin(PH);
                    double CP = FastMath.cos(PH);

                    double RI = R1 * CP - S1 * SP + PI - P1 * CP + Q1 * SP;
                    double SI = R1 * SP + S1 * CP + QI - P1 * SP - Q1 * CP;

                    double test1 = (R1 - R0) * (R1 - R0)
                                 + (S1 - S0) * (S1 - S0);
                    // IF(TEST1.LT.0.1D-9) TEST1=0.1D-9
                    if (test1 < 1.0e-10) {
                        test1 = 1.0e-10;
                    }

                    double num = test1
                               - (RI - R0) * (RI - R0)
                               - (SI - S0) * (SI - S0);

                    // IF(DABS((TEST1-(RI-R0)**2-(SI-S0)**2)/TEST1).LT.0.1D-2)GOTO53
                    if (FastMath.abs(num / test1) < 1.0e-2) {
                        // 53 CALCX / CALCY
                        double CALCX = XPT[0] * CP - YPT[0] * SP
                                     + PI - P1 * CP + Q1 * SP;
                        double CALCY = XPT[0] * SP + YPT[0] * CP
                                     + QI - P1 * SP - Q1 * CP;

                        int idx = i - 1; // Fortran XPT(I) → XPT[i-1]
                        sum += (CALCX - XPT[idx]) * (CALCX - XPT[idx])
                             + (CALCY - YPT[idx]) * (CALCY - YPT[idx]);
                        break; // esce dal while, va al prossimo I
                    }

                    // IF(AJ.EQ.0.2D+1) GOTO 51
                    if (AJ == 2.0) {
                        // label 51: errore
                        return 1.0e20;
                    }

                    // TEST=-TEST; AJ=0.2D+1; GOTO 52
                    TEST = -TEST;
                    AJ = 2.0;
                    PH = FastMath.asin(TEST) - FastMath.atan(B / A);
                }
            }

            // SQL=(R1-R0)**2+(S1-S0)**2+(R1-P1)**2+(S1-Q1)**2+(P1-P0)**2+(Q1-Q0)**2
            double sql = (R1 - R0) * (R1 - R0)
                       + (S1 - S0) * (S1 - S0)
                       + (R1 - P1) * (R1 - P1)
                       + (S1 - Q1) * (S1 - Q1)
                       + (P1 - P0) * (P1 - P0)
                       + (Q1 - Q0) * (Q1 - Q0);

            // FX=SUM/0.1D+3+SQL/0.625D+5
            return sum / 100.0 + sql / 62500.0;
        }

        @Override
        public double value(RealVector x) {
            return f(x.toArray());
        }

        @Override
        public RealVector gradient(RealVector x) {
            // Gradient numerico (Fortran usa differenze finite in NLPQL1)
            double[] X = x.toArray();
            double[] g = new double[DIM];

            double f0 = f(X);

            for (int k = 0; k < DIM; k++) {
                double xk = X[k];
                double h = 1e-6 * FastMath.max(1.0, FastMath.abs(xk));

                X[k] = xk + h;
                double fp = f(X);

                X[k] = xk - h;
                double fm = f(X);

                X[k] = xk;

                g[k] = (fp - fm) / (2.0 * h);
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // lascia a BFGS
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS357() {

        // X(1)=0.136D+3, X(2)=0, X(3)=0.748D+2, X(4)=0.757D+2
        double[] x0 = new double[]{136.0, 0.0, 74.8, 75.7};

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        // Bounds from MODE=1:
        // 0 <= x1 <= 150; 0 <= x2 <= 50; 0 <= x3 <= 100; 0 <= x4 <= 100
        SimpleBounds bounds = new SimpleBounds(
                new double[]{0.0, 0.0, 0.0, 0.0},
                new double[]{150.0, 50.0, 100.0, 100.0}
        );

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS357Obj()),
                null,   // no equalities
                null,   // no inequalities
                bounds
        );

        double f = sol.getValue();

        double fExpected = 0.35845660;
        double tol = 1e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
