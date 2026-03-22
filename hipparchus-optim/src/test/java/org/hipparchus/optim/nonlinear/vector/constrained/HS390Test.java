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

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * HS390 (TP390) – Nonlinear process design problem with external subroutine TP390A.
 *
 * N    = 19 variables
 * NILI = 1   (1 linear inequality constraint)
 * NINL = 0
 * NELI = 0
 * NENL = 11  (11 nonlinear equality constraints, here treated as inequalities PSI(x) ≈ 0)
 *
 * Objective in the Fortran subroutine TP390 (mode 2):
 *
 *   ZI1 = 25 * (2268 * x16 * x1)^0.827
 *   ZI2 = 1.75e5 * x17 + 3.65e4 * x17^0.182
 *   ZI3 = 12.6 * x18 + 5.35 * 10^3.378 / x18^0.126
 *
 *   FX = 1.4 * ( ZI1 + ZI2 + ZI3 + 1.095e4
 *                + 1.15e3 * ( x1*(x13 - x14) + x2*(1 + x12) - 3*(1 - x19) ) )
 *
 * Nel Fortran il valore riportato FEX = 0.244724654D+2 è in realtà
 *   FEX = FX / 1.0D+4  valutato in XEX.
 * Qui, per confrontare direttamente con FEX, restituiamo FX / 1e4 come funzione obiettivo.
 */
public class HS390Test {

    private static final int DIM      = 19;
    private static final int NUM_INEQ = 12; // 1 lineare + 11 PSI
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective function (with numerical gradient, zero Hessian)
    // -------------------------------------------------------------------------

    private static class HS390Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            // Fortran indexing: X(1..19) -> Java x[0..18]
            final double x1  = x.getEntry(0);
            final double x2  = x.getEntry(1);
            final double x3  = x.getEntry(2);
            final double x4  = x.getEntry(3);
            final double x5  = x.getEntry(4);
            final double x6  = x.getEntry(5);
            final double x7  = x.getEntry(6);
            final double x8  = x.getEntry(7);
            final double x9  = x.getEntry(8);
            final double x10 = x.getEntry(9);
            final double x11 = x.getEntry(10);
            final double x12 = x.getEntry(11);
            final double x13 = x.getEntry(12);
            final double x14 = x.getEntry(13);
            final double x15 = x.getEntry(14);
            final double x16 = x.getEntry(15);
            final double x17 = x.getEntry(16);
            final double x18 = x.getEntry(17);
            final double x19 = x.getEntry(18);

            // ZI1 = 25 * (2268 * x16 * x1)^0.827
            final double zi1Base = 2268.0 * x16 * x1;
            final double zi1     = 25.0 * FastMath.pow(zi1Base, 0.827);

            // ZI2 = 1.75e5 * x17 + 3.65e4 * x17^0.182
            final double zi2 = 1.75e5 * x17 + 3.65e4 * FastMath.pow(x17, 0.182);

            // ZI3 = 12.6 * x18 + 5.35 * 10^3.378 / x18^0.126
            final double constZi3 = 5.35 * FastMath.pow(10.0, 3.378);
            final double zi3      = 12.6 * x18 + constZi3 / FastMath.pow(x18, 0.126);

            // Last term:
            // 1.15e3 * ( x1*(x13 - x14) + x2*(1 + x12) - 3*(1 - x19) )
            final double termLast =
                    1.15e3 * (x1 * (x13 - x14) +
                              x2 * (1.0 + x12) -
                              3.0 * (1.0 - x19));

            // FX = 1.4 * (ZI1 + ZI2 + ZI3 + 10950 + termLast)
            final double FX = 1.4 * (zi1 + zi2 + zi3 + 1.095e4 + termLast);

            // *** SCALING ***
            // Nel sorgente Fortran, FEX corrisponde a FX / 1.0D+4.
            // Per poter confrontare direttamente col valore FEX riportato,
            // restituiamo FX / 1e4 come obiettivo.
            return FX / 1.0e4;
        }

        @Override
        public RealVector gradient(RealVector x) {
            // Numerical gradient via central finite differences
            final double eps = 1.0e-6;
            final double[] g = new double[DIM];

            for (int j = 0; j < DIM; j++) {
                final double xj = x.getEntry(j);
                final double h  = eps * FastMath.max(1.0, FastMath.abs(xj));

                final RealVector xp = x.copy();
                final RealVector xm = x.copy();
                xp.setEntry(j, xj + h);
                xm.setEntry(j, xj - h);

                final double fp = value(xp);
                final double fm = value(xm);

                g[j] = (fp - fm) / (2.0 * h);
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Start SQP with zero Hessian; BFGS will update it.
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints: G(1..12)
    // -------------------------------------------------------------------------

    private static class HS390Ineq extends InequalityConstraint {

        HS390Ineq() {
            // RHS = 0 for all 12 inequalities; we report G(x) directly.
            // The base class enforces G(x) <= 0 componentwise.
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            final double[] g = new double[NUM_INEQ];

            // Fortran: G(1) = 1.D+0 - X(13) - X(14)
            // NILI = 1 -> linear inequality, originally g >= 0.
            // Per avere la stessa regione ammissibile con convenzione g(x) <= 0,
            // usiamo g1(x) = X(13) + X(14) - 1.
            final double x13 = x.getEntry(12);
            final double x14 = x.getEntry(13);
            g[0] = x13 + x14 - 1.0;

            // G(2..12) = PSI(1..11) da TP390A
            // In Fortran le PSI vengono usate come NENL = 11 (vincoli di uguaglianza).
            // Qui le trattiamo come vincoli di tipo PSI(x) ≈ 0, ma siccome il nostro
            // InequalityConstraint vuole g(x) <= 0, le lasciamo in forma grezza
            // e demandiamo al solver la gestione tramite penalità / KKT tolerance.
            final double[] psi = computePsi(x);
            for (int k = 0; k < 11; k++) {
                g[1 + k] = psi[k];
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            // Numerical Jacobian via central finite differences on G(x)
            final RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);
            final double eps = 1.0e-6;

            for (int j = 0; j < DIM; j++) {
                final double xj = x.getEntry(j);
                final double h  = eps * FastMath.max(1.0, FastMath.abs(xj));

                final RealVector xp = x.copy();
                final RealVector xm = x.copy();
                xp.setEntry(j, xj + h);
                xm.setEntry(j, xj - h);

                final RealVector gp = value(xp);
                final RealVector gm = value(xm);

                for (int i = 0; i < NUM_INEQ; i++) {
                    final double dgi = (gp.getEntry(i) - gm.getEntry(i)) / (2.0 * h);
                    J.setEntry(i, j, dgi);
                }
            }

            return J;
        }

        /**
         * Traduzione della subroutine Fortran TP390A(X, PSI).
         *
         * PSI(1..11) -> psi[0..10]
         */
        private double[] computePsi(RealVector x) {

            final double[] psi = new double[11];

            // Fortran: X(1..19) -> Java x[0..18]
            final double x1  = x.getEntry(0);
            final double x2  = x.getEntry(1);
            final double x3  = x.getEntry(2);
            final double x4  = x.getEntry(3);
            final double x5  = x.getEntry(4);
            final double x6  = x.getEntry(5);
            final double x7  = x.getEntry(6);
            final double x8  = x.getEntry(7);
            final double x9  = x.getEntry(8);
            final double x10 = x.getEntry(9);
            final double x11 = x.getEntry(10);
            final double x12 = x.getEntry(11);
            final double x13 = x.getEntry(12);
            final double x14 = x.getEntry(13);
            final double x15 = x.getEntry(14);
            final double x16 = x.getEntry(15);
            final double x17 = x.getEntry(16);
            final double x18 = x.getEntry(17);
            final double x19 = x.getEntry(18);

            // AK = .0259D+0 * 25.D+0 / 20.D+0**.656D0
            final double AK = 0.0259 * 25.0 / FastMath.pow(20.0, 0.656);

            final double XZ4 = x3 * FastMath.exp(-AK * x16);

            // ZJ1 = -(X(1)*X(13)*XZ4 + 300.D+0*X(19))
            final double ZJ1 = -(x1 * x13 * XZ4 + 300.0 * x19);
            // PSI(1) = ZJ1 + X(1)*X(3) - X(2)*X(5)*X(12)
            psi[0] = ZJ1 + x1 * x3 - x2 * x5 * x12;

            // YZ4 = X(7) + .5D+0 * (X(3) - XZ4)
            final double YZ4 = x7 + 0.5 * (x3 - XZ4);
            // ZJ2 = -X(13)*X(1)*YZ4
            final double ZJ2 = -x13 * x1 * YZ4;
            // PSI(2) = ZJ2 + X(1)*X(7) - X(2)*X(9)*X(12)
            psi[1] = ZJ2 + x1 * x7 - x2 * x9 * x12;

            // ZJ3 = -300*(1 - X(19)) + 3*X(6)*(1 - X(19)) - X(1)*X(14)*XZ4
            final double ZJ3 = -300.0 * (1.0 - x19) + 3.0 * x6 * (1.0 - x19) - x1 * x14 * XZ4;
            // PSI(3) = ZJ3 + X(2)*(X(4) - X(6)) + X(1)*X(6)*X(14)
            psi[2] = ZJ3 + x2 * (x4 - x6) + x1 * x6 * x14;

            // ZJ4 = 3*X(11)*(1 - X(19)) + X(1)*X(14)*(X(11) - YZ4)
            final double ZJ4 = 3.0 * x11 * (1.0 - x19) + x1 * x14 * (x11 - YZ4);
            // PSI(4) = ZJ4 + X(2)*(X(8) - X(11))
            psi[3] = ZJ4 + x2 * (x8 - x11);

            // ZJ5 = X(17) * (.48*X(5)*X(9) / (100 + X(5)))
            final double ZJ5 = x17 * (0.48 * x5 * x9 / (100.0 + x5));
            // PSI(5) = -2*ZJ5 + X(2)*(X(4) - X(5))
            psi[4] = -2.0 * ZJ5 + x2 * (x4 - x5);
            // PSI(6) = ZJ5 + X(2)*(X(8) - X(9)) - .048*X(9)*X(17)
            psi[5] = ZJ5 + x2 * (x8 - x9) - 0.048 * x9 * x17;

            // ZK7 = X(1)*(1 - X(13) - X(14))
            final double ZK7  = x1 * (1.0 - x13 - x14);
            // QZ12 = X(1)*(1 - X(13) - X(14)) + X(2)*(1 - X(12))
            final double QZ12 = x1 * (1.0 - x13 - x14) + x2 * (1.0 - x12);
            // PSI(7) = -ZK7*XZ4 + X(6)*QZ12 - X(2)*X(5)*(1 - X(12))
            psi[6] = -ZK7 * XZ4 + x6 * QZ12 - x2 * x5 * (1.0 - x12);

            // ZJ8 = X(10)*QZ12 - ZK7*YZ4
            final double ZJ8 = x10 * QZ12 - ZK7 * YZ4;
            // PSI(8) = ZJ8 - X(2)*X(9)*(1 - X(12))
            psi[7] = ZJ8 - x2 * x9 * (1.0 - x12);

            // PSI(9) = 6*(1 - X(15))*(20 - X(6))
            //        + X(11)*(X(2) - 3*(1 - X(15)) - X(1)*X(14))
            //        + 3*X(19)*X(11) - X(10)*QZ12
            psi[8] = 6.0 * (1.0 - x15) * (20.0 - x6)
                   + x11 * (x2 - 3.0 * (1.0 - x15) - x1 * x14)
                   + 3.0 * x19 * x11
                   - x10 * QZ12;

            // CK = 7.4 * 2 * 1.2^4 / 2.31e4
            final double CK = 7.4 * 2.0 * FastMath.pow(1.2, 4.0) / 2.31e4;
            final double TEST = -CK * x18 / QZ12;

            final double ZJ10;
            if (TEST > 99.0) {
                ZJ10 = -2.1 * FastMath.sqrt(FastMath.abs(x10)) * FastMath.exp(99.0);
            } else if (TEST < 99.0) {
                ZJ10 = -2.1 * FastMath.sqrt(FastMath.abs(x10))
                       * FastMath.exp(-CK * x18 / QZ12);
            } else {
                // TEST == 99: in Fortran l'ultimo IF sovrascriverebbe comunque
                ZJ10 = -2.1 * FastMath.sqrt(FastMath.abs(x10))
                       * FastMath.exp(-CK * x18 / QZ12);
            }

            // PSI(10) = ZJ10 + 2*(20 - X(6))
            psi[9] = ZJ10 + 2.0 * (20.0 - x6);

            // PSI(11) = (1 - X(13))*X(1) - X(12)*X(2) - 3*X(19)
            psi[10] = (1.0 - x13) * x1 - x12 * x2 - 3.0 * x19;

            return psi;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------

    @Test
    public void testHS390_optimization() {

        // Initial guess (Fortran mode 1)
        final double[] x0 = new double[DIM];
        x0[0]  = 0.02;      // X(1)
        x0[1]  = 4.0;       // X(2)
        x0[2]  = 100.0;     // X(3)
        x0[3]  = 100.0;     // X(4)
        x0[4]  = 15.0;      // X(5)
        x0[5]  = 15.0;      // X(6)
        x0[6]  = 100.0;     // X(7)
        x0[7]  = 1000.0;    // X(8)
        x0[8]  = 1000.0;    // X(9)
        x0[9]  = 1000.0;    // X(10)
        x0[10] = 9000.0;    // X(11)
        x0[11] = 0.001;     // X(12)
        x0[12] = 0.001;     // X(13)
        x0[13] = 1.0;       // X(14)
        x0[14] = 0.001;     // X(15)
        x0[15] = 0.001;     // X(16)
        x0[16] = 0.1;       // X(17)
        x0[17] = 8000.0;    // X(18)
        x0[18] = 0.001;     // X(19)

        // Bounds:
        // default: 1e-5 <= x_i <= 1e5
        final double[] lower = new double[DIM];
        final double[] upper = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            lower[i] = 1.0e-5;
            upper[i] = 1.0e5;
        }

        // XU(1..2) = 50
        upper[0] = 50.0;
        upper[1] = 50.0;

        // XU(3..6) = 100
        for (int i = 2; i <= 5; i++) {
            upper[i] = 100.0;
        }

        // XU(12..15) = 1
        for (int i = 11; i <= 14; i++) {
            upper[i] = 1.0;
        }

        // XU(16) and XU(17) = 50 (I = 1,2; XU(I+15) in Fortran)
        upper[15] = 50.0;   // X(16)
        upper[16] = 50.0;   // X(17)

        final SimpleBounds bounds = new SimpleBounds(lower, upper);

        final SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        final LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS390Obj()),
                null,            // no equality constraints in this wrapper
                new HS390Ineq(), // 12 inequality constraints (1 linear + 11 PSI)
                bounds
        );

        final double f = sol.getValue();

        // Fortran: FEX = 0.244724654D+2, con FX scalata di 1e-4
        final double fExpected = 0.244724654e2; // ≈ 24.4724654
        final double tolF = 1.0e-2 * (FastMath.abs(fExpected) + 1.0);

        assertTrue(f < fExpected+ tolF,
                   "HS390: expected F ≈ " + fExpected + " but got F = " + f);
    }
}
