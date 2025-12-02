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
 * HS372 (TP372) – 9-variable nonlinear least-squares-type problem
 * formulated here as:
 *
 *   f(x) = sum_{i=4..9} x_i^2
 *
 * with 12 nonlinear inequality constraints.
 *
 * From the Fortran subroutine TP372:
 *
 *   N     = 9
 *   NILI  = 0
 *   NINL  = 12
 *   NELI  = 0
 *   NENL  = 0
 *
 * Decision variables:
 *   x1..x9
 *
 * Objective:
 *   FX = sum_{i=4..9} x(i)**2
 *
 * Bounds:
 *   x1, x2: free (no bounds)
 *   x3:     -1 <= x3 <= 0
 *   x4..x9: xj >= 0, no upper bounds
 *
 * Nonlinear inequality constraints (G <= 0 convention in Java):
 *
 *   For convenience define k_i = 2*i - 7 for i = 1..6:
 *       k = (-5, -3, -1, 1, 3, 5)
 *
 *   For i = 1..6:
 *
 *     G_i(x) =
 *       x1 + x2 * exp(k_i * x3) + x_{i+3} - c_i
 *
 *     where
 *       c = (127, 151, 379, 421, 460, 426)
 *
 *   For i = 1..6:
 *
 *     G_{i+6}(x) =
 *       -x1 - x2 * exp(k_i * x3) + x_{i+3} + c_i
 *
 * These correspond 1:1 alle espressioni Fortran:
 *
 *   G(1)  = X(1)+X(2)*EXP(-5*X(3))+X(4) -127
 *   ...
 *   G(6)  = X(1)+X(2)*EXP( 5*X(3))+X(9) -426
 *
 *   G(7)  =-X(1)-X(2)*EXP(-5*X(3))+X(4)+127
 *   ...
 *   G(12) =-X(1)-X(2)*EXP( 5*X(3))+X(9)+426
 *
 * Reference value (LEX = .FALSE. in TP372):
 *
 *   FEX = 0.13390093D+5 = 13390.093
 *
 * In TP372 LEX = .FALSE., quindi FEX è solo un valore di riferimento
 * (upper bound): nel test verifichiamo f <= FEX + tol.
 */
public class HS372Test {

    private static final int DIM      = 9;
    private static final int NUM_INEQ = 12;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS372Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        /**
         * f(x) = sum_{i=4..9} x_i^2
         *
         * In Fortran (MODE = 2):
         *
         *   FX = 0
         *   DO I = 4,9
         *       F(I-3) = X(I)
         *       FX     = FX + X(I)**2
         *   END DO
         */
        @Override
        public double value(RealVector x) {
            double fx = 0.0;
            for (int i = 3; i < DIM; i++) { // indices 3..8 correspond to x4..x9
                double xi = x.getEntry(i);
                fx += xi * xi;
            }
            return fx;
        }

        /**
         * Gradient:
         *
         *   ∂f/∂x_i = 0   for i = 1..3
         *   ∂f/∂x_i = 2*x_i for i = 4..9
         *
         * In Fortran (MODE = 3):
         *
         *   DO I = 4,9
         *       GF(I) = 2*X(I)
         *   END DO
         *   GF(1..3) remain 0.
         */
        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[DIM];

            // x1, x2, x3 do not appear in the objective
            g[0] = 0.0;
            g[1] = 0.0;
            g[2] = 0.0;

            // x4..x9: derivative 2*x_i
            for (int i = 3; i < DIM; i++) {
                double xi = x.getEntry(i);
                g[i] = 2.0 * xi;
            }

            return new ArrayRealVector(g, false);
        }

        /**
         * Hessian:
         *
         * f(x) = sum_{i=4..9} x_i^2
         *
         * ⇒ ∂²f/∂x_i² = 0 for i = 1..3,
         *    ∂²f/∂x_i² = 2 for i = 4..9,
         *    cross-derivatives = 0.
         *
         * In TP372 non viene usata esplicitamente, ma la forniamo per completezza.
         */
        @Override
        public RealMatrix hessian(RealVector x) {
            RealMatrix H = new Array2DRowRealMatrix(DIM, DIM);

            // x1..x3: second derivative 0
            // x4..x9: second derivative 2
            for (int i = 3; i < DIM; i++) {
                H.setEntry(i, i, 2.0);
            }

            return H;
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints (12 nonlinear constraints, G <= 0)
    // -------------------------------------------------------------------------
    private static class HS372Ineq extends InequalityConstraint {

        // Right-hand sides all 0 (G(x) <= 0)
        HS372Ineq() {
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        /**
         * Inequality constraints G(x) <= 0.
         *
         * For i = 1..6, define:
         *   k_i = 2*i - 7  →  [-5, -3, -1, 1, 3, 5]
         *   c_i = [127, 151, 379, 421, 460, 426]
         *
         * Then:
         *
         *   G_i(x)   = x1 + x2 * exp(k_i * x3) + x_{i+3} - c_i
         *   G_{i+6}(x) = -x1 - x2 * exp(k_i * x3) + x_{i+3} + c_i
         *
         * Direct translation of MODE = 4 in TP372.
         */
        @Override
        public RealVector value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);

            final double[] c = {127.0, 151.0, 379.0, 421.0, 460.0, 426.0};
            final int[] k    = {-5, -3, -1, 1, 3, 5};

            double[] g = new double[NUM_INEQ];

            for (int i = 0; i < 6; i++) {
                int ki   = k[i];
                double ci = c[i];
                double expTerm = FastMath.exp(ki * x3);
                double xip3    = x.getEntry(3 + i); // x4..x9

                // G(i+1)   = x1 + x2*exp(k_i*x3) + x_{i+3} - c_i
                g[i] = x1 + x2 * expTerm + xip3 - ci;

                // G(i+7)   = -x1 - x2*exp(k_i*x3) + x_{i+3} + c_i
                g[i + 6] = -x1 - x2 * expTerm + xip3 + ci;
            }

            return new ArrayRealVector(g, false);
        }

        /**
         * Jacobian J (NUM_INEQ x DIM):
         *
         * For i = 1..6 (0-based idx i = 0..5):
         *
         *   G_i(x)   = x1 + x2 * e^{k_i x3} + x_{i+3} - c_i
         *   G_{i+6}(x) = -x1 - x2 * e^{k_i x3} + x_{i+3} + c_i
         *
         * Derivatives:
         *
         *   For G_i:
         *     dG_i/dx1       = 1
         *     dG_i/dx2       = e^{k_i x3}
         *     dG_i/dx3       = x2 * k_i * e^{k_i x3}
         *     dG_i/dx_{i+3}  = 1
         *     others         = 0
         *
         *   For G_{i+6}:
         *     dG_{i+6}/dx1     = -1
         *     dG_{i+6}/dx2     = -e^{k_i x3}
         *     dG_{i+6}/dx3     = -x2 * k_i * e^{k_i x3}
         *     dG_{i+6}/x_{i+3} = 1
         *     others           = 0
         *
         * Questo corrisponde alla combinazione del setup di GG in MODE=1
         * (colonne 1 e 4..9) più gli aggiornamenti in MODE=5 (colonne 2 e 3).
         */
        @Override
        public RealMatrix jacobian(RealVector x) {
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);

            final int[] k = {-5, -3, -1, 1, 3, 5};

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            for (int i = 0; i < 6; i++) {
                int ki       = k[i];
                double expKi = FastMath.exp(ki * x3);
                double commonDx3 = x2 * ki * expKi;
                int rowUpper = i;       // G(i+1)
                int rowLower = i + 6;   // G(i+7)
                int colXip3  = 3 + i;   // x4..x9

                // G_i
                J.setEntry(rowUpper, 0, 1.0);             // dG_i/dx1
                J.setEntry(rowUpper, 1, expKi);           // dG_i/dx2
                J.setEntry(rowUpper, 2, commonDx3);       // dG_i/dx3
                J.setEntry(rowUpper, colXip3, 1.0);       // dG_i/dx_{i+3}

                // G_{i+6}
                J.setEntry(rowLower, 0, -1.0);            // dG_{i+6}/dx1
                J.setEntry(rowLower, 1, -expKi);          // dG_{i+6}/dx2
                J.setEntry(rowLower, 2, -commonDx3);      // dG_{i+6}/dx3
                J.setEntry(rowLower, colXip3, 1.0);       // dG_{i+6}/dx_{i+3}
            }

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS372_optimization() {

        // Initial guess from MODE = 1:
        //
        //   X(1)=  3.0E+1
        //   X(2)= -1.0E+2
        //   X(3)= -0.1997
        //   X(4)= 127
        //   X(5)= 151
        //   X(6)= 379
        //   X(7)= 421
        //   X(8)= 460
        //   X(9)= 426
        //
        double[] x0 = new double[DIM];
        x0[0] = 30.0;
        x0[1] = -100.0;
        x0[2] = -0.1997;
        x0[3] = 127.0;
        x0[4] = 151.0;
        x0[5] = 379.0;
        x0[6] = 421.0;
        x0[7] = 460.0;
        x0[8] = 426.0;

        // Bounds from MODE = 1:
        //
        //   For i = 4..9:
        //      LXL(i) = .TRUE., LXU(i) = .FALSE., XL(i) = 0
        //      ⇒ x_i >= 0, no upper bound
        //
        //   For i = 1,2:
        //      LXL(i) = .FALSE., LXU(i) = .FALSE.
        //      ⇒ no bounds
        //
        //   For i = 3:
        //      XL(3) = -1, XU(3) = 0
        //      LXL(3) = .TRUE., LXU(3) = .TRUE.
        //      ⇒ -1 <= x3 <= 0
        //
        double[] lower = new double[DIM];
        double[] upper = new double[DIM];

        // x1, x2: free
        lower[0] = Double.NEGATIVE_INFINITY;
        upper[0] = Double.POSITIVE_INFINITY;
        lower[1] = Double.NEGATIVE_INFINITY;
        upper[1] = Double.POSITIVE_INFINITY;

        // x3: -1 <= x3 <= 0
        lower[2] = -1.0;
        upper[2] = 0.0;

        // x4..x9: >= 0, no upper bound
        for (int i = 3; i < DIM; i++) {
            lower[i] = 0.0;
            upper[i] = Double.POSITIVE_INFINITY;
        }

        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS372Obj()),
                new HS372Ineq(),  // 12 nonlinear inequality constraints
                bounds
        );

        double f = sol.getValue();

        // From TP372:
        //   LEX = .FALSE.
        //   FEX = 0.13390093D+5 = 13390.093
        //
        // Quindi usiamo FEX come upper bound:
        final double fExpected = 0.13390093e5;
        final double tolF      = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        assertTrue(fExpected + tolF >= f,
                   "HS372: expected F <= " + (fExpected + tolF) + " but got F = " + f);
    }
}
