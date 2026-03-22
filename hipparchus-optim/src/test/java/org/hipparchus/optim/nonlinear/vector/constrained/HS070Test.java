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

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

/** HS TP70 (Schittkowski). Objective = sum of squares, 1 nonlinear inequality. */
public class HS070Test {

    // --- Costanti e dati come in TP70 ---------------------------------------

    /** C(1..19): 0.1, 1, 2, ..., 18 */
    private static final double[] C = buildC();

    /** YO(1..19) */
    private static final double[] YO = {
        1.89e-3, 0.1038, 0.268, 0.506, 0.577, 0.604, 0.725, 0.898, 0.947, 0.845,
        0.702, 0.528, 0.385, 0.257, 0.159, 0.0869, 0.0453, 0.01509, 1.89e-3
    };

    /** H3 = 1 / 7.658D0 */
    private static final double H3 = 1.0 / 7.658;

    /** 2*pi “Fortran-like” costante usata nel codice (6.2832D0). */
    private static final double TWO_PI_FORTRAN = 6.2832;

    /** Bounds (MODE=1): XL(i)=1e-5; XU(i)=100; eccetto XU(3)=1. */
    private static final double[] LB = { 1e-5, 1e-5, 1e-5, 1e-5 };
    private static final double[] UB = { 100.0, 100.0, 1.0, 100.0 };

    private static double[] buildC() {
        double[] c = new double[19];
        c[0] = 0.1;
        for (int i = 1; i < 19; i++) {
            c[i] = i; // 1,2,...,18
        }
        return c;
    }

    // --- Objective -----------------------------------------------------------

    /** f(x) = sum_i (YC(i) - YO(i))^2, con protezione LOG come nel Fortran. */
    private static class HS070Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 4; }

        @Override public double value(RealVector X) {
            final double[] x = X.toArray();
            final double x1 = x[0], x2 = x[1], x3 = x[2], x4 = x[3];

            // B = x3 + (1 - x3)*x4
            final double B = x3 + (1.0 - x3) * x4;

            // Variabili intermedie e fattori come in TP70
            final double H1 = x1 - 1.0;
            final double H2 = x2 - 1.0;
            final double H5 = B * H3;
            final double H4 = H5 / x4;
            final double H6 = 12.0 * x1 / (12.0 * x1 + 1.0);
            final double H7 = 12.0 * x2 / (12.0 * x2 + 1.0);

            final double V10 = x2 / TWO_PI_FORTRAN;
            final double V11 = B / x4;
            final double V12 = x1 / TWO_PI_FORTRAN;

            // LOG-guard: se qualunque base di potenza/divisione è negativa o non valida,
            // si usa la penalità SUM (x-5)^2 (MODE=8 in TP70).
            boolean badDomain =
                    (B < 0.0) || (V10 < 0.0) || (V11 < 0.0) || (V12 < 0.0)
                 || (!Double.isFinite(H4)) || (!Double.isFinite(H5));

            if (badDomain) {
                double sum = 0.0;
                for (int i = 0; i < 4; i++) sum += FastMath.pow(x[i] - 5.0, 2);
                return sum;
            }

            // Z1..Z7
            final double Z1 = x3 * FastMath.pow(B, x2);
            final double Z2 = FastMath.sqrt(V10);
            final double Z5 = 1.0 - x3;
            final double Z6 = FastMath.pow(V11, x1);
            final double Z7 = FastMath.sqrt(V12);

            double fx = 0.0;
            for (int i = 0; i < 19; i++) {
                final double ci = C[i];
                // V3 = (C*H3)^H2, V4 = exp(x2*(1 - C*H5))
                final double V3 = FastMath.pow(ci * H3, H2);
                final double V4 = FastMath.exp(x2 * (1.0 - ci * H5));
                // V8 = (C*H3)^H1, V9 = exp(x1*(1 - C*H4))
                final double V8 = FastMath.pow(ci * H3, H1);
                final double V9 = FastMath.exp(x1 * (1.0 - ci * H4));

                final double U1 = Z1 * Z2 * V3 * V4 * H7;
                final double U2 = Z5 * Z6 * Z7 * V8 * V9 * H6;

                final double YC = U1 + U2;
                final double Fi = YC - YO[i];
                fx += Fi * Fi;
            }
            return fx;
        }

        @Override public RealVector gradient(RealVector x) {
            throw new UnsupportedOperationException("Analytical gradient not provided");
        }

        @Override public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Analytical Hessian not provided");
        }
    }

    // --- Vincolo di disuguaglianza: g(x) - 1 <= 0 ----------------------------

    /** g(x) = x3 + (1 - x3)*x4 - 1 <= 0 (coerente con il gradiente TP70). */
     private static class HS070Ineq extends InequalityConstraint {
        HS070Ineq() { super(new ArrayRealVector(new double[] { 0.0 })); }

        @Override public int dim() { return 4; }

        @Override public RealVector value(RealVector X) {
            final double[] x = X.toArray();
            final double g = x[2] + (1.0 - x[2]) * x[3] - 1.0;
            return new ArrayRealVector(new double[]{ g });
        }

        @Override public RealMatrix jacobian(RealVector x) {
            throw new UnsupportedOperationException("Analytical Jacobian not provided");
        }
    }

    // --- Test ----------------------------------------------------------------

    @Test
    public void testHS070() {
        final InitialGuess guess = new InitialGuess(new double[]{ 2.0, 4.0, 0.04, 2.0 });
        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        optimizer.setDebugPrinter(System.out::println);

        final LagrangeSolution sol = optimizer.optimize(
                guess,
                new ObjectiveFunction(new HS070Obj()),
                new HS070Ineq(),
                bounds
        );

        // FEX dal Fortran: 0.749846356143D-02
        final double expected = 0.00749846356143;
        assertEquals(expected, sol.getValue(), 1e-6);
    }
}
