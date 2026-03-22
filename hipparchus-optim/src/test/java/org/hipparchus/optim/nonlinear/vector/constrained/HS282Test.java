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
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/** HS282 (10-dim): least-squares of chained residuals with two anchor terms. */
public class HS282Test {

    /** Objective f(x) = sum_{i=1..9} [H_i (x_i^2 - x_{i+1})]^2 + (x_1-1)^2 + (x_{10}-1)^2.
     *  Gradient computed as 2 * Jᵀ * F with sparse Jacobian.
     */
    static final class HS282Objective extends TwiceDifferentiableFunction {
        @Override public int dim() { return 10; }

        private static double H(int i) {
            // i = 0..8 corresponds to Fortran i=1..9: H = sqrt((11 - i)*11)
            return Math.sqrt((10 - i) * 11.0);
        }

        @Override public double value(RealVector x) {
            // Build residuals F[0..10]
            double[] F = new double[11];

            // F[0..8] = H_i * (x_i^2 - x_{i+1})
            for (int i = 0; i < 9; i++) {
                double hi = H(i);
                F[i] = hi * (x.getEntry(i) * x.getEntry(i) - x.getEntry(i + 1));
            }
            // F[9] = x1 - 1 ; F[10] = x10 - 1
            F[9]  = x.getEntry(0) - 1.0;
            F[10] = x.getEntry(9) - 1.0;

            double fx = 0.0;
            for (double v : F) fx += v * v;
            return fx;
        }

        @Override public RealVector gradient(RealVector x) {
            double[] F = new double[11];
            for (int i = 0; i < 9; i++) {
                double hi = H(i);
                F[i] = hi * (x.getEntry(i) * x.getEntry(i) - x.getEntry(i + 1));
            }
            F[9]  = x.getEntry(0) - 1.0;
            F[10] = x.getEntry(9) - 1.0;

            double[] g = new double[10];

            // k = 0
            {
                double h0 = H(0);
                double dF0_dx0 = 2.0 * h0 * x.getEntry(0);
                double dF9_dx0 = 1.0;
                g[0] = 2.0 * (dF0_dx0 * F[0] + dF9_dx0 * F[9]);
            }

            // k = 1..8
            for (int k = 1; k <= 8; k++) {
                double hm1 = H(k - 1);
                double hk  = H(k);
                double dFkm1_dxk = -hm1;                 // ∂F_{k-1}/∂x_k
                double dFk_dxk   =  2.0 * hk * x.getEntry(k); // ∂F_k/∂x_k
                g[k] = 2.0 * (dFkm1_dxk * F[k - 1] + dFk_dxk * F[k]);
            }

            // k = 9
            {
                double h8 = H(8);
                double dF8_dx9  = -h8;
                double dF10_dx9 = 1.0;
                g[9] = 2.0 * (dF8_dx9 * F[8] + dF10_dx9 * F[10]);
            }

            return new ArrayRealVector(g, false);
        }

        @Override public RealMatrix hessian(RealVector x) {
            // Not required by the solver; provide zeros.
            return new Array2DRowRealMatrix(10, 10);
        }
    }

    private LagrangeSolution solve() {
        final SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        optimizer.setDebugPrinter(System.out::println);

        final double[] start = new double[10];
        start[0] = -1.2; // as in Fortran; others 0

        final double[] lo = new double[10];
        final double[] up = new double[10];
        for (int i = 0; i < 10; i++) {
            lo[i] = Double.NEGATIVE_INFINITY;
            up[i] = Double.POSITIVE_INFINITY;
        }

        return optimizer.optimize(
            new InitialGuess(start),
            new ObjectiveFunction(new HS282Objective()),
            new SimpleBounds(lo, up)
        );
    }

    @Test
    public void testHS282() {
        LagrangeSolution sol = solve();
        double f = sol.getValue();
        double fEx = 0.0;
        assertEquals(fEx, f, 1.0e-6 * (Math.abs(fEx) + 1.0), "objective mismatch at optimum");
    }
}

