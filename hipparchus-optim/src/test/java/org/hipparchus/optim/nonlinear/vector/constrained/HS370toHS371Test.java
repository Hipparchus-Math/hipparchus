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
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Disabled;

/**
 * HS370 / HS371 – Least-squares polynomial recursion problems.
 *
 * TP370: N = 6
 * TP371: N = 9  (same objective structure, different dimension and expected FEX)
 *
 * Objective:
 *   For a given dimension N (6 or 9):
 *
 *   F(1) = x1
 *   F(2) = x2 - x1^2 - 1
 *   For i = 2..30:
 *       basis = (i - 1) / 29
 *       sum1  = Σ_{j=2..N} x_j * (j-1) * basis^(j-2)
 *       sum   = Σ_{j=1..N} x_j * basis^(j-1)
 *       F(i+1) = sum1 - sum^2 - 1
 *
 *   FX = Σ_{k=1..31} F(k)^2
 *
 * No constraints, pure unconstrained least-squares minimization.
 */
public class HS370toHS371Test {

    /**
     * Shared objective for HS370/HS371, parametrized by dimension N.
     */
    private static class HS370371Obj extends TwiceDifferentiableFunction {

        private final int dim;

        HS370371Obj(int dim) {
            this.dim = dim;
        }

        @Override
        public int dim() {
            return dim;
        }

        @Override
        public double value(RealVector x) {

            final double[] xv = x.toArray();

            // F(1..31) stored as F[0..30]
            final double[] fVals = new double[31];

            // F(1) = x1
            fVals[0] = xv[0];

            // F(2) = x2 - x1^2 - 1
            fVals[1] = xv[1] - xv[0] * xv[0] - 1.0;

            // For i = 2..30:
            for (int i = 2; i <= 30; i++) {
                final double basis = (double) (i - 1) / 29.0;

                double sum1 = 0.0;
                // sum1 = Σ_{j=2..N} x_j * (j-1) * basis^(j-2)
                for (int j = 2; j <= dim; j++) {
                    final double xj = xv[j - 1];
                    sum1 += xj * (j - 1) * FastMath.pow(basis, j - 2);
                }

                double sum = 0.0;
                // sum = Σ_{j=1..N} x_j * basis^(j-1)
                for (int j = 1; j <= dim; j++) {
                    final double xj = xv[j - 1];
                    sum += xj * FastMath.pow(basis, j - 1);
                }

                // F(i+1) => index i in zero-based array
                fVals[i] = sum1 - sum * sum - 1.0;
            }

            double fx = 0.0;
            for (int k = 0; k < 31; k++) {
                fx += fVals[k] * fVals[k];
            }

            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            // Analytic gradient is lengthy; rely on optimizer finite differences if available.
            throw new UnsupportedOperationException("HS370/HS371: analytic gradient not implemented");
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Not used; return empty matrix.
            return new Array2DRowRealMatrix(dim, dim);
        }
    }

    /**
     * TP370 – N=6, FEX ≈ 2.28767005355e-3, LEX = .FALSE ⇒ FEX >= f.
     */
    @Test
    public void testHS370() {

        final int n = 6;

        // Initial guess X(I) = 0
        double[] x0 = new double[n];

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
            new InitialGuess(x0),
            new ObjectiveFunction(new HS370371Obj(n))
        );

        final double f = sol.getValue();
        final double val = 0.228767005355e-2; // FEX from TP370
       
       HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
    /**
     * TP371 – N=9, FEX ≈ 1.3997601e-6, LEX = .FALSE ⇒ FEX >= f.
     */
    //@Disabled
    @Test
    public void testHS371() {

        final int n = 9;

        // Initial guess X(I) = 0
        double[] x0 = new double[n];

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }
        SQPOption option=new SQPOption();
        option.setGradientMode(GradientMode.CENTRAL);
        LagrangeSolution sol = opt.optimize(
            option,    
            new InitialGuess(x0),
            new ObjectiveFunction(new HS370371Obj(n))
        );

        final double f = sol.getValue();
        final double val = 0.13997601e-5; // corrected FEX for TP371
        

         HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
