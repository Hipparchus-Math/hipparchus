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
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;

import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

/**
 * HS368 – Unconstrained polynomial problem (8 variables, simple bounds 0 ≤ x ≤ 1).
 *
 * Objective:
 *    S2 = Σ x_i^2
 *    S3 = Σ x_i^3
 *    S4 = Σ x_i^4
 *    f(x) = -S2*S4 + S3^2
 *
 * Gradient:
 *    df/dx_i = -2*x_i*S4 - 4*x_i^3*S2 + 6*x_i^2*S3
 *
 * Bounds:
 *    0 ≤ x_i ≤ 1
 *
 * Reference optimal solution:
 *    XEX = [1, 0.5, 0.5, 1, 1, 1, 0.5, 0.5]
 *    FEX = -0.74997564
 */
public class HS368Test {

    private static final int DIM = 8;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS368Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {

            double S2 = 0.0;
            double S3 = 0.0;
            double S4 = 0.0;

            for (int i = 0; i < DIM; i++) {
                double xi = x.getEntry(i);
                double x2 = xi * xi;
                S2 += x2;
                S3 += x2 * xi;
                S4 += x2 * x2;
            }

            return -S2 * S4 + S3 * S3;
        }

        @Override
        public RealVector gradient(RealVector x) {

            double S2 = 0.0;
            double S3 = 0.0;
            double S4 = 0.0;

            for (int i = 0; i < DIM; i++) {
                double xi = x.getEntry(i);
                double x2 = xi * xi;
                S2 += x2;
                S3 += x2 * xi;
                S4 += x2 * x2;
            }

            double[] g = new double[DIM];

            for (int i = 0; i < DIM; i++) {
                double xi = x.getEntry(i);

                g[i] =
                    -2.0 * xi * S4
                    - 4.0 * FastMath.pow(xi, 3.0) * S2
                    + 6.0 * FastMath.pow(xi, 2.0) * S3;
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Hessian not required
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // TEST
    // -------------------------------------------------------------------------
    @Test
    public void testHS368() {

        // Start values: X(i) = 1 - 1/i
        double[] x0 = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            x0[i] = 1.0 - 1.0 / (i + 1.0);
        }
        x0[6] = 0.7;
        x0[7] = 0.7;

        // Bounds: 0 ≤ x ≤ 1
        double[] lower = new double[DIM];
        double[] upper = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            lower[i] = 0.0;
            upper[i] = 1.0;
        }
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        SQPOption option=new SQPOption();
        option.setGradientMode(GradientMode.FORWARD);
        LagrangeSolution sol = opt.optimize(
            option,    
            new InitialGuess(x0),
            new ObjectiveFunction(new HS368Obj()),
            bounds
        );

        final double f = sol.getValue();

        // Best known solution (LEX = .TRUE.)
        final double fExpected = -0.74997564;

        final double tolF = 1.0e-4 * (FastMath.abs(fExpected) + 1.0);

        // For LEX = TRUE → must match accurately
        assertEquals(fExpected, f, tolF,
                "HS368: objective function mismatch");
    }
}
