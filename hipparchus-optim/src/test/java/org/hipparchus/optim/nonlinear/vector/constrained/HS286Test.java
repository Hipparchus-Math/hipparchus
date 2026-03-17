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
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

/** TP286 — 20-dimensional unconstrained least-squares. */
public class HS286Test {

    /** Objective: f(x) = sum_i F_i(x)^2 with
     *  F_i   = x_i - 10            for i=1..10
     *  F_{i+10} = 100*(x_i^2 - x_{i+10})  for i=1..10
     */
    static final class TP286Objective extends TwiceDifferentiableFunction {
        @Override public int dim() { return 20; }

        @Override
        public double value(RealVector x) {
            double fx = 0.0;
            // First 10 residuals
            for (int i = 0; i < 10; i++) {
                double Fi = x.getEntry(i) - 1.0;
                fx += Fi * Fi;
            }
            // Next 10 residuals
            for (int i = 0; i < 10; i++) {
                double xi  = x.getEntry(i);
                double xip = x.getEntry(10 + i);
                double Fi2 = 10.0 * (xi * xi - xip);
                fx += Fi2 * Fi2;
            }
            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[20];

            // Precompute residuals needed in gradient
            double[] F1 = new double[10];   // x_i - 10
            double[] F2 = new double[10];   // 100*(x_i^2 - x_{i+10})
            for (int i = 0; i < 10; i++) {
                double xi  = x.getEntry(i);
                double xip = x.getEntry(10 + i);
                F1[i] = xi - 1.0;
                F2[i] = 10.0 * (xi * xi - xip);
            }

            // Grad w.r.t. x_1..x_10
            for (int i = 0; i < 10; i++) {
                // df/dx_i = 2*F1_i*1 + 2*F2_i*(d/dx_i of F2_i)
                // dF2_i/dx_i = 100*(2*xi) = 200*xi
                g[i] = 2.0 * F1[i] + 2.0 * F2[i] * (200.0 * x.getEntry(i));
                // = 2*(x_i-10) + 40000 * x_i * (x_i^2 - x_{i+10})
            }

            // Grad w.r.t. x_11..x_20
            for (int i = 0; i < 10; i++) {
                // dF2_i/dx_{i+10} = -100
                g[10 + i] = 2.0 * F2[i] * (-100.0);
                // = -20000 * (x_i^2 - x_{i+10})
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public org.hipparchus.linear.RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private static LagrangeSolution solve(double[] x0) {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer(); // keep debug output;
         SQPOption sqpOption=new SQPOption();
        //sqpOption.setGradientMode(GradientMode.CENTRAL);
        // keep default settings; no bounds; no constraints
        return opt.optimize(
                sqpOption,
                new InitialGuess(x0),
                new ObjectiveFunction(new TP286Objective())
        );
    }

    @Test
    public void testTP286() {
        // Fortran start: x(1..10) = -1.2 ; x(11..20) = 1
        double[] x0 = new double[20];
        for (int i = 0; i < 10; i++) x0[i] = -1.2;
        for (int i = 10; i < 20; i++) x0[i] = 1.0;

        LagrangeSolution sol = solve(x0);

        double f = sol.getValue();
        double expected = 0.0;
        assertEquals(expected, f, 1.0e-2 * (Math.abs(expected) + 1.0), "objective mismatch");
    }
}
