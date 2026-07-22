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

package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

/** HS288 — 20D unconstrained (sum of squared residuals) translated from TP288. */
public class HS288Test {

    /** f is sum of squares of these residuals (for i=0..4):
     *  f1 = x[i]   + 10*x[i+5]
     *  f2 = sqrt(5)*(x[i+10] - x[i+15])
     *  f3 = (x[i+5] - 2*x[i+10])^2
     *  f4 = sqrt(10)*(x[i] - x[i+15])^2
     *  FEX = 0 at x* = 0 (all zeros)
     */
    static final class HS288Objective extends TwiceDifferentiableFunction {

        private static final double SQRT5  = Math.sqrt(5.0);
        private static final double SQRT10 = Math.sqrt(10.0);

        @Override
        public int dim() { return 20; }

        @Override
        public double value(RealVector x) {
            double fx = 0.0;
            for (int i = 0; i < 5; i++) {
                double xi    = x.getEntry(i);
                double xi5   = x.getEntry(i + 5);
                double xi10  = x.getEntry(i + 10);
                double xi15  = x.getEntry(i + 15);

                double f1 = xi + 10.0 * xi5;
                double f2 = SQRT5 * (xi10 - xi15);
                double f3 = (xi5 - 2.0 * xi10);
                f3 = f3 * f3;
                double f4 = SQRT10 * (xi - xi15);
                f4 = f4 * f4;

                fx += f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4;
            }
            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[20];

            for (int i = 0; i < 5; i++) {
                int i0 = i, i5 = i + 5, i10 = i + 10, i15 = i + 15;

                double xi   = x.getEntry(i0);
                double xi5  = x.getEntry(i5);
                double xi10 = x.getEntry(i10);
                double xi15 = x.getEntry(i15);

                double f1 = xi + 10.0 * xi5;
                double f2 = SQRT5 * (xi10 - xi15);
                double t3 = (xi5 - 2.0 * xi10);         // inner of f3
                double f3 = t3 * t3;
                double t4 = (xi - xi15);               // inner of f4 / sqrt(10)
                double f4 = SQRT10 * t4;
                f4 = f4 * f4;                           // = 10 * (xi - xi15)^2

                // Chain rule: grad = sum_k 2*f_k * df_k/dx
                // df1/dxi  = 1
                // df1/dxi5 = 10
                g[i0]  += 2.0 * f1 * 1.0;
                g[i5]  += 2.0 * f1 * 10.0;

                // df2/dxi10 = sqrt(5), df2/dxi15 = -sqrt(5)
                g[i10] += 2.0 * f2 * SQRT5;
                g[i15] += 2.0 * f2 * (-SQRT5);

                // f3 = (xi5 - 2*xi10)^2 -> df3/dxi5 = 2*(xi5 - 2*xi10)
                //                           df3/dxi10 = -4*(xi5 - 2*xi10)
                double df3_d_xi5  =  2.0 * t3;
                double df3_d_xi10 = -4.0 * t3;
                g[i5]  += 2.0 * f3 * df3_d_xi5;
                g[i10] += 2.0 * f3 * df3_d_xi10;

                // In value():
                //   f4 = (sqrt(10) * (xi - xi15))^2 = 10 * (xi - xi15)^2
                // So:
                //   df4/dxi   =  20 * (xi - xi15)
                //   df4/dxi15 = -20 * (xi - xi15)
                double df4_d_xi   =  20.0 * t4;
                double df4_d_xi15 = -20.0 * t4;
                g[i0]  += 2.0 * f4 * df4_d_xi;
                g[i15] += 2.0 * f4 * df4_d_xi15;
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public org.hipparchus.linear.RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private static LagrangeSolution solve(double[] x0) {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer() ;// obbligatorio per il tracing;
        SQPOption sqpOption=new SQPOption();
        sqpOption.setGradientMode(GradientMode.FORWARD);
        return opt.optimize(
                sqpOption,
                new InitialGuess(x0),
                new ObjectiveFunction(new HS288Objective())
        );
    }

    @Test
    public void testHS288() {
        // Fortran init:
        // x(1..5)= 3, x(6..10)=-1, x(11..15)=0, x(16..20)=10
        double[] x0 = new double[20];
        for (int i = 0; i < 5; i++) {
            x0[i]      = 3.0;
            x0[i + 5]  = -1.0;
            x0[i + 10] = 0.0;
            x0[i + 15] = 10.0;
        }

        LagrangeSolution sol = solve(x0);
        double f = sol.getValue();
        double expected = 0.0; // FEX
        assertEquals(expected, f, 1.0e-2 * (Math.abs(expected) + 1.0),
                     "objective mismatch at optimum");
    }
}