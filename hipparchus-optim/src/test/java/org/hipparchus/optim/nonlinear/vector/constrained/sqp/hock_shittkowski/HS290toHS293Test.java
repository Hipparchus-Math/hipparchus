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

/** HS290–HS293: F(x) = ( Σ_{i=1..n} i * x_i^2 )^2, with x* = 0, f* = 0. */
public class HS290toHS293Test {

    static final class WeightedSquareSum extends TwiceDifferentiableFunction {
        private final int n;
        WeightedSquareSum(int n) { this.n = n; }
        @Override public int dim() { return n; }

        @Override
        public double value(RealVector x) {
            double s = 0.0;
            for (int i = 0; i < n; i++) {
                double xi = x.getEntry(i);
                s += (i + 1) * xi * xi;              // F(1) in Fortran
            }
            return s * s;                              // FX = F(1)^2
        }

        @Override
        public RealVector gradient(RealVector x) {
            double s = 0.0;
            for (int i = 0; i < n; i++) {
                double xi = x.getEntry(i);
                s += (i + 1) * xi * xi;               // F(1)
            }
            // DF(1,i) = 2*(i+1)*x_i
            // GF(i)   = 2*F(1)*DF(1,i) = 4*(i+1)*F(1)*x_i
            double[] g = new double[n];
            for (int i = 0; i < n; i++) {
                g[i] = 4.0 * (i + 1) * s * x.getEntry(i);
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public org.hipparchus.linear.RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private static LagrangeSolution solve(int n) {
        double[] x0 = new double[n];
        for (int i = 0; i < n; i++) x0[i] = 10.0;      // Fortran init: X(i)=.1D+1
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer()  ;
         SQPOption sqpOption=new SQPOption();
        sqpOption.setGradientMode(GradientMode.EXTERNAL);// richiesto
        return opt.optimize(
                sqpOption,
                new InitialGuess(x0),
                new ObjectiveFunction(new WeightedSquareSum(n))
        );
    }

    private static void runCase(int n) {
        LagrangeSolution sol = solve(n);
        double f = sol.getValue();
        double expected = 0.0;                         // FEX
        assertEquals(expected, f, 1.0e-3 * (Math.abs(expected) + 1.0),
                     "objective mismatch");
    }

    @Test public void testHS290() { runCase(2);  }
    @Test public void testHS291() { runCase(10); }
    @Test public void testHS292() { runCase(30); }
    @Test public void testHS293() { runCase(50); }
}
