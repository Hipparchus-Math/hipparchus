/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS303to305Test {

    
    static final class HS303Obj extends TwiceDifferentiableFunction {
        private final int n;
        HS303Obj(int n) { this.n = n; }
        @Override public int dim() { return n; }

        private static double pom(RealVector x) {
            int n = x.getDimension();
            double p = 0.0;
            for (int i = 1; i <= n; i++) {
                p += 0.5 * i * x.getEntry(i - 1);
            }
            return p;
        }

        @Override public double value(RealVector x) {
            double p = pom(x);
            double fx = p*p + p*p*p*p;
            double sum = 0.0;
            for (int i = 0; i < n; i++) sum += x.getEntry(i) * x.getEntry(i);
            return fx + sum;
        }

        @Override public RealVector gradient(RealVector x) {
            double p = pom(x);
            double factor = (2.0 * p + 4.0 * p*p*p); // d/dPOM (POM^2 + POM^4)
            double[] g = new double[n];
            for (int i = 1; i <= n; i++) {
                double dpi = 0.5 * i;
                g[i - 1] = 2.0 * x.getEntry(i - 1) + factor * dpi; // = 2*x_i + POM*i + 2*i*POM^3
            }
            return new ArrayRealVector(g, false);
        }

        @Override public org.hipparchus.linear.RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private static double[] initX(int n) {
        double[] x = new double[n];
        for (int i = 0; i < n; i++) x[i] = 0.1; // X(I)=.1D+0
        return x;
    }

    private void runCase(int n) {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(initX(n)),
                new ObjectiveFunction(new HS303Obj(n))
        );

        double f = sol.getValue();
        double fExpected = 0.0; // FEX
        assertEquals(fExpected, f, 1.0e-6 * (Math.abs(fExpected) + 1.0),
                "objective mismatch (n=" + n + ")");
    }

    @Test public void testHS303() { runCase(20); }  // TP303
    @Test public void testHS304() { runCase(50); }  // TP304
    @Test public void testHS305() { runCase(100); } // TP305
}
