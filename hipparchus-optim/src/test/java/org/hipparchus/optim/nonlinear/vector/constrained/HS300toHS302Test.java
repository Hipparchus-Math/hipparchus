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

package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;


public class HS300toHS302Test {

   
    static final class HS300Obj extends TwiceDifferentiableFunction {
        private final int n;
        HS300Obj(int n) { this.n = n; }
        @Override public int dim() { return n; }

        @Override
        public double value(RealVector x) {
            double fx = x.getEntry(0) * x.getEntry(0) - 2.0 * x.getEntry(0);
            for (int i = 1; i < n; i++) {
                double xi   = x.getEntry(i);
                double xim1 = x.getEntry(i - 1);
                fx += 2.0 * xi * xi - 2.0 * xim1 * xi;
            }
            return fx;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[n];
            // i = 0
            g[0] = 2.0 * x.getEntry(0) - 2.0 * x.getEntry(1) - 2.0;
            // i = 1 .. n-2
            for (int i = 1; i < n - 1; i++) {
                g[i] = 4.0 * x.getEntry(i) - 2.0 * x.getEntry(i - 1) - 2.0 * x.getEntry(i + 1);
            }
            // i = n-1
            if (n >= 2) {
                g[n - 1] = 4.0 * x.getEntry(n - 1) - 2.0 * x.getEntry(n - 2);
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public org.hipparchus.linear.RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private static double[] initX(int n) {
        double[] x = new double[n];
        for (int i = 0; i < n; i++) x[i] = 0.0;
        return x;
    }

    private static void runCase(int n, double fExpected) {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(initX(n)),
                new ObjectiveFunction(new HS300Obj(n))
        );

        double f = sol.getValue();
        assertEquals(fExpected, f, 1.0e-6 * (Math.abs(fExpected) + 1.0),
                "objective mismatch for n=" + n);
    }

    @Test public void testHS300() { runCase(20,  -20.0); }
    @Test public void testHS301() { runCase(50,  -50.0); }
    @Test public void testHS302() { runCase(100, -100.0); }
}
