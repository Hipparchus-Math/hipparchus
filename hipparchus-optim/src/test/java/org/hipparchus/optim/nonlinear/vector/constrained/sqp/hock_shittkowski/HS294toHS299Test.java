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

public class HS294toHS299Test {

   
    static final class HS294Obj extends TwiceDifferentiableFunction {
        private final int n;
        HS294Obj(int n) { this.n = n; }
        @Override public int dim() { return n; }

        @Override
        public double value(RealVector x) {
            double f = 0.0;
            for (int i = 0; i < n - 1; i++) {
                double xi = x.getEntry(i);
                double xi1 = x.getEntry(i + 1);
                double t1 = 100.0 * (xi1 - xi * xi);
                double t2 = 1.0 - xi;
                f += t1 * t1 + t2 * t2;
            }
            return f * 1.0e-4;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[n];
            int k = n - 1;
            double[] F = new double[2 * k];
            double[][] DF = new double[2 * k][n];

           
            for (int i = 0; i < k; i++) {
                F[i]     = 100.0 * (x.getEntry(i + 1) - x.getEntry(i) * x.getEntry(i));
                F[i + k] = 10.0 * (1.0 - x.getEntry(i));

                for (int j = 0; j < n; j++) {
                    DF[i][j] = 0.0;
                    DF[i + k][j] = 0.0;
                }
                DF[i][i]     = -200.0 * x.getEntry(i);
                DF[i][i + 1] = 100.0;
                DF[i + k][i] = -10.0;
            }

            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int i = 0; i < 2 * k; i++) {
                    sum += 2.0 * F[i] * DF[i][j];
                }
                g[j] = sum * 1.0e-4;
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
        for (int i = 0; i < n; i++) x[i] = -1.2;   
        for (int i = 0; i < n ; i += 2) x[i] = 1.0; 
        return x;
    }

    private static void runCase(int n) {
        double[] x0 = initX(n);
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS294Obj(n))
        );

        double f = sol.getValue();
        double expected = 0.0;
        assertEquals(expected, f, 1.0e-2 * (Math.abs(expected) + 1.0),
                     "objective mismatch for n=" + n);
    }

    @Test public void testHS294() { runCase(6);   }
    @Test public void testHS295() { runCase(10);  }
    @Test public void testHS296() { runCase(16);  }
    @Test public void testHS297() { runCase(30);  }
    @Test public void testHS298() { runCase(50);  }
    @Test public void testHS299() { runCase(100); }
}
