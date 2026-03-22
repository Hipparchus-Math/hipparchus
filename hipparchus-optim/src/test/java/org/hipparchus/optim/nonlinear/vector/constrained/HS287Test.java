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

/** HS287 — 20-dimensional unconstrained problem. */
public class HS287Test {

    /**
     * Objective function:
     *
     * f(x) = 1.0e-5 * sum_{i=1..5} [
     *     1.0e3*(x_i^2 - x_{i+5})^2 + (x_i - 10)^2
     *   + 0.9e2*(x_{i+10}^2 - x_{i+15})^2 + (x_{i+10} - 10)^2
     *   + 1.01e2*((x_{i+5}-10)^2 + (x_{i+15}-10)^2)
     *   + 1.98e2*(x_{i+5}-10)*(x_{i+15}-10)
     * ]
     *
     * Expected optimum: x_i = 10 for all i, f* = 0.
     */
    static final class HS287Objective extends TwiceDifferentiableFunction {

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

                fx += 1.0e3 * Math.pow(xi * xi - xi5, 2)
                    + Math.pow(xi - 10.0, 2)
                    + 0.9e2 * Math.pow(xi10 * xi10 - xi15, 2)
                    + Math.pow(xi10 - 10.0, 2)
                    + 1.01e2 * (Math.pow(xi5 - 10.0, 2) + Math.pow(xi15 - 10.0, 2))
                    + 1.98e2 * (xi5 - 10.0) * (xi15 - 10.0);
            }
            return fx * 1.0e-5;
        }

        @Override
        public RealVector gradient(RealVector x) {
            double[] g = new double[20];

            for (int i = 0; i < 5; i++) {
                double xi    = x.getEntry(i);
                double xi5   = x.getEntry(i + 5);
                double xi10  = x.getEntry(i + 10);
                double xi15  = x.getEntry(i + 15);

                // First group (x_i)
                g[i] = (4.0e3 * (xi * xi - xi5) * xi + 2.0 * (xi - 10.0)) * 1.0e-5;

                // Second group (x_{i+5})
                g[i + 5] = (-2.0e3 * (xi * xi - xi5)
                            + 2.02e2 * (xi5 - 10.0)
                            + 1.98e2 * (xi15 - 10.0)) * 1.0e-5;

                // Third group (x_{i+10})
                g[i + 10] = (3.6e3 * xi10 * (xi10 * xi10 - xi15)
                             + 2.0 * (xi10 - 10.0)) * 1.0e-5;

                // Fourth group (x_{i+15})
                g[i + 15] = (-1.8e3 * (xi10 * xi10 - xi15)
                             + 2.02e2 * (xi15 - 10.0)
                             + 1.98e2 * (xi5 - 10.0)) * 1.0e-5;
            }

            return new ArrayRealVector(g, false);
        }

        @Override
        public org.hipparchus.linear.RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private static LagrangeSolution solve(double[] x0) {
        SQPOptimizerS2 opt = new SQPOptimizerS2();
        opt.setDebugPrinter(System.out::println);
        SQPOption sqpOption=new SQPOption();
        sqpOption.setGradientMode(GradientMode.EXTERNAL);
        return opt.optimize(
                sqpOption,
                new InitialGuess(x0),
                new ObjectiveFunction(new HS287Objective())
        );
    }

    @Test
    public void testHS287() {
        // Fortran start: groupwise initialization
        double[] x0 = new double[20];
        for (int i = 0; i < 5; i++) {
            x0[i] = -3.0;      // x(1..5)
            x0[i + 5] = -1.0;  // x(6..10)
            x0[i + 10] = -3.0; // x(11..15)
            x0[i + 15] = -1.0; // x(16..20)
        }

        LagrangeSolution sol = solve(x0);
        double f = sol.getValue();
        double expected = 0.0;
        assertEquals(expected, f, 1.0e-2 * (Math.abs(expected) + 1.0), "objective mismatch at optimum");
    }
}
