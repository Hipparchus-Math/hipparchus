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
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS369Test {

    private static class HS369Obj extends TwiceDifferentiableFunction {
        @Override
        public int dim() { return 8; }

        @Override
        public double value(RealVector x) {
            return x.getEntry(0) + x.getEntry(1) + x.getEntry(2);
        }

        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { return null; }
    }

    private static class HS369Ineq extends InequalityConstraint {
        // 22 inequality constraints total:
        // - 3 linear
        // - 3 nonlinear
        // - 16 bounds (8 lower + 8 upper)
        HS369Ineq() { super(new ArrayRealVector(22)); }

        private static double safeDiv(double numerator, double denominator) {
            return numerator / Math.max(Math.abs(denominator), 1e-8);
        }

        @Override
        public RealVector value(RealVector x) {
            double[] g = new double[22];
            double[] c = {
                833.33252, 100.0, -83333.333,
                1250.0, 1.0, -1250.0,
                1250000.0, 1.0, -2500.0,
                0.0025, 0.0025, 0.0025, 0.0025,
                -0.0025, 0.01, -0.01
            };

            // Linear inequality constraints
            g[0] = 1.0 - c[9]  * x.getEntry(3) - c[10] * x.getEntry(5);
            g[1] = 1.0 - c[11] * x.getEntry(4) - c[12] * x.getEntry(6) - c[13] * x.getEntry(3);
            g[2] = 1.0 - c[14] * x.getEntry(7) - c[15] * x.getEntry(4);

            // Nonlinear inequality constraints with safe division
            g[3] = 1.0 - safeDiv(c[0] * x.getEntry(3), x.getEntry(0) * x.getEntry(5))
                       - safeDiv(c[1], x.getEntry(5))
                       - safeDiv(c[2], x.getEntry(0) * x.getEntry(5));

            g[4] = 1.0 - safeDiv(c[3] * x.getEntry(4), x.getEntry(1) * x.getEntry(6))
                       - safeDiv(c[4] * x.getEntry(3), x.getEntry(6))
                       - safeDiv(c[5] * x.getEntry(3), x.getEntry(1) * x.getEntry(6));

            g[5] = 1.0 - safeDiv(c[6], x.getEntry(2) * x.getEntry(7))
                       - safeDiv(c[7] * x.getEntry(4), x.getEntry(7))
                       - safeDiv(c[8] * x.getEntry(4), x.getEntry(2) * x.getEntry(7));

            // Bound constraints
            double[] lb = {100.0, 1000.0, 1000.0, 10.0, 10.0, 10.0, 10.0, 10.0};
            double[] ub = {10000.0, 10000.0, 10000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0};
            for (int i = 0; i < 8; i++) {
                g[6 + i]      = x.getEntry(i) - lb[i];     // lower bound: x_i >= lb_i
                g[14 + i]     = ub[i] - x.getEntry(i);     // upper bound: x_i <= ub_i → ub_i - x_i >= 0
            }

            return new ArrayRealVector(g, false);
        }

        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 8; }
    }
    
    private static class HS369IneqNoBounds extends InequalityConstraint {
        // 22 inequality constraints total:
        // - 3 linear
        // - 3 nonlinear
        // - 16 bounds (8 lower + 8 upper)
        HS369IneqNoBounds() { super(new ArrayRealVector(22-16)); }

       

        @Override
        public RealVector value(RealVector x) {
            double[] g = new double[22-16];
            double[] c = {
                833.33252, 100.0, -83333.333,
                1250.0, 1.0, -1250.0,
                1250000.0, 1.0, -2500.0,
                0.0025, 0.0025, 0.0025, 0.0025,
                -0.0025, 0.01, -0.01
            };

            // Linear inequality constraints
            g[0] = 1.0 - c[9]  * x.getEntry(3) - c[10] * x.getEntry(5);
            g[1] = 1.0 - c[11] * x.getEntry(4) - c[12] * x.getEntry(6) - c[13] * x.getEntry(3);
            g[2] = 1.0 - c[14] * x.getEntry(7) - c[15] * x.getEntry(4);

            // Nonlinear inequality constraints with safe division
            g[3] = 1.0 - c[0] * x.getEntry(3)/ (x.getEntry(0) * x.getEntry(5))
                       - c[1]/ x.getEntry(5)
                       - c[2]/( x.getEntry(0) * x.getEntry(5));

            g[4] = 1.0 - c[3] * x.getEntry(4)/ (x.getEntry(1) * x.getEntry(6))
                       - c[4] * x.getEntry(3)/ x.getEntry(6)
                       - c[5] * x.getEntry(3)/( x.getEntry(1) * x.getEntry(6));

            g[5] = 1.0 - c[6]/ (x.getEntry(2) * x.getEntry(7))
                       - c[7] * x.getEntry(4)/ x.getEntry(7)
                       - c[8] * x.getEntry(4)/ (x.getEntry(2) * x.getEntry(7));


            return new ArrayRealVector(g, false);
        }
         @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 8; }
    }

    @Test
    public void testHS369() {
        SQPOption sqpOption = new SQPOption();
        sqpOption.setMaxLineSearchIteration(20);
//        sqpOption.setB(0.5);
//        sqpOption.setMu(1.0e-4);
//        sqpOption.setEps(1e-10);

        double[] start = {5000, 5000, 5000, 200, 350, 150, 225, 425};
        InitialGuess guess = new InitialGuess(start);

        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        double val = 7049.2480;
        LagrangeSolution sol = optimizer.optimize(
            sqpOption,
            guess,
            new ObjectiveFunction(new HS369Obj()),
            new HS369Ineq()
        );

        HSProblemTestUtils.assertExpectedObjective(val, sol); // tolleranza rilassata per robustezza
    }
    
    @Test
    public void testHS369Bound() {
        SQPOption sqpOption = new SQPOption();
        sqpOption.setMaxLineSearchIteration(50);
        sqpOption.setB(0.5);
        sqpOption.setMu(1.0e-4);
        sqpOption.setEps(1e-10);

        double[] start = {5000, 5000, 5000, 200, 350, 150, 225, 425};
        InitialGuess guess = new InitialGuess(start);
        double[] lb = {100.0, 1000.0, 1000.0, 10.0, 10.0, 10.0, 10.0, 10.0};
         double[] ub = {10000.0, 10000.0, 10000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0};
        SimpleBounds bounds=new SimpleBounds(lb,ub);
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        optimizer.setDebugPrinter(s -> {});

        double val = 7049.2480;
        LagrangeSolution sol = optimizer.optimize(
            sqpOption,
            guess,
            new ObjectiveFunction(new HS369Obj()),
            new HS369IneqNoBounds(),bounds
        );

        HSProblemTestUtils.assertExpectedObjective(val, sol);  
    }
}
