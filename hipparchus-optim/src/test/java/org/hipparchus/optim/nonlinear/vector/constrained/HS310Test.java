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
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;


public class HS310Test {

    
    static final class HS310Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }

        @Override
        public double value(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double A  = x1 * x2;
            final double B  = 10.0 - x1;
            final double B2 = B * B;
            final double B4 = B2 * B2;
            final double B5 = B4 * B;
            final double C  = B - x2 * B5;

            return (A * A) * (B2) * (C * C);
        }

        @Override
        public RealVector gradient(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double A  = x1 * x2;
            final double B  = 10.0 - x1;
            final double B2 = B * B;
            final double B4 = B2 * B2;
            final double B5 = B4 * B;
            final double C  = B - x2 * B5;

            // 
            final double common = 2.0 * A * B * C;
            final double g1 = common * (x2 - 1.0 - 5.0 * x2 * B4);
            final double g2 = common * (x1 - B5);

            return new ArrayRealVector(new double[] { g1, g2 }, false);
        }

        @Override
        public RealMatrix hessian(final RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private LagrangeSolution solve() {
        final SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer(); // richiesto;
         SQPOption sqpOption=new SQPOption();
         sqpOption.setGradientMode(GradientMode.CENTRAL);
        // 
        final double[] x0 = { -1.2, 1.0 };

       
        return opt.optimize(
                sqpOption,
                new InitialGuess(x0),
                new ObjectiveFunction(new HS310Obj())
        );
    }

    @Test
    public void testHS310_valueOnly() {
        final double val = 0.0;
        final LagrangeSolution sol = solve();
        final double f = sol.getValue();
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
