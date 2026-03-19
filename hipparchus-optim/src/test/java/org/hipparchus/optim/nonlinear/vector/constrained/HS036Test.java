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
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.jupiter.api.Test;

public class HS036Test {

    private static class HS036Obj extends TwiceDifferentiableFunction {
        @Override
        public int dim() {
            return 3;
        }

        @Override
        public double value(final RealVector x) {
            // TP36 MODE=2: FX = -X(1)*X(2)*X(3)
            return -x.getEntry(0) * x.getEntry(1) * x.getEntry(2);
        }

        @Override
        public RealVector gradient(final RealVector x) {
            throw new UnsupportedOperationException();
        }

        @Override
        public RealMatrix hessian(final RealVector x) {
            throw new UnsupportedOperationException();
        }
    }

    private static class HS036Ineq extends InequalityConstraint {
        HS036Ineq() {
            // TP36 has only 1 nonlinear inequality: 72 - x1 - 2*x2 - 2*x3 >= 0
            super(new ArrayRealVector(new double[] {0.0}));
        }

        @Override
        public RealVector value(final RealVector x) {
            final double g1 = 72.0 - x.getEntry(0) - 2.0 * x.getEntry(1) - 2.0 * x.getEntry(2);
            return new ArrayRealVector(new double[] {g1});
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {
            throw new UnsupportedOperationException();
        }

        @Override
        public int dim() {
            return 3;
        }
    }

    @Test
    public void testHS036() {
        final InitialGuess guess = new InitialGuess(new double[] {10.0, 10.0, 10.0});

        // TP36 MODE=1 bounds: 0 <= xi, x1<=20, x2<=11, x3<=42
        final SimpleBounds bounds = new SimpleBounds(
                new double[] {0.0, 0.0, 0.0},
                new double[] {20.0, 11.0, 42.0}
        );

        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        final LagrangeSolution sol = optimizer.optimize(
                guess,
                new ObjectiveFunction(new HS036Obj()),
                new HS036Ineq(),
                bounds
        );

        HSProblemTestUtils.assertExpectedObjective(-3300.0, sol);
    }
}