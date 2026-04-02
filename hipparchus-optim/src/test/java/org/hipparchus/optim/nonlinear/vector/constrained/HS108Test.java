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

public class HS108Test {

    private static final int DIM = 9;

    /** TP108 MODE=2/3 objective. */
    private static final class HS108Obj extends TwiceDifferentiableFunction {
        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(final RealVector x) {
            return -0.5 * (x.getEntry(0) * x.getEntry(3)
                         - x.getEntry(1) * x.getEntry(2)
                         + x.getEntry(2) * x.getEntry(8)
                         - x.getEntry(4) * x.getEntry(8)
                         + x.getEntry(4) * x.getEntry(7)
                         - x.getEntry(5) * x.getEntry(6));
        }

        @Override
        public RealVector gradient(final RealVector x) {
            final double[] g = new double[DIM];
            g[0] = -0.5 * x.getEntry(3);
            g[1] =  0.5 * x.getEntry(2);
            g[2] =  0.5 * (x.getEntry(1) - x.getEntry(8));
            g[3] = -0.5 * x.getEntry(0);
            g[4] =  0.5 * (x.getEntry(8) - x.getEntry(7));
            g[5] =  0.5 * x.getEntry(6);
            g[6] =  0.5 * x.getEntry(5);
            g[7] = -0.5 * x.getEntry(4);
            g[8] = -g[7] - g[1];
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(final RealVector x) {
            throw new UnsupportedOperationException();
        }
    }

    /** TP108 MODE=4 inequalities G(1..13) >= 0. */
    private static final class HS108Ineq extends InequalityConstraint {
        HS108Ineq() {
            super(new ArrayRealVector(new double[] {
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            }));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);
            final double x6 = x.getEntry(5);
            final double x7 = x.getEntry(6);
            final double x8 = x.getEntry(7);
            final double x9 = x.getEntry(8);
            return new ArrayRealVector(new double[] {
                    1.0 - x3 * x3 - x4 * x4,
                    1.0 - x9 * x9,
                    1.0 - x5 * x5 - x6 * x6,
                    1.0 - x1 * x1 - (x2 - x9) * (x2 - x9),
                    1.0 - (x1 - x5) * (x1 - x5) - (x2 - x6) * (x2 - x6),
                    1.0 - (x1 - x7) * (x1 - x7) - (x2 - x8) * (x2 - x8),
                    1.0 - (x3 - x5) * (x3 - x5) - (x4 - x6) * (x4 - x6),
                    1.0 - (x3 - x7) * (x3 - x7) - (x4 - x8) * (x4 - x8),
                    1.0 - x7 * x7 - (x8 - x9) * (x8 - x9),
                    x1 * x4 - x2 * x3,
                    x3 * x9,
                    -x5 * x9,
                    x5 * x8 - x6 * x7
            }, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {
            throw new UnsupportedOperationException();
        }
    }

    @Test
    public void testHS108() {
        final InitialGuess guess = new InitialGuess(new double[] {1, 1, 1, 1, 1, 1, 1, 1, 1});
        final double[] lb = new double[DIM];
        final double[] ub = new double[DIM];
        for (int i = 0; i < DIM; i++) {
            lb[i] = 0.0;
            ub[i] = 1.0;
        }
        final SQPOption option = new SQPOption();
        

        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        final LagrangeSolution sol = optimizer.optimize(
                option,
                guess,
                new ObjectiveFunction(new HS108Obj()),
                new HS108Ineq(),
                new SimpleBounds(lb, ub)
        );

        HSProblemTestUtils.assertExpectedObjective(-0.866025403841, sol);
    }
}