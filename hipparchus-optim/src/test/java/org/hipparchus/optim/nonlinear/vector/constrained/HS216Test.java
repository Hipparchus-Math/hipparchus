/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
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
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

public class HS216Test {

    private static final int DIM = 2;
    private static final int NUM_EQ = 1;

    private static class HS216Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            return 100.0 * FastMath.pow(x1 * x1 - x2, 2) +
                   FastMath.pow(x1 - 1.0, 2);
        }

        @Override
        public RealVector gradient(RealVector x) {
            throw new UnsupportedOperationException();
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException();
        }
    }

    private static class HS216Eq extends EqualityConstraint {

        HS216Eq() {
            super(new ArrayRealVector(new double[NUM_EQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            double g1 = x1 * (x1 - 4.0) - 2.0 * x2 + 12.0;
            return new ArrayRealVector(new double[] { g1 }, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            throw new UnsupportedOperationException();
        }
    }

    @Test
    public void testHS216_optimization() {

        double[] x0 = new double[] { -1.2, 1.0 };

        double[] lower = new double[] { -3.0, -3.0 };
        double[] upper = new double[] { 10.0, 10.0 };
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();
        SQPOption option = HSProblemTestUtils.newCentralDifferenceOption();

        LagrangeSolution sol = opt.optimize(
                
                new InitialGuess(x0),
                new ObjectiveFunction(new HS216Obj()),
                new HS216Eq(),
                bounds
        );

        HSProblemTestUtils.assertExpectedObjective(1.0, sol);
    }
}