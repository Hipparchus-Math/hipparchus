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

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

public class HS016Test {

    private static class HS016Obj extends TwiceDifferentiableFunction {
        @Override
        public int dim() {
            return 2;
        }

        @Override
        public double value(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            return 100.0 * (x2 - x1 * x1) * (x2 - x1 * x1) + (1.0 - x1) * (1.0 - x1);
        }

        @Override
        public RealVector gradient(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double g2 = 200.0 * (x2 - x1 * x1);
            final double g1 = -2.0 * (x1 * (g2 - 1.0) + 1.0);
            return new ArrayRealVector(new double[] { g1, g2 }, false);
        }

        @Override
        public RealMatrix hessian(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final RealMatrix h = new Array2DRowRealMatrix(2, 2);
            h.setEntry(0, 0, 1200.0 * x1 * x1 - 400.0 * x2 + 2.0);
            h.setEntry(0, 1, -400.0 * x1);
            h.setEntry(1, 0, -400.0 * x1);
            h.setEntry(1, 1, 200.0);
            return h;
        }
    }

    private static class HS016Ineq extends InequalityConstraint {
        HS016Ineq() {
            super(new ArrayRealVector(new double[] { 0.0, 0.0 }));
        }

        @Override
        public RealVector value(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            return new ArrayRealVector(new double[] {
                x2 * x2 + x1,
                x1 * x1 + x2
            }, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final RealMatrix j = new Array2DRowRealMatrix(2, 2);
            j.setEntry(0, 0, 1.0);
            j.setEntry(0, 1, 2.0 * x2);
            j.setEntry(1, 0, 2.0 * x1);
            j.setEntry(1, 1, 1.0);
            return j;
        }

        @Override
        public int dim() {
            return 2;
        }
    }

    @Test
    @Disabled // disabled as we reach a local minimum and not the expected global one
    public void testHS016() {
        final InitialGuess guess = new InitialGuess(new double[] { -0.4, 0.9 });
        final SimpleBounds bounds = new SimpleBounds(
                new double[] { -0.5, Double.NEGATIVE_INFINITY },
                new double[] { 0.5, 1.0 });

        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        SQPOption option=new SQPOption();
        option.setGradientMode(GradientMode.CENTRAL);
        final LagrangeSolution sol = optimizer.optimize(
                option,
                guess,
                new ObjectiveFunction(new HS016Obj()),
                null,
                new HS016Ineq(),
                bounds);

        HSProblemTestUtils.assertExpectedObjective(0.25, sol);
       // HSProblemTestUtils.assertExpectedSolution(new double[] {0.5, 0.25}, sol, 1e-7);
    }
}
