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
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

public class HS059Test {

    private static class HS059Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return 2;
        }

        @Override
        public double value(final RealVector x) {

            final double x1  = x.getEntry(0);
            final double x2  = x.getEntry(1);
            final double x12 = x1 * x1;
            final double x13 = x12 * x1;
            final double x14 = x13 * x1;
            final double x22 = x2 * x2;
            final double x23 = x22 * x2;
            final double x24 = x23 * x2;

            return -75.196
                   + 3.8112 * x1
                   - 0.12694 * x12
                   + 2.0567e-3 * x13
                   - 1.0345e-5 * x14
                   + 6.8306 * x2
                   - 3.0234e-2 * x1 * x2
                   + 1.28134e-3 * x12 * x2
                   - 3.5256e-5 * x13 * x2
                   + 2.266e-7 * x14 * x2
                   - 0.25645 * x22
                   + 3.4604e-3 * x23
                   - 1.3514e-5 * x24
                   + 28.106 / (x2 + 1.0)
                   + 5.2375e-6 * x12 * x22
                   + 6.3e-8 * x13 * x22
                   - 7.0e-10 * x13 * x23
                   - 3.4054e-4 * x1 * x22
                   + 1.6638e-6 * x1 * x23
                   + 2.8673 * FastMath.exp(5.0e-4 * x1 * x2);
        }

        @Override
        public RealVector gradient(final RealVector x) {
            final double eps = 1.0e-6;
            final double[] g = new double[2];
            for (int j = 0; j < 2; ++j) {
                final double xj = x.getEntry(j);
                final double h = eps * FastMath.max(1.0, FastMath.abs(xj));
                final RealVector xp = x.copy();
                final RealVector xm = x.copy();
                xp.setEntry(j, xj + h);
                xm.setEntry(j, xj - h);
                g[j] = (value(xp) - value(xm)) / (2.0 * h);
            }
            return new ArrayRealVector(g, false);
        }

        @Override
        public RealMatrix hessian(final RealVector x) {
            return new Array2DRowRealMatrix(2, 2);
        }
    }

    private static class HS059Ineq extends InequalityConstraint {

        HS059Ineq() {
            super(new ArrayRealVector(new double[] { 0.0, 0.0, 0.0 }));
        }

        @Override
        public int dim() {
            return 2;
        }

        @Override
        public RealVector value(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            return new ArrayRealVector(new double[] {
                x1 * x2 - 700.0,
                x2 - 8.0e-3 * x1 * x1,
                (x2 - 50.0) * (x2 - 50.0) - 5.0 * (x1 - 55.0)
            }, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final RealMatrix j = new Array2DRowRealMatrix(3, 2);
            j.setEntry(0, 0, x2);
            j.setEntry(0, 1, x1);
            j.setEntry(1, 0, -1.6e-2 * x1);
            j.setEntry(1, 1, 1.0);
            j.setEntry(2, 0, -5.0);
            j.setEntry(2, 1, 2.0 * (x2 - 50.0));
            return j;
        }
    }
//    @Disabled
    @Test
    public void testHS059() {
        final InitialGuess guess = new InitialGuess(new double[] { 37.5, 32.5 });
        final SimpleBounds bounds = new SimpleBounds(new double[] { 0.0, 0.0 }, new double[] { 75.0, 65.0 });
        SQPOption option=new SQPOption();
        option.setGradientMode(GradientMode.CENTRAL);
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        final LagrangeSolution sol = optimizer.optimize(
                option,
                guess,
                new ObjectiveFunction(new HS059Obj()),
                null,
                new HS059Ineq(),
                bounds);

        HSProblemTestUtils.assertExpectedObjective(-7.80422632408, sol);
    }
}